#include "DefaultPgMatcher.h"

#include "../pseudogenome/DefaultPseudoGenome.h"
#include "../pseudogenome/PackedPseudoGenome.h"
#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../pseudogenome/readslist/SeparatedExtendedReadsList.h"

namespace PgTools {

    DefaultPgMatcher::DefaultPgMatcher(const string& srcPgPrefix, const string& targetPgPrefix, bool revComplMatching)
            :srcPgPrefix(srcPgPrefix), targetPgPrefix(targetPgPrefix), revComplMatching(revComplMatching) {
        destPgIsSrcPg = srcPgPrefix == targetPgPrefix;
        if (destPgIsSrcPg)
            cout << "Reading pseudogenome..." << endl;
        else
            cout << "Reading source pseudogenome..." << endl;
        PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(srcPgPrefix, srcPgh, srcRsProp, plainTextReadMode);
        srcPg = PgTools::SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(srcPgPrefix);
        cout << "Pseudogenome length: " << srcPgh->getPseudoGenomeLength() << endl;
        readLength = srcPgh->getMaxReadLength();

        if (targetPgPrefix != srcPgPrefix) {
            cout << "Reading target pseudogenome..." << endl;
            PseudoGenomeHeader *pgh = 0;
            ReadsSetProperties* rsProp = 0;
            bool plainTextReadMode = false;
            PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(targetPgPrefix, pgh, rsProp,
                                                                                 plainTextReadMode);
            destPg = PgTools::SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(targetPgPrefix);
            cout << "Pseudogenome length: " << pgh->getPseudoGenomeLength() << endl;
            delete(pgh);
            delete(rsProp);
        } else
            destPg = srcPg;

        if (revComplMatching)
            destPg = PgSAHelpers::reverseComplement(destPg);
    }

    DefaultPgMatcher::~DefaultPgMatcher() {
        delete(srcPgh);
        delete(srcRsProp);
    }

    void DefaultPgMatcher::exactMatchPg(uint32_t minMatchLength) {
        TextMatcher* matcher;
        matcher = new DefaultTextMatcher(srcPg, minMatchLength);
        matcher->matchTexts(textMatches, destPg, destPgIsSrcPg, revComplMatching, minMatchLength);
        delete(matcher);

        if (revComplMatching)
            correctDestPositionDueToRevComplMatching();
    }

    void DefaultPgMatcher::writeMatchesInfo(const string &dumpFilePrefix) {
        cout << "Sorry, unimplemented matches info dump feature." << endl;
    }

    using namespace PgTools;

    void DefaultPgMatcher::transferMatchedReads(const string &destPgFilePrefix) {
        srcRl = ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(srcPgPrefix,
                                                                                       srcPgh->getPseudoGenomeLength());

        fillPgMatches();
        mapPgMatches2SrcReadsList();
        if (destPgIsSrcPg)
            reverseDestWithSrcForBetterMatchesMappingInTheSamePg();
        resolveMatchesOverlapInSrc();
        if (destPgIsSrcPg)
            resolveDestSrcReadsOverlapConflictsInTheSamePg();

        transferMatchesFromSrcToDest(destPgFilePrefix);

        delete(srcRl);
    }

    void DefaultPgMatcher::mapPgMatches2SrcReadsList() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.mapping.posSrcText < pgMatch2.mapping.posSrcText; });

        uint32_t cancelledMatchesCount = 0;
        uint_pg_len_max totalNetMatchCharsWithOverlapCount = 0;
        int64_t i = 0;
        for(PgMatch& pgMatch: pgMatches) {
            pgMatch.alignToReads(srcRl->pos, readLength, i);
            if (pgMatch.inactive(srcRl->pos, readLength))
                cancelledMatchesCount++;
            totalNetMatchCharsWithOverlapCount += pgMatch.netSrcLength(srcRl->pos, readLength);
        }
        cout << "Cancelled " << cancelledMatchesCount << " short matches due to external reads overlap." << endl;
        cout << "Stage 2: Net-matched " << getTotalMatchStat(totalNetMatchCharsWithOverlapCount) << " Pg src->dest sum length of all matches." << endl;
    }

    void DefaultPgMatcher::reverseDestWithSrcForBetterMatchesMappingInTheSamePg() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.mapping.posDestText < pgMatch2.mapping.posDestText; });

        uint_pg_len_max totalNetReverseMatchCharsWithOverlapCount = 0;
        uint_pg_len_max totalNetMatchCharsWithOverlapCount = 0;
        uint32_t restoredMatchesCount = 0;
        int64_t i = 0;
        for (PgMatch& pgMatch: pgMatches) {
            TextMatch revTextMatch(pgMatch.mapping.posDestText, pgMatch.mapping.length, pgMatch.mapping.posSrcText);
            PgMatch revMatch(revTextMatch);
            revMatch.alignToReads(srcRl->pos, readLength, i);
            uint_pg_len_max revNetSrcLength = revMatch.netSrcLength(srcRl->pos, readLength);
            totalNetReverseMatchCharsWithOverlapCount += revNetSrcLength;
            if (revNetSrcLength > pgMatch.netSrcLength(srcRl->pos, readLength)) {
                if (pgMatch.netSrcLength(srcRl->pos, readLength) == 0)
                    restoredMatchesCount++;
                pgMatch = revMatch;
            }
            totalNetMatchCharsWithOverlapCount += pgMatch.netSrcLength(srcRl->pos, readLength);
        }

        cout << "Restored " << restoredMatchesCount << " short matches (cancelled earlier due to external reads overlap)." << endl;
        cout << "Stage 2a: Net-matched " << getTotalMatchStat(totalNetReverseMatchCharsWithOverlapCount) << " Pg dest->src sum length of all matches." << endl;
        cout << "Stage 2b: Optimized net-matched " << getTotalMatchStat(totalNetMatchCharsWithOverlapCount) << " Pg sum length of all matches." << endl;
    }

    void DefaultPgMatcher::resolveDestSrcReadsOverlapConflictsInTheSamePg() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.netEndPosSrcPg(srcRl->pos) < pgMatch2.netEndPosSrcPg(srcRl->pos); });

        vector<uint32_t> matchDestOrderedIdx(pgMatches.size());
        for(uint32_t i = 0; i < pgMatches.size(); i++)
            matchDestOrderedIdx[i] = i;
        sort(matchDestOrderedIdx.begin(), matchDestOrderedIdx.end(),
                [this](const uint32_t& pgMatchIdx1, const uint32_t& pgMatchIdx2) -> bool
            { return pgMatches[pgMatchIdx1].netPosDestPg(srcRl->pos, readLength, revComplMatching)
                < pgMatches[pgMatchIdx2].netPosDestPg(srcRl->pos,readLength, revComplMatching); });

        vector<uint_pg_len_max> maxNetEndPosDestPgUpTo(pgMatches.size());
        maxNetEndPosDestPgUpTo[0] = pgMatches[matchDestOrderedIdx[0]].netEndPosDestPg(srcRl->pos, readLength,
                revComplMatching);
        for(uint32_t i = 1; i < pgMatches.size(); i++) {
            maxNetEndPosDestPgUpTo[i] = pgMatches[matchDestOrderedIdx[i]].netEndPosDestPg(srcRl->pos, readLength,
                    revComplMatching);
            if (maxNetEndPosDestPgUpTo[i] < maxNetEndPosDestPgUpTo[i - 1])
                maxNetEndPosDestPgUpTo[i] = maxNetEndPosDestPgUpTo[i - 1];
        }

        uint_pg_len_max totalNetMatchCharsCount = 0;
        uint_pg_len_max collidedCharsCount = 0;
        uint32_t conflictsCount = 0;
        uint32_t resolvedConflictsCount = 0;
        int64_t dOrdIdx = 0;
        for (PgMatch& srcMatch: pgMatches) {
            totalNetMatchCharsCount += srcMatch.netSrcLength(srcRl->pos, readLength);
            while(srcMatch.netPosSrcPg(srcRl->pos, readLength) < maxNetEndPosDestPgUpTo[dOrdIdx] && --dOrdIdx > 0);
            while (++dOrdIdx < pgMatches.size()
                && pgMatches[matchDestOrderedIdx[dOrdIdx]].netPosDestPg(srcRl->pos, readLength, revComplMatching)
                    < srcMatch.netEndPosSrcPg(srcRl->pos)) {
                if (srcMatch.inactive(srcRl->pos, readLength))
                    break;
                PgMatch& destMatch = pgMatches[matchDestOrderedIdx[dOrdIdx]];
                if (!destMatch.inactive(srcRl->pos, readLength)
                    && destMatch.netEndPosDestPg(srcRl->pos, readLength, revComplMatching) > srcMatch.netPosSrcPg(srcRl->pos, readLength)) {
                    if (resolveDestSrcCollision(destMatch, srcMatch, collidedCharsCount))
                        resolvedConflictsCount++;
                    else
                        conflictsCount++;
                }
            }
        }
        totalNetMatchCharsCount -= collidedCharsCount;
        cout << conflictsCount << " src<->dest overlap conflicts result in cancelling matches (resolved " << resolvedConflictsCount << ")." << endl;
        cout << "Stage 3a: src<->dest conflicts filtered net-matched " << getTotalMatchStat(totalNetMatchCharsCount) << " Pg sum length of all matches." << endl;
    }

    bool DefaultPgMatcher::resolveDestSrcCollision(PgMatch &destMatch, PgMatch &srcMatch,
                                                   uint_pg_len_max &collidedCharsCount) {
        if (destMatch.netEndPosDestPg(srcRl->pos, readLength, revComplMatching) < srcMatch.netEndPosSrcPg(srcRl->pos) &&
            destMatch.netPosDestPg(srcRl->pos, readLength, revComplMatching) < srcMatch.netPosSrcPg(srcRl->pos, readLength)) {
            uint_pg_len_max overlapLength = destMatch.netEndPosDestPg(srcRl->pos, readLength, revComplMatching)
                    - srcMatch.netPosSrcPg(srcRl->pos, readLength);
            srcMatch.trimLeft(overlapLength, srcRl->pos, readLength);
            collidedCharsCount += overlapLength;
            return !srcMatch.inactive(srcRl->pos, readLength);
        }
        if (destMatch.netEndPosDestPg(srcRl->pos, readLength, revComplMatching) > srcMatch.netEndPosSrcPg(srcRl->pos) &&
            destMatch.netPosDestPg(srcRl->pos, readLength, revComplMatching) > srcMatch.netPosSrcPg(srcRl->pos, readLength)) {
            uint_pg_len_max overlapLength = srcMatch.netEndPosSrcPg(srcRl->pos)
                    - destMatch.netPosDestPg(srcRl->pos, readLength, revComplMatching);
            srcMatch.trimRight(overlapLength, srcRl->pos, readLength);
            collidedCharsCount += overlapLength;
            return !srcMatch.inactive(srcRl->pos, readLength);
        }
        if (destMatch.netSrcLength(srcRl->pos, readLength) > srcMatch.netSrcLength(srcRl->pos, readLength)) {
            collidedCharsCount += srcMatch.netSrcLength(srcRl->pos, readLength);
            srcMatch.endRlIdx = srcMatch.startRlIdx - 1;
        } else {
            if (&destMatch <= &srcMatch)
                collidedCharsCount += destMatch.netSrcLength(srcRl->pos, readLength);
            destMatch.endRlIdx = destMatch.startRlIdx - 1;
        }
        return false;
    }

    void DefaultPgMatcher::resolveMatchesOverlapInSrc() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.netPosSrcPg(srcRl->pos, readLength) < pgMatch2.netPosSrcPg(srcRl->pos, readLength); });

        vector<uint_pg_len_max> maxNetEndPosSrcPgUpTo(pgMatches.size());
        maxNetEndPosSrcPgUpTo[0] = pgMatches[0].netEndPosSrcPg(srcRl->pos);
        for(uint32_t i = 1; i < pgMatches.size(); i++) {
            maxNetEndPosSrcPgUpTo[i] = pgMatches[i].netEndPosSrcPg(srcRl->pos);
            if (maxNetEndPosSrcPgUpTo[i] < maxNetEndPosSrcPgUpTo[i - 1])
                maxNetEndPosSrcPgUpTo[i] = maxNetEndPosSrcPgUpTo[i - 1];
        }

        uint_pg_len_max totalNetMatchCharsCount = 0;
        uint_pg_len_max collidedCharsCount = 0;
        uint32_t conflictsCount = 0;
        uint32_t resolvedConflictsCount = 0;
        int64_t lIdx = 0;
        for (PgMatch& rightMatch: pgMatches) {
            totalNetMatchCharsCount += rightMatch.netSrcLength(srcRl->pos, readLength);
            while(rightMatch.netPosSrcPg(srcRl->pos, readLength) < maxNetEndPosSrcPgUpTo[lIdx] && --lIdx > 0);
            while (++lIdx < pgMatches.size()
                   && pgMatches[lIdx].netPosSrcPg(srcRl->pos, readLength) < rightMatch.netEndPosSrcPg(srcRl->pos)) {
                if (rightMatch.inactive(srcRl->pos, readLength) || &pgMatches[lIdx] == &rightMatch)
                    break;
                if (!pgMatches[lIdx].inactive(srcRl->pos, readLength)
                    && pgMatches[lIdx].netEndPosSrcPg(srcRl->pos) > rightMatch.netPosSrcPg(srcRl->pos, readLength)) {
                    if (resolveSrcSrcCollision(pgMatches[lIdx], rightMatch, collidedCharsCount))
                        resolvedConflictsCount++;
                    else
                        conflictsCount++;
                }
            }
        }
        totalNetMatchCharsCount -= collidedCharsCount;
        cout << conflictsCount << " src<->src overlap conflicts result in cancelling matches (resolved " << resolvedConflictsCount << ")." << endl;
        cout << "Stage 3: All conflicts filtered net-matched " << getTotalMatchStat(totalNetMatchCharsCount) << " Pg sum length of all matches." << endl;
    }

    bool DefaultPgMatcher::resolveSrcSrcCollision(PgMatch &leftMatch, PgMatch &rightMatch,
                                                  uint_pg_len_max &collidedCharsCount) {
        if (leftMatch.netEndPosSrcPg(srcRl->pos) < rightMatch.netEndPosSrcPg(srcRl->pos) &&
            leftMatch.netPosSrcPg(srcRl->pos, readLength) < rightMatch.netPosSrcPg(srcRl->pos, readLength)) {
            uint_pg_len_max overlapLength = leftMatch.netEndPosSrcPg(srcRl->pos) - rightMatch.netPosSrcPg(srcRl->pos, readLength);
            rightMatch.trimLeft(overlapLength, srcRl->pos, readLength);
            collidedCharsCount += overlapLength;
            return !rightMatch.inactive(srcRl->pos, readLength);
        }
        if (leftMatch.netSrcLength(srcRl->pos, readLength) > rightMatch.netSrcLength(srcRl->pos, readLength)) {
            collidedCharsCount += rightMatch.netSrcLength(srcRl->pos, readLength);
            rightMatch.endRlIdx = rightMatch.startRlIdx - 1;
        } else {
            collidedCharsCount += leftMatch.netSrcLength(srcRl->pos, readLength);
            leftMatch.endRlIdx = leftMatch.startRlIdx - 1;
        }
        return false;
    }

    void DefaultPgMatcher::fillPgMatches() {
        pgMatches.clear();
        pgMatches.reserve(textMatches.size());
        for (TextMatch& textMatch: textMatches) {
            pgMatches.push_back(textMatch);
        }
    }

    void DefaultPgMatcher::correctDestPositionDueToRevComplMatching() {
        for (TextMatch& match: textMatches)
            match.posDestText = srcPgh->getPseudoGenomeLength() - (match.posDestText + match.length);
    }

    string DefaultPgMatcher::getTotalMatchStat(uint_pg_len_max totalMatchLength) {
        return toString(totalMatchLength) + " (" + toString((totalMatchLength * 100.0) / srcPgh->getPseudoGenomeLength(), 1)+ "%)";
    }

    void DefaultPgMatcher::transferMatchesFromSrcToDest(const string &destPgPrefix) {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.netPosSrcPg(srcRl->pos, readLength) < pgMatch2.netPosSrcPg(srcRl->pos, readLength); });

        vector<uint32_t> matchDestOrderedIdx;
        if (destPgIsSrcPg) {
            matchDestOrderedIdx.resize(pgMatches.size());
            for (uint32_t i = 0; i < pgMatches.size(); i++)
                matchDestOrderedIdx[i] = i;
            sort(matchDestOrderedIdx.begin(), matchDestOrderedIdx.end(),
                 [this](const uint32_t &pgMatchIdx1, const uint32_t &pgMatchIdx2) -> bool {
                     return pgMatches[pgMatchIdx1].netPosDestPg(srcRl->pos, readLength, revComplMatching)
                     < pgMatches[pgMatchIdx2].netPosDestPg(srcRl->pos, readLength, revComplMatching);
                 });
        }
        string newSrcPg;
        vector<bool> isReadRemapped(srcPgh->getReadsCount() + 1, false);
        newRlPos.resize(srcPgh->getReadsCount() + 1);

        uint64_t dOrdIdx = 0;
        uint_pg_len_max pos = 0;
        uint_pg_len_max lastReadPos = 0;
        uint_reads_cnt_max i = 0;
        uint_pg_len_max removedCount = 0;
        for (PgMatch& pgMatch: pgMatches) {
            if (pgMatch.inactive(srcRl->pos, readLength))
                continue;
            if (i >= pgMatch.startRlIdx) {
                removedCount += pgMatch.netEndPosSrcPg(srcRl->pos) - pos;
            } else {
                uint_read_len_max offset = srcRl->pos[i] - removedCount - lastReadPos;
                newSrcPg.append(srcPg, pos, pgMatch.netPosSrcPg(srcRl->pos, readLength) - pos);
                while (i < pgMatch.startRlIdx) {
                    newRlPos[i] = srcRl->pos[i] - removedCount;
                    lastReadPos = newRlPos[i++];
                }
                if (destPgIsSrcPg) {
                    while (dOrdIdx < pgMatches.size() &&
                           pgMatches[matchDestOrderedIdx[dOrdIdx]].netPosDestPg(srcRl->pos, readLength, revComplMatching) <
                           pgMatch.netPosSrcPg(srcRl->pos, readLength)) {
                        PgMatch &destMatch = pgMatches[matchDestOrderedIdx[dOrdIdx++]];
                        if (!destMatch.inactive(srcRl->pos, readLength))
                            destMatch.mapping.posDestText -= removedCount;
                    }
                }
                removedCount += pgMatch.netSrcLength(srcRl->pos, readLength);
            }
            pos = pgMatch.netEndPosSrcPg(srcRl->pos);
            i = pgMatch.endRlIdx + 1;
        }
        newSrcPg.append(srcPg, pos, srcPg.length() - pos);
        cout << "Source Pg reduced to " << newSrcPg.length() << " symbols (removed: " <<
            getTotalMatchStat(srcPgh->getPseudoGenomeLength() - newSrcPg.length()) << ")." << endl;
        if (destPgIsSrcPg)
            while (dOrdIdx < pgMatches.size()) {
                PgMatch& destMatch = pgMatches[matchDestOrderedIdx[dOrdIdx++]];
                if (!destMatch.inactive(srcRl->pos, readLength))
                    destMatch.mapping.posDestText -= removedCount;
            }

        for (PgMatch& pgMatch: pgMatches) {
            if (pgMatch.inactive(srcRl->pos, readLength))
                continue;
            for(uint_reads_cnt_max i = pgMatch.startRlIdx; i <= pgMatch.endRlIdx; i++) {
                isReadRemapped[i] = true;
                newRlPos[i] = pgMatch.mapSrcReadToDest(srcRl->pos[i], readLength, revComplMatching);
            }
        }
        srcRl->pos.clear();

        vector<uint_reads_cnt_max> rlPosOrd(srcPgh->getReadsCount());
        for (uint_reads_cnt_max i = 0; i < srcPgh->getReadsCount(); i++)
            rlPosOrd[i] = i;
        sort(rlPosOrd.begin(), rlPosOrd.end(),
             [this](const uint_reads_cnt_max &rlPosIdx1, const uint_reads_cnt_max &rlPosIdx2) -> bool {
                 return newRlPos[rlPosIdx1] < newRlPos[rlPosIdx2];
             });

        SeparatedPseudoGenomeOutputBuilder builder(destPgIsSrcPg?destPgPrefix:srcPgPrefix);
        builder.copyPseudoGenomeProperties(srcPgPrefix);

        pos = 0;
        uint_pg_len_max totalOffsetOverflow = 0;
        DefaultReadsListEntry entry;
        for(uint_reads_cnt_max iOrd = 0; iOrd < srcPgh->getReadsCount(); iOrd++) {
            uint_reads_cnt_max i = rlPosOrd[iOrd];
            if(destPgIsSrcPg || !isReadRemapped[i]) {
                entry.advanceEntryByPosition(newRlPos[i] - totalOffsetOverflow, srcRl->orgIdx[i],
                        srcRl->revComp[i] != (destPgIsSrcPg?(revComplMatching && isReadRemapped[i]):false));
                srcRl->copyMismatchesToEntry(i, entry);
                if (entry.offset > readLength) {
                    uint_pg_len_max overflow = entry.offset - readLength;
                    totalOffsetOverflow += overflow;
                    entry.offset = readLength;
                    entry.pos -= overflow;
                    builder.appendPseudoGenome(newSrcPg.substr(pos, newRlPos[i] - overflow - pos));
                    pos = newRlPos[i];
                }
                builder.writeExtraReadEntry(entry);
            }
        }
        builder.appendPseudoGenome(newSrcPg.substr(pos, newSrcPg.length() - pos));
        cout << "Final size of Pg: " << (newSrcPg.length() - totalOffsetOverflow) << " (removed "
            << totalOffsetOverflow << " chars in overflowed offsets)" << endl;
        builder.build();

        if(!destPgIsSrcPg) {
            cout << "Error: Unimplemented adding reads removed from source Pg to destination Pg.";
            exit(EXIT_FAILURE);
        }
        newRlPos.clear();
    }

    void DefaultPgMatcher::matchPgInPgFile(const string &srcPgPrefix, const string &targetPgPrefix, uint_pg_len_max targetMatchLength,
                         const string &destPgPrefix, bool revComplPg, bool dumpInfo) {

        PgTools::DefaultPgMatcher matcher(srcPgPrefix, targetPgPrefix, revComplPg);

        matcher.exactMatchPg(targetMatchLength);

        if (dumpInfo)
            matcher.writeMatchesInfo(destPgPrefix);

        matcher.transferMatchedReads(destPgPrefix);
    }
}

