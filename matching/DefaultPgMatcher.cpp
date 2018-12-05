#include "DefaultPgMatcher.h"

#include "../pseudogenome/DefaultPseudoGenome.h"
#include "../pseudogenome/PackedPseudoGenome.h"
#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../pseudogenome/persistence/SeparatedExtendedReadsListIterator.h"
#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include <list>

namespace PgTools {

    DefaultPgMatcher::DefaultPgMatcher(const string& srcPgPrefix, const string& targetPgPrefix, bool revComplMatching)
            :srcPgPrefix(srcPgPrefix), targetPgPrefix(targetPgPrefix), revComplMatching(revComplMatching) {
        textFromSamePg = srcPgPrefix == targetPgPrefix;
        if (textFromSamePg)
            cout << "Reading pseudogenome..." << endl;
        else
            cout << "Reading source pseudogenome..." << endl;
        PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(srcPgPrefix, srcPgh, plainTextReadMode);
        srcPg = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(srcPgPrefix);
        cout << "Pseudogenome length: " << srcPgh->getPseudoGenomeLength() << endl;
        readLength = srcPgh->getMaxReadLength();

        if (targetPgPrefix != srcPgPrefix) {
            cout << "Reading target pseudogenome..." << endl;
            PseudoGenomeHeader *pgh = 0;
            bool plainTextReadMode = false;
            PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(targetPgPrefix, pgh,
                                                                                 plainTextReadMode);
            destPg = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(targetPgPrefix);
            cout << "Pseudogenome length: " << pgh->getPseudoGenomeLength() << endl;
        } else
            destPg = srcPg;

        if (revComplMatching)
            destPg = PgSAHelpers::reverseComplement(destPg);
    }

    DefaultPgMatcher::~DefaultPgMatcher() {
    }

    void backMatchExpand(const string& text, uint64_t& textPos, const string& pattern, uint64_t& patternPos, uint64_t& length) {
        const char* textGuardPtr = text.data();
        const char* textPtr = text.data() + textPos;
        const char* patternGuardPtr = pattern.data();
        const char* patternPtr = pattern.data() + patternPos;
        uint64_t i = 0;
        while (textPtr-- != textGuardPtr && patternPtr-- != patternGuardPtr && *textPtr == *patternPtr)
            i++;
        textPos -= i;
        patternPos -= i;
        length += i;
    }

    void forwardMatchExpand(const string& text, uint64_t textPos, const string& pattern, uint64_t patternPos, uint64_t& length) {
        const char* textGuardPtr = text.data() + text.length();
        const char* textPtr = text.data() + textPos + length;
        const char* patternGuardPtr = pattern.data() + pattern.length();
        const char* patternPtr = pattern.data() + patternPos + length;
        while (*textPtr == *patternPtr && textPtr++ != textGuardPtr && patternPtr++ != patternGuardPtr)
            length++;
    }

    void DefaultPgMatcher::exactMatchPg(uint32_t minMatchLength) {
        clock_checkpoint();

        cout << "Feeding pattern pseudogenome parts...\n" << endl;
        const uint_pg_len_max matchingLength = minMatchLength / 2;
        DefaultConstantLengthPatternsOnTextHashMatcher hashMatcher(matchingLength);

        const char* pgPtr = srcPg.data();
        uint32_t partsCount = srcPg.length() / matchingLength;
        for (uint32_t i = 0; i < partsCount; i++)
            hashMatcher.addPattern(pgPtr + i * matchingLength, i);

        cout << "... finished in " << clock_millis() << " msec. " << endl;
        clock_checkpoint();
        cout << "Matching...\n" << endl;

        const char* textPtr = destPg.data();
        hashMatcher.iterateOver(textPtr, destPg.length());

        pgMatches.clear();
        uint_pg_len_max furthestMatchEndPos = 0;
        list<PgMatch> currentMatches;

        int i = 0;
        uint64_t matchCount = 0;
        uint64_t matchCharsCount = 0;
        uint64_t matchCharsWithOverlapCount = 0;
        uint64_t shorterMatchCount = 0;
        uint64_t falseMatchCount = 0;
        while (hashMatcher.moveNext()) {
            uint64_t matchDestPos = hashMatcher.getHashMatchTextPosition();
            const uint32_t matchPatternIndex = hashMatcher.getHashMatchPatternIndex();
            uint64_t matchSrcPos = matchPatternIndex * matchingLength;
            if(textFromSamePg && revComplMatching?destPg.length() - matchSrcPos < matchDestPos:matchDestPos >= matchSrcPos)
                continue;
            auto cmIt = currentMatches.begin();
            bool continueMatch = false;
            while (cmIt != currentMatches.end()) {
                if (matchDestPos + matchingLength > (*cmIt).endGrossPosDestPg())  {
                    currentMatches.erase(cmIt++);
                } else {
                    if ((*cmIt).posGrossSrcPg - (*cmIt).posGrossDestPg == matchSrcPos - matchDestPos) {
                        continueMatch = true;
                        break;
                    }
                    cmIt++;
                }
            }
            if (continueMatch)
                continue;

            bool confirmPatternMatch = strncmp(textPtr + matchDestPos, pgPtr + matchSrcPos, matchingLength) == 0;
            uint64_t matchLength = matchingLength;
            if (confirmPatternMatch) {
                backMatchExpand(destPg, matchDestPos, srcPg, matchSrcPos, matchLength);
                forwardMatchExpand(destPg, matchDestPos, srcPg, matchSrcPos, matchLength);
                if (matchLength >= minMatchLength) {
                    matchCount++;
                    matchCharsWithOverlapCount += matchLength;
                    matchCharsCount += matchLength - (matchDestPos < furthestMatchEndPos?
                                                      furthestMatchEndPos - matchDestPos:0);
                    const PgMatch &matchInfo = PgMatch(matchSrcPos, matchLength, matchDestPos);
                    if (furthestMatchEndPos < matchDestPos + matchLength)
                        furthestMatchEndPos = matchDestPos + matchLength;
                    pgMatches.push_back(matchInfo);
                    currentMatches.push_back(matchInfo);
                } else
                    shorterMatchCount++ ;
            } else
                falseMatchCount++;

            if ((i++ % 1000000) == 0 || matchLength > minMatchLength * 15) {
                 cout << (confirmPatternMatch?"Matched ":"False-matched ") << matchLength << " chars: <" << matchSrcPos << "; "
                     << (matchSrcPos + matchLength) << ") " << " in " << matchDestPos << " (" << matchCharsCount << " chars matched total)" << endl;
                 cout << "Elapsed time: " << ((double) destPg.length() / matchDestPos) * clock_millis() / 1000 << "[s]" << endl;
                if (matchLength < matchingLength * 10) {
                    uint8_t beforeChars = matchingLength / 4;
                    const uint64_t &minPos = min<uint64_t>(matchPatternIndex, matchDestPos);
                    beforeChars = minPos < beforeChars ? minPos : beforeChars;
                    uint8_t afterChars = matchingLength / 4;
                    const uint64_t &maxLength =
                            min<uint64_t>(srcPg.length() - matchSrcPos, destPg.length() - matchDestPos) -
                            matchLength;
                    afterChars = maxLength < afterChars ? maxLength : afterChars;
                    cout << srcPg.substr(matchSrcPos, matchLength) << endl;
                    cout << destPg.substr(matchDestPos, matchLength) << endl;
                    cout << srcPg.substr(matchSrcPos - beforeChars, beforeChars) << " ... " << srcPg.substr(matchSrcPos + matchLength, afterChars) << endl;
                    cout << destPg.substr(matchDestPos - beforeChars, beforeChars) << " ... " << destPg.substr(matchDestPos + matchLength, afterChars) << endl;
                }
            }
        }
        cout << "Finished matching in  " << clock_millis() << " msec. " << endl;
        cout << "Exact matched " << matchCount << " parts (too short matches: " << shorterMatchCount << "). False matches reported: " << falseMatchCount << "." << endl;
        cout << "Stage 1: Matched " << getTotalMatchStat(matchCharsCount) << " Pg characters. Sum length of all matches is " << matchCharsWithOverlapCount << "." << endl;

        std::sort(pgMatches.begin(), pgMatches.end(), [](const PgMatch& match1, const PgMatch& match2) -> bool
        { return match1.grossLength > match2.grossLength; });

        cout << "Largest matches:" << endl;
        for(uint32_t i = 0; i < pgMatches.size() && i < 10; i++)
            pgMatches[i].report(cout);
    }

    void DefaultPgMatcher::writeMatchesInfo(const string &dumpFilePrefix) {
        cout << "Sorry, unimplemented matches info dump feature." << endl;
    }

    using namespace PgTools;

    void DefaultPgMatcher::writeIntoPseudoGenome(const string &destPgFilePrefix) {
        fillSrcReadsList();

        if (revComplMatching)
            correctDestPositionDueToRevComplMatching();

        mapPgMatches2SrcReadsList();

        if (textFromSamePg) {
            reverseDestWithSrcForBetterMatchesMappingInTheSamePg();
            resolveDestOverlapSrcConflictsInTheSamePg();
        }
/*
        resolveMatchesOverlapInSrc()
        removeMatchesFromSrc();
        if (textFromSamePg)
            addMatchesToDestWithRemovalCorrection(destPgFilePrefix);
        else
            addMatchesToDest(destPgFilePrefix);
*/
        rlPos.clear();
        rlIdx.clear();
    }

    void DefaultPgMatcher::mapPgMatches2SrcReadsList() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.posGrossSrcPg < pgMatch2.posGrossSrcPg; });

        uint32_t cancelledMatchesCount = 0;
        uint_pg_len_max totalNetMatchCharsWithOverlapCount = 0;
        int64_t i = 0;
        for(PgMatch& pgMatch: pgMatches) {
            pgMatch.alignToReads(rlPos, readLength, i);
            if (pgMatch.inactive(rlPos, readLength))
                cancelledMatchesCount++;
            totalNetMatchCharsWithOverlapCount += pgMatch.netSrcLength(rlPos, readLength);
        }
        cout << "Cancelled " << cancelledMatchesCount << " short matches due to external reads overlap." << endl;
        cout << "Stage 2: Net-matched " << getTotalMatchStat(totalNetMatchCharsWithOverlapCount) << " Pg src->dest sum length of all matches." << endl;
    }

    void DefaultPgMatcher::reverseDestWithSrcForBetterMatchesMappingInTheSamePg() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.posGrossDestPg < pgMatch2.posGrossDestPg; });

        uint_read_len_max readLength = srcPgh->getMaxReadLength();
        uint_pg_len_max totalNetReverseMatchCharsWithOverlapCount = 0;
        uint_pg_len_max totalNetMatchCharsWithOverlapCount = 0;
        uint32_t restoredMatchesCount = 0;
        int64_t i = 0;
        for (PgMatch& pgMatch: pgMatches) {
            PgMatch revMatch(pgMatch.posGrossDestPg, pgMatch.grossLength, pgMatch.posGrossSrcPg);
            revMatch.alignToReads(rlPos, readLength, i);
            uint_pg_len_max revNetSrcLength = revMatch.netSrcLength(rlPos, readLength);
            totalNetReverseMatchCharsWithOverlapCount += revNetSrcLength;
            if (revNetSrcLength > pgMatch.netSrcLength(rlPos, readLength)) {
                if (pgMatch.netSrcLength(rlPos, readLength) == 0)
                    restoredMatchesCount++;
                pgMatch = revMatch;
            }
            totalNetMatchCharsWithOverlapCount += pgMatch.netSrcLength(rlPos, readLength);
        }

        cout << "Restored " << restoredMatchesCount << " short matches (cancelled earlier due to external reads overlap)." << endl;
        cout << "Stage 2a: Net-matched " << getTotalMatchStat(totalNetReverseMatchCharsWithOverlapCount) << " Pg dest->src sum length of all matches." << endl;
        cout << "Stage 2b: Optimal net-matched " << getTotalMatchStat(totalNetMatchCharsWithOverlapCount) << " Pg sum length of all matches." << endl;
    }

    void DefaultPgMatcher::resolveDestOverlapSrcConflictsInTheSamePg() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.netEndPosSrcPg(rlPos) < pgMatch2.netEndPosSrcPg(rlPos); });

        vector<uint32_t> matchDestOrderedIdx(pgMatches.size());
        for(uint32_t i = 0; i < pgMatches.size(); i++)
            matchDestOrderedIdx[i] = i;
        sort(matchDestOrderedIdx.begin(), matchDestOrderedIdx.end(), [this](const uint32_t& pgMatchIdx1, const uint32_t& pgMatchIdx2) -> bool
        { return pgMatches[pgMatchIdx1].netPosDestPg(rlPos) < pgMatches[pgMatchIdx2].netPosDestPg(rlPos); });

        uint_pg_len_max totalNetMatchCharsCount = 0;
        uint_pg_len_max collidedCharsCount = 0;
        uint32_t conflictsCount = 0;
        uint32_t resolvedConflictsCount = 0;
        vector<uint_pg_len_max> maxNetEndPosDestPgUpTo(pgMatches.size());
        maxNetEndPosDestPgUpTo[0] = pgMatches[matchDestOrderedIdx[0]].netEndPosDestPg(rlPos, readLength);
        for(uint32_t i = 1; i < pgMatches.size(); i++) {
            maxNetEndPosDestPgUpTo[i] = pgMatches[matchDestOrderedIdx[i]].netEndPosDestPg(rlPos, readLength);
            if (maxNetEndPosDestPgUpTo[i] < maxNetEndPosDestPgUpTo[i - 1])
                maxNetEndPosDestPgUpTo[i] = maxNetEndPosDestPgUpTo[i - 1];
        }
        int64_t dOrdIdx = 0;
        for (PgMatch& sPgMatch: pgMatches) {
            totalNetMatchCharsCount += sPgMatch.netSrcLength(rlPos, readLength);
            while(sPgMatch.netPosSrcPg(rlPos, readLength) < maxNetEndPosDestPgUpTo[dOrdIdx] && --dOrdIdx > 0);
            while (++dOrdIdx < pgMatches.size()
                && pgMatches[matchDestOrderedIdx[dOrdIdx]].netPosDestPg(rlPos) < sPgMatch.netEndPosSrcPg(rlPos)) {
                if (sPgMatch.inactive(rlPos, readLength))
                    break;
                int64_t dIdx = matchDestOrderedIdx[dOrdIdx];
                if (!pgMatches[dIdx].inactive(rlPos, readLength)
                    && pgMatches[dIdx].netEndPosDestPg(rlPos, readLength) > sPgMatch.netPosSrcPg(rlPos, readLength)) {
//                    cout << pgMatches[dIdx].netPosDestPg(rlPos) << "\t" << pgMatches[dIdx].netEndPosDestPg(rlPos, readLength) <<
//                    "\t" << sPgMatch.netPosSrcPg(rlPos, readLength) << "\t" << sPgMatch.netEndPosSrcPg(rlPos) << endl;
                    if (resolveCollision(pgMatches[dIdx], sPgMatch, collidedCharsCount))
                        resolvedConflictsCount++;
                    else
                        conflictsCount++;
                }
            }
        }
        totalNetMatchCharsCount -= collidedCharsCount;
        cout << conflictsCount << " conflicts remain in matches (resolved " << resolvedConflictsCount << ")." << endl;
        cout << "Stage 2c: Conflicts filtered optimal net-matched " << getTotalMatchStat(totalNetMatchCharsCount) << " Pg sum length of all matches." << endl;
    }

    bool DefaultPgMatcher::resolveCollision(PgMatch &destMatch, PgMatch &srcMatch, uint_pg_len_max& collidedCharsCount) {
        if (destMatch.netEndPosDestPg(rlPos, readLength) < srcMatch.netEndPosSrcPg(rlPos) &&
            destMatch.netPosDestPg(rlPos) < srcMatch.netPosSrcPg(rlPos, readLength)) {
            uint64_t overlapLength = destMatch.netEndPosDestPg(rlPos, readLength) - srcMatch.netPosSrcPg(rlPos, readLength);
            int64_t charsReduced = trimOverlap(destMatch, srcMatch, overlapLength);
            if (charsReduced != 0) {
                if (charsReduced < 0) {
                    charsReduced = -charsReduced;
                } else {

                }
                collidedCharsCount += charsReduced;
                return true;
            }
        }
        if (destMatch.netEndPosDestPg(rlPos, readLength) > srcMatch.netEndPosSrcPg(rlPos) &&
            destMatch.netPosDestPg(rlPos) > srcMatch.netPosSrcPg(rlPos, readLength)) {
            uint64_t overlapLength = srcMatch.netEndPosSrcPg(rlPos) - destMatch.netPosDestPg(rlPos);
            int64_t charsReduced = trimOverlap(srcMatch, destMatch, overlapLength);
            if (charsReduced != 0) {
                if (charsReduced < 0) {
                    charsReduced = -charsReduced;
                } else {

                }
                collidedCharsCount += charsReduced;
                return true;
            }
        }
        if (destMatch.netSrcLength(rlPos, readLength) > srcMatch.netSrcLength(rlPos, readLength)) {
            collidedCharsCount += srcMatch.netSrcLength(rlPos, readLength);
            srcMatch.endRlIdx = srcMatch.startRlIdx - 1;
        } else {
            collidedCharsCount += destMatch.netSrcLength(rlPos, readLength);
            destMatch.endRlIdx = destMatch.startRlIdx - 1;
        }
        return false;
    }

    int64_t DefaultPgMatcher::trimOverlap(PgMatch &leftMatch, PgMatch &rightMatch, uint_pg_len_max overlapLength) {
        uint64_t leftCharsReduced = 0;
        bool leftTrim = false;
        uint_reads_cnt_max leIdx = leftMatch.endRlIdx + 1;
        while (--leIdx > leftMatch.startRlIdx && rlPos[leftMatch.endRlIdx + 1] - rlPos[leIdx] < overlapLength);
        if (leIdx > leftMatch.startRlIdx) {
            leftCharsReduced = rlPos[leftMatch.endRlIdx + 1] - rlPos[leIdx];
            leftTrim = true;
        }

        uint_reads_cnt_max rsIdx = rightMatch.startRlIdx - 1;
        while (++rsIdx < rightMatch.endRlIdx && rlPos[rsIdx] - rlPos[rightMatch.startRlIdx - 1] < overlapLength);
        if (rsIdx < rightMatch.endRlIdx) {
            int64_t rightCharsReduced = rlPos[rsIdx] - rlPos[rightMatch.startRlIdx - 1];
            if ((!leftTrim && rightCharsReduced < leftMatch.netSrcLength(rlPos, readLength))
                ||(leftTrim && rightCharsReduced < leftCharsReduced)) {
                rightMatch.startRlIdx = rsIdx + 1;
                return rightCharsReduced;
            }
        }
        if (leftTrim && leftCharsReduced < rightMatch.netSrcLength(rlPos, readLength)) {
            leftMatch.endRlIdx = leIdx - 1;
        }
        return -((int64_t) leftCharsReduced);
    }

    void DefaultPgMatcher::fillSrcReadsList() {
        rlPos.clear();
        rlIdx.clear();
        rlPos.reserve(srcPgh->getReadsCount());
        rlIdx.reserve(srcPgh->getReadsCount());
        SimpleSeparatedReadsListIterator* rlIt = new SimpleSeparatedReadsListIterator(srcPgPrefix);
        while (rlIt->moveNext()) {
            rlPos.push_back(rlIt->peekReadEntry().pos);
            rlIdx.push_back(rlIt->peekReadEntry().idx);
        }
        rlPos.push_back(srcPgh->getPseudoGenomeLength());
        delete(rlIt);
    }


    void DefaultPgMatcher::correctDestPositionDueToRevComplMatching() {
        for (PgMatch& pgMatch: pgMatches)
            pgMatch.posGrossDestPg = srcPgh->getPseudoGenomeLength() - (pgMatch.posGrossDestPg + pgMatch.grossLength);
    }

    string DefaultPgMatcher::getTotalMatchStat(uint_pg_len_max totalMatchLength) {
        return toString(totalMatchLength) + " (" + toString((totalMatchLength * 100.0) / srcPgh->getPseudoGenomeLength(), 1)+ "%)";
    }
}

