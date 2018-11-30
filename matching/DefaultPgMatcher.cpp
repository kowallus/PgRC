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
                if (matchDestPos + matchingLength > (*cmIt).endPosDestPg())  {
                    currentMatches.erase(cmIt++);
                } else {
                    if ((*cmIt).posSrcPg - (*cmIt).posDestPg == matchSrcPos - matchDestPos) {
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

            if ((i++ % 100000) == 0 || matchLength > minMatchLength * 15) {
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
        cout << "Stage 1: Matched " << matchCharsCount << " Pg characters. Sum length of all matches is " << matchCharsWithOverlapCount << "." << endl;

        std::sort(pgMatches.begin(), pgMatches.end(), [](const PgMatch& match1, const PgMatch& match2) -> bool
        { return match1.length > match2.length; });

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
        { return pgMatch1.posSrcPg < pgMatch2.posSrcPg; });

        uint_read_len_max readLength = srcPgh->getMaxReadLength();
        uint_pg_len_max totalNetMatchCharsWithOverlapCount = 0;
        int64_t i = -1;
        for(PgMatch& pgMatch: pgMatches) {
            while(i > 0 && rlPos[--i] >= pgMatch.posSrcPg);
            while(rlPos[++i] < pgMatch.posSrcPg);
            pgMatch.startRlIdx = i--;
            while(rlPos[++i] + readLength <= pgMatch.posSrcPg + pgMatch.length);
            pgMatch.endRlIdx = i - 1;
            uint_pg_len_max nettoPosSrcPg = pgMatch.startRlIdx > 0?rlPos[pgMatch.startRlIdx - 1] + readLength:0;
            pgMatch.netPosAlignment = nettoPosSrcPg - pgMatch.posSrcPg;
            pgMatch.netLength = rlPos[i] - nettoPosSrcPg;

            totalNetMatchCharsWithOverlapCount += pgMatch.netLength;
        }
        cout << "Stage 2: Net-matched " << totalNetMatchCharsWithOverlapCount << " Pg src->dest sum length of all matches." << endl;
    }

    void DefaultPgMatcher::reverseDestWithSrcForBetterMatchesMappingInTheSamePg() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.posDestPg < pgMatch2.posDestPg; });

        uint_read_len_max readLength = srcPgh->getMaxReadLength();
        uint_pg_len_max totalNetReverseMatchCharsWithOverlapCount = 0;
        uint_pg_len_max totalNetMatchCharsWithOverlapCount = 0;
        int64_t i = -1;
        for (PgMatch& pgMatch: pgMatches) {
            while (i > 0 && rlPos[--i] >= pgMatch.posDestPg);
            while (rlPos[++i] < pgMatch.posDestPg);
            uint_reads_cnt_max startRlIdx = i--;
            while (rlPos[++i] + readLength <= pgMatch.posDestPg + pgMatch.length);
            uint_reads_cnt_max endRlIdx = i - 1;
            uint_pg_len_max nettoPosDestPg = startRlIdx > 0?rlPos[startRlIdx - 1] + readLength:0;
            uint_pg_len_max netLength = rlPos[i] - nettoPosDestPg;
            totalNetReverseMatchCharsWithOverlapCount += netLength;
            if (netLength > pgMatch.netLength) {
                pgMatch.netPosAlignment = nettoPosDestPg - pgMatch.posDestPg;
                pgMatch.reverseMatch();
                pgMatch.startRlIdx = startRlIdx;
                pgMatch.endRlIdx = endRlIdx;
                pgMatch.netLength = netLength;
            }
            totalNetMatchCharsWithOverlapCount += pgMatch.netLength;
        }

        cout << "Stage 2a: Net-matched " << totalNetReverseMatchCharsWithOverlapCount << " Pg dest->src sum length of all matches." << endl;
        cout << "Stage 2b: Optimal net-matched " << totalNetMatchCharsWithOverlapCount << " Pg sum length of all matches." << endl;
    }

    void DefaultPgMatcher::resolveDestOverlapSrcConflictsInTheSamePg() {
        sort(pgMatches.begin(), pgMatches.end(), [this](const PgMatch& pgMatch1, const PgMatch& pgMatch2) -> bool
        { return pgMatch1.netEndPosSrcPg() < pgMatch2.netEndPosSrcPg(); });

        vector<uint32_t> matchDestOrderedIdx(pgMatches.size());
        for(uint32_t i = 0; i < pgMatches.size(); i++)
            matchDestOrderedIdx[i] = i;
        sort(matchDestOrderedIdx.begin(), matchDestOrderedIdx.end(), [this](const uint32_t& pgMatchIdx1, const uint32_t& pgMatchIdx2) -> bool
        { return pgMatches[pgMatchIdx1].netPosDestPg() < pgMatches[pgMatchIdx2].netPosDestPg(); });

        uint_pg_len_max totalNetMatchCharsCount = 0;
        uint_pg_len_max collidedCharsCount = 0;
        uint32_t collisionsCount = 0;
        int64_t destFirstGuardIdx = 0;
        for (PgMatch& pgMatch: pgMatches) {
            totalNetMatchCharsCount += pgMatch.netLength;
            if (pgMatch.inactiveDueToCollision)
                continue;
            bool foundFirst = false;
            int64_t dOrdIdx = destFirstGuardIdx - 1;
            while (++dOrdIdx < pgMatches.size() && pgMatches[matchDestOrderedIdx[dOrdIdx]].netPosDestPg() < pgMatch.netEndPosSrcPg()) {
                int64_t dIdx = matchDestOrderedIdx[dOrdIdx];
                if (!pgMatches[dIdx].inactiveDueToCollision && pgMatches[dIdx].netEndPosDestPg() > pgMatch.netEndPosSrcPg()) {
                    if (!foundFirst) destFirstGuardIdx = dIdx;
                    collisionsCount++;
                    if (pgMatches[dIdx].netLength > pgMatch.netLength) {
                        pgMatch.inactiveDueToCollision = true;
                        collidedCharsCount += pgMatch.netLength;
                    } else {
                        pgMatches[dIdx].inactiveDueToCollision = true;
                        collidedCharsCount += pgMatches[dIdx].netLength;
                    }
                }
            }
        }
        totalNetMatchCharsCount -= collidedCharsCount;
        cout << "Found " << collisionsCount << " collisions in matches." << endl;
        cout << "Stage 2c: Optimal collision filtered net-matched " << totalNetMatchCharsCount << " Pg sum length of all matches." << endl;

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
            pgMatch.posDestPg = srcPgh->getPseudoGenomeLength() - (pgMatch.posDestPg + pgMatch.length);
    }
}

