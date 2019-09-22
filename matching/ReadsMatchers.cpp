#include "ReadsMatchers.h"

#include "../readsset/persistance/ReadsSetPersistence.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../pseudogenome/readslist/SeparatedExtendedReadsList.h"
#include <omp.h>

namespace PgTools {

    uint8_t countMismatches(const char *pattern, const char *text, uint64_t length, uint8_t maxMismatches) {
        uint8_t res = 0;
        const char *patEnd = pattern + length;
        while (pattern != patEnd) {
            if (*pattern++ != *text++) {
                if (res++ >= maxMismatches)
                    return NOT_MATCHED_COUNT;
            }
        }
        return res;
    }

    void reportMismatches(const char *read, const char *pgPart, const uint_read_len_max length, ofstream &mismatchesDest) {
        uint64_t pos = 0;
        do {
            if (read[pos] != pgPart[pos])
                mismatchesDest << "\t" << pos << "\t" << read[pos];
        } while (++pos < length);
    }

    void fillEntryWithMismatches(const char *read, const char *pgPart,
            const uint8_t mismatchesCount, DefaultReadsListEntry &entry) {
        uint64_t pos = 0;
        uint8_t count = 0;
        while (count < mismatchesCount) {
            if (read[pos] != pgPart[pos]) {
                entry.addMismatch(mismatch2code(pgPart[pos], read[pos]), pos);
                count++;
            }
            pos++;
        };
    }

    void fillEntryWithReversedMismatches(const char *read, const char *pgPart,
                                 const uint8_t mismatchesCount, DefaultReadsListEntry &entry,
                                 uint_read_len_max readLength) {
        uint64_t pos = readLength;
        uint8_t count = 0;
        while (count < mismatchesCount) {
            if (read[--pos] != pgPart[pos]) {
                entry.addMismatch(mismatch2code(reverseComplement(pgPart[pos]), reverseComplement(read[pos])),
                        readLength - pos - 1);
                count++;
            }
        };
    }

    const uint_read_len_max DefaultReadsMatcher::DISABLED_PREFIX_MODE = (uint_read_len_max) -1;
    const uint64_t DefaultReadsMatcher::NOT_MATCHED_POSITION = UINT64_MAX;

    DefaultReadsMatcher::DefaultReadsMatcher(char* pgPtr, const uint_pg_len_max pgLength,
                                             bool revComplPg, ConstantLengthReadsSetInterface *readsSet,
                                             uint32_t matchPrefixLength) :
                                             pgPtr(pgPtr), pgLength(pgLength), revComplPg(revComplPg),
                                             readsSet(readsSet), matchPrefixLength(matchPrefixLength),
                                             readLength(readsSet->maxReadLength()),
                                             matchingLength(readLength > matchPrefixLength ? matchPrefixLength : readLength),
                                             readsCount(readsSet->readsCount()) { }

    DefaultReadsMatcher::~DefaultReadsMatcher() {
    }

    DefaultReadsExactMatcher::~DefaultReadsExactMatcher() {
        delete(hashMatcher);
    }

    AbstractReadsApproxMatcher::~AbstractReadsApproxMatcher() {}

    DefaultReadsApproxMatcher::~DefaultReadsApproxMatcher() {
        delete(hashMatcher);
    }

    InterleavedReadsApproxMatcher::~InterleavedReadsApproxMatcher() {
        delete(hashMatcher);
    }

    void DefaultReadsMatcher::initMatching() {
        readMatchPos.clear();
        readMatchPos.insert(readMatchPos.begin(), readsCount, NOT_MATCHED_POSITION);
        readMatchRC.clear();
        readMatchRC.insert(readMatchRC.begin(), readsCount, false);
        matchedReadsCount = 0;
        betterMatchCount = 0;
        falseMatchCount = 0;
    }

    void AbstractReadsApproxMatcher::initMatchingContinuation(DefaultReadsMatcher *pMatcher) {
        pMatcher->transferMatchingResults(this);
    }

    void DefaultReadsMatcher::transferMatchingResults(AbstractReadsApproxMatcher *approxMatcher) {
        approxMatcher->readMatchPos = std::move(readMatchPos);
        approxMatcher->readMatchRC = std::move(readMatchRC);
        approxMatcher->matchedReadsCount = matchedReadsCount;
        approxMatcher->betterMatchCount = betterMatchCount;
        approxMatcher->falseMatchCount = falseMatchCount;
    }

    void AbstractReadsApproxMatcher::transferMatchingResults(AbstractReadsApproxMatcher *approxMatcher) {
        DefaultReadsMatcher::transferMatchingResults(approxMatcher);
        approxMatcher->readMismatchesCount = std::move(readMismatchesCount);
        memcpy(approxMatcher->matchedCountPerMismatches, matchedCountPerMismatches, sizeof(uint_reads_cnt_max) * (NOT_MATCHED_COUNT + 1));

    }

    void DefaultReadsExactMatcher::transferMatchingResults(AbstractReadsApproxMatcher *approxMatcher) {
        approxMatcher->readMismatchesCount.resize(readsCount);
        for(uint_reads_cnt_max i = 0; i < readsCount; i++)
            approxMatcher->readMismatchesCount[i] = readMatchPos[i] == NOT_MATCHED_POSITION?NOT_MATCHED_COUNT:0;
        DefaultReadsMatcher::transferMatchingResults(approxMatcher);
        approxMatcher->matchedCountPerMismatches[0] = matchedReadsCount;
        approxMatcher->matchedCountPerMismatches[NOT_MATCHED_COUNT] = readsCount - matchedReadsCount;
    }

    void DefaultReadsMatcher::writeMatchesInfo(const string& outPrefix) {
        string offsetsFile = outPrefix + OFFSETS_SUFFIX;
        std::ofstream offsetsDest(offsetsFile, std::ios::out | std::ios::binary);
        if (offsetsDest.fail()) {
            fprintf(stderr, "cannot write to offsets file %s\n", offsetsFile.c_str());
            exit(EXIT_FAILURE);
        }
        string missedReadsFile = outPrefix + MISSED_READS_SUFFIX;
        std::ofstream missedReadsDest(missedReadsFile, std::ios::out | std::ios::binary);
        if (missedReadsDest.fail()) {
            fprintf(stderr, "cannot write to missed reads file %s\n", missedReadsFile.c_str());
            exit(EXIT_FAILURE);
        }
        string dumpFile = outPrefix + DUMP_SUFFIX;
        std::ofstream dumpDest;
        dumpDest.open(dumpFile, std::ios::out | std::ios::binary);
        if (dumpDest.fail()) {
            fprintf(stderr, "cannot write to dump file %s\n", dumpFile.c_str());
            exit(EXIT_FAILURE);
        }

        writeMatchesInfo(offsetsDest, missedReadsDest, dumpDest);

        offsetsDest.close();
        missedReadsDest.close();
        dumpDest.close();
    }
    void DefaultReadsMatcher::matchConstantLengthReads() {
        initMatching();

        this->executeMatching(false);

        if (revComplPg) {
            PgSAHelpers::reverseComplementInPlace(pgPtr, pgLength);
            this->executeMatching(true);
            PgSAHelpers::reverseComplementInPlace(pgPtr, pgLength);
        }
    }

    void AbstractReadsApproxMatcher::continueMatchingConstantLengthReads(DefaultReadsMatcher *pMatcher) {
        this->initMatchingContinuation(pMatcher);

        this->executeMatching(false);

        if (revComplPg) {
            PgSAHelpers::reverseComplementInPlace(pgPtr, pgLength);
            this->executeMatching(true);
            PgSAHelpers::reverseComplementInPlace(pgPtr, pgLength);
        }
    }

    DefaultReadsExactMatcher::DefaultReadsExactMatcher(char* pgPtr, const uint_pg_len_max pgLength, bool revComplPg,
                                                       ConstantLengthReadsSetInterface *readsSet, uint32_t matchPrefixLength)
            : DefaultReadsMatcher(pgPtr, pgLength, revComplPg, readsSet, matchPrefixLength) {}

    void DefaultReadsExactMatcher::initMatching() {
        *logout << "Feeding patterns...\n" << endl;
        this->hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(matchingLength);
        this->hashMatcher->addReadsSetOfPatterns(readsSet);
        *logout << "... checkpoint " << clock_millis() << " msec. " << endl;
        DefaultReadsMatcher::initMatching();
    }

    void DefaultReadsExactMatcher::executeMatching(bool revCompMode) {
        clock_checkpoint();
        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;
        hashMatcher->iterateOver(pgPtr, pgLength);

        while (hashMatcher->moveNext()) {
            const uint64_t matchPosition = hashMatcher->getHashMatchTextPosition();
            const uint_reads_cnt_max matchReadIndex = hashMatcher->getHashMatchPatternIndex();

            bool exactMatch = readsSet->compareReadWithPattern(matchReadIndex, pgPtr + matchPosition) == 0;
            if (exactMatch) {
                if (readMatchPos[matchReadIndex] == NOT_MATCHED_POSITION) {
                    readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                    if (revCompMode) readMatchRC[matchReadIndex] = true;
                    matchedReadsCount++;
                } else
                    betterMatchCount++;
            } else
                falseMatchCount++;
/*            if (i++ < 1) {
                cout << "Matched: " << matchReadIndex << "; " << matchPosition
                    << "; " << (exactMatch?"positive":"false") << "; " << (revCompMode?"pair strand (RC)":"") << endl;
                cout << matchedRead << endl;
                const string pgPart(txt + (matchPosition), matchingLength);
                cout << pgPart << endl;
            }*/
        }

        cout << "... exact matching procedure completed in " << clock_millis() << " msec. " << endl;
        cout << "Exact matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << betterMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;
    }

    AbstractReadsApproxMatcher::AbstractReadsApproxMatcher(char* pgPtr, const uint_pg_len_max pgLength, bool revComplPg,
                                                           ConstantLengthReadsSetInterface *readsSet, uint32_t matchPrefixLength,
                                                           uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches)
            : DefaultReadsMatcher(pgPtr, pgLength, revComplPg, readsSet, matchPrefixLength),
            targetMismatches(readsSet->maxReadLength() / readsExactMatchingChars - 1),
            maxMismatches(maxMismatches), minMismatches(minMismatches){
        currentRead.resize(readsSet->maxReadLength());
        matchedCountPerMismatches[NOT_MATCHED_COUNT] = readsCount;
    }

    InterleavedReadsApproxMatcher::InterleavedReadsApproxMatcher(char* pgPtr, const uint_pg_len_max pgLength, bool revComplPg,
                                                                 ConstantLengthReadsSetInterface *readsSet, uint32_t matchPrefixLength,
                                                uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches)
            : AbstractReadsApproxMatcher(pgPtr, pgLength, revComplPg, readsSet, matchPrefixLength, readsExactMatchingChars,
                    maxMismatches, minMismatches), partLength(readsExactMatchingChars)  {}

    DefaultReadsApproxMatcher::DefaultReadsApproxMatcher(char* pgPtr, const uint_pg_len_max pgLength, bool revComplPg,
                                                         ConstantLengthReadsSetInterface *readsSet, uint32_t matchPrefixLength,
                                                         uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches)
             :AbstractReadsApproxMatcher(pgPtr, pgLength, revComplPg, readsSet, matchPrefixLength, readsExactMatchingChars,
                     maxMismatches, minMismatches), partLength(readsExactMatchingChars) {}

    CopMEMReadsApproxMatcher::CopMEMReadsApproxMatcher(char* pgPtr, const uint_pg_len_max pgLength, bool revComplPg,
                                                       ConstantLengthReadsSetInterface *readsSet,
                                                       uint32_t matchPrefixLength, uint16_t readsExactMatchingChars,
                                                       uint8_t maxMismatches, uint8_t minMismatches)
            : AbstractReadsApproxMatcher(pgPtr, pgLength, revComplPg, readsSet, matchPrefixLength, readsExactMatchingChars,
                                         maxMismatches, minMismatches), partLength(readsExactMatchingChars) {

    }

    CopMEMReadsApproxMatcher::~CopMEMReadsApproxMatcher() {
    }


    void AbstractReadsApproxMatcher::printApproxMatchingStats() {
        cout << "Matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << betterMatchCount << " better-matches) in " << clock_millis() << " msec.  False matches reported: " << falseMatchCount << "."
             << endl;

        for (uint8_t i = 0; i <= maxMismatches; i++)
            *logout << "Matched " << matchedCountPerMismatches[i] << " reads with " << (int) i << " mismatches." << endl;
    }

    void DefaultReadsApproxMatcher::initMatching() {
        *logout << "Feeding patterns...\n" << endl;
        hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(partLength);
        this->hashMatcher->addReadsSetOfPatterns(readsSet, targetMismatches + 1);
        *logout << "... checkpoint " << clock_millis() << " msec. " << endl;

        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, NOT_MATCHED_COUNT);
    }

    void DefaultReadsApproxMatcher::initMatchingContinuation(DefaultReadsMatcher *pMatcher) {
        cout << "Feeding unmatched patterns...\n" << endl;
        hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(partLength);
        this->hashMatcher->addReadsSetOfPatterns(readsSet, targetMismatches + 1,
                                                 pMatcher->getMatchedReadsBitmap(minMismatches));
        *logout << "... checkpoint " << clock_millis() << " msec. " << endl;

        AbstractReadsApproxMatcher::initMatchingContinuation(pMatcher);
    }

    void DefaultReadsApproxMatcher::executeMatching(bool revCompMode) {
        clock_checkpoint();
        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;
        hashMatcher->iterateOver(pgPtr, pgLength);
        while (hashMatcher->moveNext()) {
            const uint32_t matchPatternIndex = hashMatcher->getHashMatchPatternIndex();
            uint32_t matchReadIndex = matchPatternIndex / (targetMismatches + 1);
            if (readMismatchesCount[matchReadIndex] <= minMismatches)
                continue;
            uint64_t matchPosition = hashMatcher->getHashMatchTextPosition();
            const uint_read_len_max positionShift = ((matchPatternIndex % (targetMismatches + 1)) * partLength);
            if (positionShift > matchPosition)
                continue;
            matchPosition -= positionShift;
            if (matchPosition + readLength > pgLength)
                continue;
            if (readMatchPos[matchReadIndex] == (revCompMode?pgLength-(matchPosition+matchingLength):matchPosition))
                continue;
            uint8_t currentMismatchesLimit = readMismatchesCount[matchReadIndex]==NOT_MATCHED_COUNT?maxMismatches
                    :(readMismatchesCount[matchReadIndex] - 1);
            const uint8_t mismatchesCount = readsSet->countMismatchesVsPattern(matchReadIndex, pgPtr + matchPosition,
                                                            matchingLength, currentMismatchesLimit);
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == NOT_MATCHED_COUNT)
                    matchedReadsCount++;
                 else
                    betterMatchCount++;
                matchedCountPerMismatches[readMismatchesCount[matchReadIndex]]--;
                matchedCountPerMismatches[mismatchesCount]++;
                readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                readMatchRC[matchReadIndex] = revCompMode;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            } else
                falseMatchCount++;

/*            if (mismatchesCount < UINT8_MAX && i++ < 2) {
                cout << "Matched: " << matchReadIndex << " (" << matchPatternIndex << "); "
                     << matchPosition << "; " << (int) mismatchesCount << "; " << (revCompMode?"pair strand (RC)":"") << endl;
                cout << matchedRead << endl;
                const string pgPart(txt + matchPosition, matchingLength);
                cout << pgPart << endl;
            }*/
        }
        this->printApproxMatchingStats();
    }

    void InterleavedReadsApproxMatcher::initMatching() {
        cout << "Feeding patterns...\n" << endl;
        hashMatcher = new InterleavedConstantLengthPatternsOnTextHashMatcher(partLength, targetMismatches + 1);
        this->hashMatcher->addPackedPatterns(readsSet, targetMismatches + 1);
        *logout << "... checkpoint " << clock_millis() << " msec. " << endl;

        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, NOT_MATCHED_COUNT);
    }

    void InterleavedReadsApproxMatcher::initMatchingContinuation(DefaultReadsMatcher *pMatcher) {
        cout << "Feeding unmatched patterns...\n" << endl;
        hashMatcher = new InterleavedConstantLengthPatternsOnTextHashMatcher(partLength, targetMismatches + 1);
        this->hashMatcher->addPackedPatterns(readsSet, targetMismatches + 1,
                pMatcher->getMatchedReadsBitmap(minMismatches));
        *logout << "... checkpoint " << clock_millis() << " msec. " << endl;

        AbstractReadsApproxMatcher::initMatchingContinuation(pMatcher);
    }

    void InterleavedReadsApproxMatcher::executeMatching(bool revCompMode) {
        clock_checkpoint();
        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;
        hashMatcher->iterateOver(pgPtr, pgLength);
        while (hashMatcher->moveNext()) {
            const uint32_t matchPatternIndex = hashMatcher->getHashMatchPatternIndex();
            uint32_t matchReadIndex = matchPatternIndex / (targetMismatches + 1);
            if (readMismatchesCount[matchReadIndex] <= minMismatches)
                continue;
            uint64_t matchPosition = hashMatcher->getHashMatchTextPosition();
            const uint8_t positionShift = matchPatternIndex % (targetMismatches + 1);
            if (positionShift > matchPosition)
                continue;
            matchPosition -= positionShift;
            if (matchPosition + readLength > pgLength)
                continue;
            if (readMatchPos[matchReadIndex] == (revCompMode?pgLength-(matchPosition+matchingLength):matchPosition))
                continue;
            uint8_t currentMatchesLimit = readMismatchesCount[matchReadIndex]==NOT_MATCHED_COUNT?maxMismatches
                    :(readMismatchesCount[matchReadIndex] - 1);
            const uint8_t mismatchesCount = readsSet->countMismatchesVsPattern(matchReadIndex, pgPtr + matchPosition,
                    matchingLength, currentMatchesLimit);
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == NOT_MATCHED_COUNT)
                    matchedReadsCount++;
                else
                    betterMatchCount++;
                matchedCountPerMismatches[readMismatchesCount[matchReadIndex]]--;
                matchedCountPerMismatches[mismatchesCount]++;
                readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                readMatchRC[matchReadIndex] = revCompMode;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            } else
                falseMatchCount++;

            //if (mismatchesCount < UINT8_MAX && i++ < 2) {
/*            if (mismatchesCount < UINT8_MAX && mismatchesCount > 8) {
                cout << "Matched: " << matchReadIndex << " (" << matchPatternIndex << "); "
                     << matchPosition << "; " << (int) mismatchesCount << "; " << (revCompMode?"pair strand (RC)":"") << endl;
                cout << matchedRead << endl;
                const string pgPart(txt + matchPosition, matchingLength);
                cout << pgPart << endl;
            }*/
        }
        this->printApproxMatchingStats();
    }

    void CopMEMReadsApproxMatcher::initMatching() {
        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, NOT_MATCHED_COUNT);
    }

    void CopMEMReadsApproxMatcher::initMatchingContinuation(DefaultReadsMatcher *pMatcher) {
        AbstractReadsApproxMatcher::initMatchingContinuation(pMatcher);
    }

    void CopMEMReadsApproxMatcher::executeMatching(bool revCompMode) {
        clock_checkpoint();
        cout << "Feeding " << (revCompMode?"rc of ":"") << "pseudogenome sequence... " << endl;
        CopMEMMatcher* copMEMMatcher = new CopMEMMatcher(pgPtr, pgLength, partLength);
        *logout << "... checkpoint " << clock_millis() << " msec. " << endl;
        uint_reads_cnt_std count = 0;
        for(uint_reads_cnt_max matchReadIndex = 0; matchReadIndex < readsCount; matchReadIndex++) {
            char_pg currentReadPtr[UINT8_MAX];
            if (readMismatchesCount[matchReadIndex] <= minMismatches)
                continue;
            readsSet->getRead(matchReadIndex, currentReadPtr);
            uint8_t mismatchesCount = readMismatchesCount[matchReadIndex];
            uint64_t matchPosition = copMEMMatcher->approxMatchPattern(currentReadPtr, matchingLength,
                                                                       maxMismatches, minMismatches, mismatchesCount,
                                                                       betterMatchCount, falseMatchCount);
            if (matchPosition == UINT64_MAX)
                continue;
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == NOT_MATCHED_COUNT)
                    count++;
                matchedCountPerMismatches[readMismatchesCount[matchReadIndex]]--; // non-thread safe
                matchedCountPerMismatches[mismatchesCount]++; // non-thread safe
                readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                readMatchRC[matchReadIndex] = revCompMode;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            }
        }
        matchedReadsCount += count;
        delete(copMEMMatcher);
        printApproxMatchingStats();
    }

    const string DefaultReadsMatcher::OFFSETS_SUFFIX = "_matched_offsets.txt";
    const string DefaultReadsMatcher::DUMP_SUFFIX = "_matches_info.txt";
    const string DefaultReadsMatcher::MISSED_READS_SUFFIX = "_missed.txt";

    void DefaultReadsExactMatcher::writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &dumpDest) {
        clock_checkpoint();

        for (uint_reads_cnt_max i = 0; i < readMatchPos.size(); i++) {
            if (readMatchPos[i] == NOT_MATCHED_POSITION)
                missedPatternsDest << readsSet->getRead(i) << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << (readMatchRC[i]?"\tRC":"") << "\n";
                if (matchingLength < readLength)
                    dumpDest << readsSet->getRead(i).substr(matchingLength);
            }
        }

        cout << "... writing info dump files completed in  " << clock_millis() << " msec. " << endl;
    }

    void AbstractReadsApproxMatcher::writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &dumpDest) {
        clock_checkpoint();

        vector<uint_reads_cnt_max> mismatchedReadsCount(maxMismatches + 1, 0);
        vector<uint64_t> mismatchesPerReadPositionCount(readLength, 0);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            if (readMatchPos[i] == NOT_MATCHED_POSITION)
                missedPatternsDest << readsSet->getRead(i) << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << (readMatchRC[i]?"\tRC":"");
                const string read = readMatchRC[i]?reverseComplement(readsSet->getRead(i)):readsSet->getRead(i);
                reportMismatches(read.data(), pgPtr + readMatchPos[i], matchingLength, offsetsDest);
                offsetsDest << "\n";
                uint_read_len_max mismatchesCount = 0;
                if (matchingLength < readLength) {
                    dumpDest << readsSet->getRead(i).substr(matchingLength);
                } else {
                    dumpDest << i << "\t" << readMatchPos[i] << (readMatchRC[i]?"\tRC":"") << endl;
                    dumpDest << read << endl;
                    const string pgPart(pgPtr + readMatchPos[i], readLength);
                    dumpDest << pgPart << endl;
                    for(int j = 0; j < read.length(); j++) {
                        const bool symbolMatch = read[j] == pgPart[j];
                        dumpDest.put(symbolMatch ? '-' : 'x');
                        if (!symbolMatch) {
                            mismatchesCount++;
                            if (readMatchRC[i])
                                mismatchesPerReadPositionCount[readLength - j - 1]++;
                            else
                                mismatchesPerReadPositionCount[j]++;
                        }
                    }
                    dumpDest << endl;
                }
                mismatchedReadsCount[mismatchesCount]++;
            }
        }
        offsetsDest << endl << "Stats:" << endl;
        for (uint8_t i = 0; i <= maxMismatches; i++)
            offsetsDest << "Matched " << mismatchedReadsCount[i] << " reads with " << (int) i << " mismatches." << endl;

        offsetsDest << endl << "Mismatches per reads position:" << endl;
        for (uint8_t i = 0; i < readLength; i++)
            offsetsDest << (int) i << ":\t" << mismatchesPerReadPositionCount[i] << endl;

        cout << endl << "... writing info dump files completed in " << clock_millis() << " msec. " << endl;
    }

    SeparatedPseudoGenomeOutputBuilder *AbstractReadsApproxMatcher::createSeparatedPseudoGenomeOutputBuilder(
            SeparatedPseudoGenome* sPg) {
        DefaultReadsListIteratorInterface* rlIt = sPg->getReadsList();
        bool isRevCompEnabled = sPg->getReadsList()->isRevCompEnabled();
        bool areMismatchesEnabled = sPg->getReadsList()->areMismatchesEnabled();
        SeparatedPseudoGenomeOutputBuilder* builder = new SeparatedPseudoGenomeOutputBuilder(
                !isRevCompEnabled && !this->revComplPg, !areMismatchesEnabled && this->maxMismatches == 0);
        builder->setReadsSourceIterator(rlIt);
        builder->copyPseudoGenomeProperties(sPg);
        return builder;
    }

    SeparatedPseudoGenomeOutputBuilder *DefaultReadsExactMatcher::createSeparatedPseudoGenomeOutputBuilder(
            SeparatedPseudoGenome *sPg){
        DefaultReadsListIteratorInterface* rlIt = sPg->getReadsList();
        bool isRevCompEnabled = sPg->getReadsList()->isRevCompEnabled();
        bool areMismatchesEnabled = sPg->getReadsList()->areMismatchesEnabled();
        SeparatedPseudoGenomeOutputBuilder* builder = new SeparatedPseudoGenomeOutputBuilder(
                !isRevCompEnabled && !this->revComplPg, !areMismatchesEnabled);
        builder->setReadsSourceIterator(rlIt);
        builder->copyPseudoGenomeProperties(sPg);
        return builder;
    }

    void AbstractReadsApproxMatcher::initEntryUpdating() { }

    void AbstractReadsApproxMatcher::updateEntry(DefaultReadsListEntry &entry, uint_reads_cnt_max matchIdx,
            bool revComplPairFile) {
        readsSet->getRead(matchIdx, (char_pg*) currentRead.data());
        if (readMatchRC[matchIdx])
            reverseComplementInPlace(currentRead);
        if (revComplPairFile?(readMatchRC[matchIdx] != (bool) (entry.idx % 2)):readMatchRC[matchIdx])
            fillEntryWithReversedMismatches(currentRead.data(), pgPtr + entry.pos, readMismatchesCount[matchIdx], entry,
                    readLength);
        else
            fillEntryWithMismatches(currentRead.data(), pgPtr + entry.pos, readMismatchesCount[matchIdx], entry);

    }

    void AbstractReadsApproxMatcher::closeEntryUpdating() { }

    void DefaultReadsMatcher::exportMatchesInPgOrder(SeparatedPseudoGenome* sPg, ostream &pgrcOut, uint8_t compressionLevel,
                                                     const string &outPgPrefix, IndexesMapping *orgIndexesMapping,
                                                     bool pairFileMode, bool revComplPairFile) {
        clock_checkpoint();
        vector<uint_reads_cnt_max> idxs(matchedReadsCount);
        uint64_t counter = 0;
        for(uint_reads_cnt_max i = 0; i < readsCount; i++)
            if (readMatchPos[i] != NOT_MATCHED_POSITION)
                idxs[counter++] = i;

        std::sort(idxs.begin(), idxs.end(), [this](const uint_reads_cnt_max& idx1, const uint_reads_cnt_max& idx2) -> bool
        { return readMatchPos[idx1] < readMatchPos[idx2]; });

        initEntryUpdating();
        SeparatedPseudoGenomeOutputBuilder* builder = this->createSeparatedPseudoGenomeOutputBuilder(sPg);
        for(uint_reads_cnt_max i = 0; i < matchedReadsCount; i++) {
            uint_reads_cnt_max matchIdx = idxs[i];
            uint64_t currPos = builder->writeReadsFromIterator(readMatchPos[matchIdx]);
            DefaultReadsListEntry entry(currPos);
            const uint_reads_cnt_max orgIdx = orgIndexesMapping->getReadOriginalIndex(matchIdx);
            entry.advanceEntryByPosition(readMatchPos[matchIdx], orgIdx, readMatchRC[matchIdx]);
            this->updateEntry(entry, matchIdx, revComplPairFile);
            builder->writeExtraReadEntry(entry);
        }
        builder->writeReadsFromIterator();
        builder->build(outPgPrefix);
        builder->compressedBuild(pgrcOut, compressionLevel);
        if (pairFileMode)
            builder->updateOriginalIndexesIn(sPg);
        delete(builder);
        closeEntryUpdating();
        *logout << "... exporting matches (" << outPgPrefix << ") completed in " << clock_millis() << " msec. " << endl << endl;
    }

    void DefaultReadsMatcher::exportMatchesInOriginalOrder(SeparatedPseudoGenome *sPg, ostream &pgrcOut,
                                                           uint8_t compressionLevel, const string &outPgPrefix,
                                                           IndexesMapping *orgIndexesMapping, bool pairFileMode,
                                                           bool revComplPairFile) {
        clock_checkpoint();

        uint_reads_cnt_std readsTotalCount = orgIndexesMapping->getReadsTotalCount();
        vector<uint_pg_len_max> orgIdx2pgPos(readsTotalCount, -1);
        ExtendedReadsListWithConstantAccessOption *const pgRl = sPg->getReadsList();
        uint_pg_len_max pos = 0;
        for(uint_reads_cnt_std i = 0; i < pgRl->readsCount; i++) {
            pos += pgRl->off[i];
            orgIdx2pgPos[pgRl->orgIdx[i]] = pos;
        }
        pgRl->off.clear();
        pgRl->orgIdx.clear();
        initEntryUpdating();
        SeparatedPseudoGenomeOutputBuilder* builder = this->createSeparatedPseudoGenomeOutputBuilder(sPg);
        uint_reads_cnt_max nI_start = readsCount;
        int64_t curOrgIdx = 0;
        for(uint_reads_cnt_std i = 0; i < readsCount; i++) {
            const uint_reads_cnt_max oIdx = orgIndexesMapping->getReadOriginalIndex(i);
            if (curOrgIdx > oIdx) {
                nI_start = i;
                break;
            }
            curOrgIdx = oIdx;
        }
        const uint8_t parts = pairFileMode?2:1;
        const int inc = parts;
        for(uint8_t p = 0; p < parts; p++) {
            curOrgIdx = p - inc;
            uint_reads_cnt_max lqI = 0;
            uint_reads_cnt_max nI = nI_start;
            while (lqI < nI_start || nI < readsCount) {
                uint_reads_cnt_max matchIdx;
                uint_reads_cnt_max oIdx;
                do {
                    matchIdx = readsCount;
                    if (lqI < nI_start)
                        matchIdx = lqI;
                    if (nI < readsCount && (lqI == nI_start || orgIndexesMapping->getReadOriginalIndex(nI) <
                                                               orgIndexesMapping->getReadOriginalIndex(lqI)))
                        matchIdx = nI++;
                    else
                        lqI++;
                    if (matchIdx == readsCount)
                        break;
                    oIdx = orgIndexesMapping->getReadOriginalIndex(matchIdx);
                } while (parts != 1 && (oIdx % parts != p));
                if (matchIdx == readsCount)
                    break;
                while ((curOrgIdx += inc) < oIdx) {
                    DefaultReadsListEntry entry(0);
                    entry.advanceEntryByPosition(0, curOrgIdx, false);
                    builder->writeExtraReadEntry(entry);
                }
                if (readMatchPos[matchIdx] != NOT_MATCHED_POSITION) {
                    DefaultReadsListEntry entry(0);
                    entry.advanceEntryByPosition(readMatchPos[matchIdx], oIdx, readMatchRC[matchIdx]);
                    this->updateEntry(entry, matchIdx, revComplPairFile);
                    builder->writeExtraReadEntry(entry);
                    orgIdx2pgPos[oIdx] = readMatchPos[matchIdx];
                }
            }
            while ((curOrgIdx += inc) < readsTotalCount) {
                DefaultReadsListEntry entry(0);
                entry.advanceEntryByPosition(0, curOrgIdx, false);
                builder->writeExtraReadEntry(entry);
            }
        }
        builder->build(outPgPrefix);
        builder->compressedBuild(pgrcOut, compressionLevel, true);
        delete(builder);
        closeEntryUpdating();

        sPg->getReadsList()->pos = std::move(orgIdx2pgPos);
        *logout << "... exporting matches (" << outPgPrefix << ") completed in " << clock_millis() << " msec. " << endl << endl;
    }

    const vector<bool> DefaultReadsMatcher::getMatchedReadsBitmap(uint8_t maxMismatches) {
        vector<bool> res;
        res.reserve(readsCount);
        for(const uint64_t matchPos: readMatchPos)
            res.push_back(matchPos != DefaultReadsMatcher::NOT_MATCHED_POSITION);
        return res;
    }

    const vector<bool> AbstractReadsApproxMatcher::getMatchedReadsBitmap(uint8_t maxMismatches) {
        vector<bool> res;
        res.reserve(readsCount);
        for(const uint8_t mismatchesCount: readMismatchesCount)
            res.push_back(mismatchesCount <= maxMismatches);
        return res;
    }

    const vector<bool> mapReadsIntoPg(SeparatedPseudoGenome* sPg, bool revComplPg, bool preserveOrderMode,
                        ConstantLengthReadsSetInterface *readsSet, bool pairFileMode, bool revComplPairFile,
                        uint_read_len_max matchPrefixLength, uint16_t preReadsExactMatchingChars,
                        uint16_t readsExactMatchingChars, uint16_t minCharsPerMismatch, char preMatchingMode,
                        char matchingMode, bool dumpInfo, ostream &pgrcOut, uint8_t compressionLevel,
                        const string &pgDestFilePrefix, IndexesMapping* orgIndexesMapping) {
        uint_read_len_max readLength = readsSet->maxReadLength();
        uint8_t maxMismatches = readLength / minCharsPerMismatch;
        if (readsExactMatchingChars > readLength)
            readsExactMatchingChars = readLength;
        if (preReadsExactMatchingChars > readLength)
            preReadsExactMatchingChars = readLength;
        uint16_t currentExactMatchingChars = readsExactMatchingChars;
        char currentMatchingMode = matchingMode;
        if (preReadsExactMatchingChars > 0) {
            currentExactMatchingChars = preReadsExactMatchingChars;
            currentMatchingMode = preMatchingMode;
        }
        bool shortcutMode = toupper(currentMatchingMode) == currentMatchingMode;
        uint8_t currentMinMismatches = shortcutMode?maxMismatches:0;
        uint8_t targetMismatches = readLength / currentExactMatchingChars - 1;
        DefaultReadsMatcher* matcher;
        if (readLength == currentExactMatchingChars)
            switch (tolower(currentMatchingMode)) {
                case 'c': matcher = new CopMEMReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                        currentExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                default: matcher = new DefaultReadsExactMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength);
            }
        else switch (tolower(currentMatchingMode)) {
                case 'd': matcher = new DefaultReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                        currentExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                case 'i': matcher = new InterleavedReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                        currentExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                case 'c': matcher = new CopMEMReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                        currentExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                default:
                    fprintf(stderr, "Unknown matching mode: %c.\n", currentMatchingMode);
                    exit(EXIT_FAILURE);
            }
        cout << "Target pseudogenome length: " << sPg->getPgSequence().length() << endl;
        *logout << endl;
        cout << "readsAlignmentSeedLength (minCharsPerMismatch, matchingMode): " << (int) currentExactMatchingChars <<
             " (" << (int) minCharsPerMismatch << ", " << currentMatchingMode << ")" << endl;
        *logout << "targetMismatches (maxMismatches, minMismatches): " << (int) targetMismatches <<
             " (" << (int) maxMismatches << ", " << (int) currentMinMismatches << ")" << endl;
        matcher->matchConstantLengthReads();

        if (preReadsExactMatchingChars > 0) {
            AbstractReadsApproxMatcher* approxMatcher;
            bool shortcutMode = toupper(matchingMode) == matchingMode;
            uint8_t currentMinMismatches = shortcutMode?maxMismatches:targetMismatches + 1;
            switch (tolower(matchingMode)) {
                case 'd': approxMatcher = new DefaultReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                    sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                    readsExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                case 'i': approxMatcher = new InterleavedReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                        readsExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                case 'c': approxMatcher = new CopMEMReadsApproxMatcher((char*) sPg->getPgSequence().data(),
                        sPg->getPgSequence().length(), revComplPg, readsSet, matchPrefixLength,
                        readsExactMatchingChars, maxMismatches, currentMinMismatches);
                    break;
                default:
                    fprintf(stderr, "Unknown mismatches mode: %c.\n", matchingMode);
                    exit(EXIT_FAILURE);
            }
            targetMismatches = readLength / readsExactMatchingChars - 1;
            cout << endl << "Reads matching 2nd PHASE." << endl;
            cout << "readsExactMatchingChars (minCharsPerMismatch, matchingMode): " << (int) readsExactMatchingChars <<
                 " (" << (int) minCharsPerMismatch << ", " << matchingMode << ")" << endl;
            cout << "targetMismatches (maxMismatches, minMismatches): " << (int) targetMismatches <<
                 " (" << (int) maxMismatches << ", " << (int) currentMinMismatches << ")" << endl;
            approxMatcher->continueMatchingConstantLengthReads(matcher);
            delete(matcher);
            matcher = approxMatcher;
        }
        if (dumpInfo)
            matcher->writeMatchesInfo(pgDestFilePrefix);

        const vector<bool> res = matcher->getMatchedReadsBitmap();

        if (matchPrefixLength == DefaultReadsMatcher::DISABLED_PREFIX_MODE) {
            if (preserveOrderMode)
                matcher->exportMatchesInOriginalOrder(sPg, pgrcOut, compressionLevel, pgDestFilePrefix, orgIndexesMapping,
                                                      pairFileMode, revComplPairFile);
            else
                matcher->exportMatchesInPgOrder(sPg, pgrcOut, compressionLevel, pgDestFilePrefix, orgIndexesMapping,
                                                pairFileMode, revComplPairFile);
        }

        delete(matcher);
        return res;
    }

    uint8_t matchingCharsCorrection(size_t pgLength) {
        int x = pgLength /= 10000000;
        return x?32 - __builtin_clz(x):0;
    }

}