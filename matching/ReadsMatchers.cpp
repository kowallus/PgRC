#include "ReadsMatchers.h"

#include "../readsset/PackedConstantLengthReadsSet.h"
#include "../readsset/persistance/ReadsSetPersistence.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../pseudogenome/readslist/SeparatedExtendedReadsList.h"

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

    const uint_read_len_max DefaultReadsMatcher::DISABLED_PREFIX_MODE = (uint_read_len_max) -1;
    const uint64_t DefaultReadsMatcher::NOT_MATCHED_POSITION = UINT64_MAX;

    DefaultReadsMatcher::DefaultReadsMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                                             uint32_t matchPrefixLength) :
                                             sPg(sPg), pgPtr(sPg->getPgSequence().data()),
                                             pgLength(sPg->getPgSequence().length()), revComplPg(revComplPg),
                                             readsSet(readsSet), matchPrefixLength(matchPrefixLength),
                                             readLength(readsSet->minReadLength()),
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
        multiMatchCount = 0;
        falseMatchCount = 0;
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

        string suffixesFile = outPrefix + SUFFIXES_SUFFIX;
        std::ofstream suffixesDest;
        if (matchPrefixLength != DISABLED_PREFIX_MODE) {
            suffixesDest.open(suffixesFile, std::ios::out | std::ios::binary);
            if (suffixesDest.fail()) {
                fprintf(stderr, "cannot write to suffixes file %s\n", suffixesFile.c_str());
                exit(EXIT_FAILURE);
            }
        }

        writeMatchesInfo(offsetsDest, missedReadsDest, suffixesDest);

        offsetsDest.close();
        missedReadsDest.close();
    }
    void DefaultReadsMatcher::matchConstantLengthReads() {
        clock_checkpoint();
        initMatching();

        this->executeMatching(false);

        if (revComplPg) {
            PgSAHelpers::reverseComplementInPlace(sPg->getPgSequence());
            this->executeMatching(true);
            PgSAHelpers::reverseComplementInPlace(sPg->getPgSequence());
        }
    }

    DefaultReadsExactMatcher::DefaultReadsExactMatcher(SeparatedPseudoGenome* sPg, bool revComplPg,
                                                       PackedConstantLengthReadsSet *readsSet, uint32_t matchPrefixLength)
            : DefaultReadsMatcher(sPg, revComplPg, readsSet, matchPrefixLength) {}

    void DefaultReadsExactMatcher::initMatching() {
        cout << "Feeding patterns...\n" << endl;
        this->hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(matchingLength);
        this->hashMatcher->addPackedPatterns(readsSet);
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        DefaultReadsMatcher::initMatching();
    }

    void DefaultReadsExactMatcher::executeMatching(bool revCompMode) {
        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;

        hashMatcher->iterateOver(pgPtr, pgLength);

        while (hashMatcher->moveNext()) {
            const uint64_t matchPosition = hashMatcher->getHashMatchTextPosition();
            const uint_reads_cnt_max matchReadIndex = hashMatcher->getHashMatchPatternIndex();

            bool exactMatch = readsSet->comparePackedReadWithPattern(matchReadIndex, pgPtr + matchPosition) == 0;
            if (exactMatch) {
                if (readMatchPos[matchReadIndex] == NOT_MATCHED_POSITION) {
                    readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                    if (revCompMode) readMatchRC[matchReadIndex] = true;
                    matchedReadsCount++;
                } else
                    multiMatchCount++;
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
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;
    }

    AbstractReadsApproxMatcher::AbstractReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg,
                                                           PackedConstantLengthReadsSet *readsSet, uint32_t matchPrefixLength,
                                                           uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches)
            : DefaultReadsMatcher(sPg, revComplPg, readsSet, matchPrefixLength),
            targetMismatches(readsSet->maxReadLength() / readsExactMatchingChars - 1),
            maxMismatches(maxMismatches), minMismatches(minMismatches){
        currentRead.resize(readsSet->maxReadLength());
        matchedCountPerMismatches[NOT_MATCHED_COUNT] = readsCount;
    }

    InterleavedReadsApproxMatcher::InterleavedReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg,
                                                                 PackedConstantLengthReadsSet *readsSet, uint32_t matchPrefixLength,
                                                uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches)
            : AbstractReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength, readsExactMatchingChars,
                    maxMismatches, minMismatches), partLength(matchingLength / (targetMismatches + 1))  {}

    DefaultReadsApproxMatcher::DefaultReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg,
                                                         PackedConstantLengthReadsSet *readsSet, uint32_t matchPrefixLength,
                                                         uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches)
             :AbstractReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength, readsExactMatchingChars,
                     maxMismatches, minMismatches), partLength(matchingLength / (targetMismatches + 1)) {}

    void AbstractReadsApproxMatcher::printApproxMatchingStats() {
        cout << "... approximate matching procedure completed in " << clock_millis() << " msec. " << endl;
        cout << "Matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;

        for (uint8_t i = 0; i <= maxMismatches; i++)
            cout << "Matched " << matchedCountPerMismatches[i] << " reads with " << (int) i << " mismatches." << endl;
    }

    void DefaultReadsApproxMatcher::initMatching() {
        cout << "Feeding patterns...\n" << endl;
        hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(partLength);
        this->hashMatcher->addPackedPatterns(readsSet, targetMismatches + 1);
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;

        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, NOT_MATCHED_COUNT);
    }

    void DefaultReadsApproxMatcher::executeMatching(bool revCompMode) {

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
            uint8_t currentMatchesLimit = readMismatchesCount[matchReadIndex]==NOT_MATCHED_COUNT?maxMismatches
                    :(readMismatchesCount[matchReadIndex] - 1);
            const uint8_t mismatchesCount = readsSet->countMismatchesVsPattern(matchReadIndex, pgPtr + matchPosition,
                                                            matchingLength, currentMatchesLimit);
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == NOT_MATCHED_COUNT)
                    matchedReadsCount++;
                 else
                    multiMatchCount++;
                matchedCountPerMismatches[readMismatchesCount[matchReadIndex]]--;
                matchedCountPerMismatches[mismatchesCount]++;
                readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                if (revCompMode) readMatchRC[matchReadIndex] = true;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            } else if (mismatchesCount == NOT_MATCHED_COUNT)
                falseMatchCount++;
            else
                multiMatchCount++;
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
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;

        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, NOT_MATCHED_COUNT);
    }

    void InterleavedReadsApproxMatcher::executeMatching(bool revCompMode) {

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
                    multiMatchCount++;
                matchedCountPerMismatches[readMismatchesCount[matchReadIndex]]--;
                matchedCountPerMismatches[mismatchesCount]++;
                readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                if (revCompMode) readMatchRC[matchReadIndex] = true;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            } else if (mismatchesCount == NOT_MATCHED_COUNT)
                falseMatchCount++;
            else
                multiMatchCount++;
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

    CopMEMReadsApproxMatcher::CopMEMReadsApproxMatcher(SeparatedPseudoGenome *sPg, bool revComplPg,
                                                       PackedConstantLengthReadsSet *readsSet,
                                                       uint32_t matchPrefixLength, uint16_t readsExactMatchingChars,
                                                       uint8_t maxMismatches, uint8_t minMismatches)
            : AbstractReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength, readsExactMatchingChars,
                    maxMismatches, minMismatches), partLength(readsExactMatchingChars) {

    }

    CopMEMReadsApproxMatcher::~CopMEMReadsApproxMatcher() {
    }

    void CopMEMReadsApproxMatcher::initMatching() {
        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, NOT_MATCHED_COUNT);
    }

    void CopMEMReadsApproxMatcher::executeMatching(bool revCompMode) {
        clock_t ref_start = clock();
        CopMEMMatcher* copMEMMatcher = new CopMEMMatcher(sPg->getPgSequence(), partLength);
        cout << "Feeding " << (revCompMode?"rc of ":"") << "reference pseudogenome finished in " << clock_millis(ref_start) << " msec. " << endl;
        for(uint_reads_cnt_max matchReadIndex = 0; matchReadIndex < readsCount; matchReadIndex++) {
            if (readMismatchesCount[matchReadIndex] <= minMismatches)
                continue;
            readsSet->getRead(matchReadIndex, (char_pg*) currentRead.data());
            uint8_t mismatchesCount = readMismatchesCount[matchReadIndex];
            uint64_t matchPosition = copMEMMatcher->approxMatchPattern(currentRead.data(), matchingLength,
                                                                       maxMismatches, minMismatches, mismatchesCount,
                                                                       multiMatchCount, falseMatchCount);
            if (matchPosition == UINT64_MAX)
                continue;
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == NOT_MATCHED_COUNT)
                    matchedReadsCount++;
                matchedCountPerMismatches[readMismatchesCount[matchReadIndex]]--;
                matchedCountPerMismatches[mismatchesCount]++;
                readMatchPos[matchReadIndex] = revCompMode?pgLength-(matchPosition+matchingLength):matchPosition;
                if (revCompMode) readMatchRC[matchReadIndex] = true;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            }
        }
        delete(copMEMMatcher);
        printApproxMatchingStats();
    }

    const string DefaultReadsMatcher::OFFSETS_SUFFIX = "_matched_offsets.txt";
    const string DefaultReadsMatcher::SUFFIXES_SUFFIX = "_matched_suffixes.txt";
    const string DefaultReadsMatcher::MISSED_READS_SUFFIX = "_missed.txt";

    void DefaultReadsExactMatcher::writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest) {
        clock_checkpoint();

        for (uint_reads_cnt_max i = 0; i < readMatchPos.size(); i++) {
            if (readMatchPos[i] == NOT_MATCHED_POSITION)
                missedPatternsDest << readsSet->getRead(i) << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << (readMatchRC[i]?"\tRC":"") << "\n";
                if (matchingLength < readLength)
                    suffixesDest << readsSet->getRead(i).substr(matchingLength);
            }
        }

        cout << "... writing info dump files completed in  " << clock_millis() << " msec. " << endl;
    }

    void AbstractReadsApproxMatcher::writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest) {
        clock_checkpoint();

        const string &text = sPg->getPgSequence();
        vector<uint_reads_cnt_max> mismatchedReadsCount(maxMismatches + 1, 0);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            if (readMatchPos[i] == NOT_MATCHED_POSITION)
                missedPatternsDest << readsSet->getRead(i) << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << (readMatchRC[i]?"\tRC":"");
                const string read = readMatchRC[i]?reverseComplement(readsSet->getRead(i)):readsSet->getRead(i);
                reportMismatches(read.data(), text.data() + readMatchPos[i], matchingLength, offsetsDest);
                offsetsDest << "\n";
                if (matchingLength < readLength)
                    suffixesDest << readsSet->getRead(i).substr(matchingLength);
                mismatchedReadsCount[readMismatchesCount[i]]++;
            }
        }

        cout << "... writing info dump files completed in " << clock_millis() << " msec. " << endl;
    }

    const vector<uint_reads_cnt_max> DefaultReadsMatcher::getMatchedReadsIndexes() const {
        vector<uint_reads_cnt_max> matchedReads(matchedReadsCount);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            if (readMatchPos[i] != NOT_MATCHED_POSITION)
                matchedReads.push_back(i);
        }
        return matchedReads;
    }

    const vector<uint64_t> &DefaultReadsMatcher::getReadMatchPos() const {
        return readMatchPos;
    }

    const vector<uint8_t> &DefaultReadsMatcher::getReadMismatches() const {
        return readMismatchesCount;
    }

    SeparatedPseudoGenomeOutputBuilder *AbstractReadsApproxMatcher::createSeparatedPseudoGenomeOutputBuilder(
            const string &outPgPrefix, bool enableRevComp, bool enableMismatches) {
        return new SeparatedPseudoGenomeOutputBuilder(outPgPrefix,
                !enableRevComp && !this->revComplPg, !enableMismatches && this->targetMismatches == 0);
    }

    SeparatedPseudoGenomeOutputBuilder *DefaultReadsExactMatcher::createSeparatedPseudoGenomeOutputBuilder(
            const string &outPgPrefix, bool enableRevComp, bool enableMismatches){
        return new SeparatedPseudoGenomeOutputBuilder(outPgPrefix,
                !enableRevComp && !this->revComplPg, !enableMismatches);
    }

    void AbstractReadsApproxMatcher::initEntryUpdating() { }

    void AbstractReadsApproxMatcher::updateEntry(DefaultReadsListEntry &entry, uint_reads_cnt_max matchIdx) {
        readsSet->getRead(matchIdx, (char_pg*) currentRead.data());
        if (readMatchRC[matchIdx])
            reverseComplementInPlace(currentRead);
        fillEntryWithMismatches(currentRead.data(), pgPtr + entry.pos, readMismatchesCount[matchIdx], entry);
    }

    void AbstractReadsApproxMatcher::closeEntryUpdating() { }

    void DefaultReadsMatcher::writeIntoPseudoGenome(const string &outPgPrefix, IndexesMapping* orgIndexesMapping) {
        clock_checkpoint();
        vector<uint_reads_cnt_max> idxs(matchedReadsCount);
        uint64_t counter = 0;
        for(uint_reads_cnt_max i = 0; i < readsCount; i++)
            if (readMatchPos[i] != NOT_MATCHED_POSITION)
                idxs[counter++] = i;

        std::sort(idxs.begin(), idxs.end(), [this](const uint_reads_cnt_max& idx1, const uint_reads_cnt_max& idx2) -> bool
        { return readMatchPos[idx1] < readMatchPos[idx2]; });

        initEntryUpdating();
        DefaultReadsListIteratorInterface* rlIt = sPg->getReadsList();
        bool isRevCompEnabled = sPg->getReadsList()->isRevCompEnabled();
        bool areMismatchesEnabled = sPg->getReadsList()->areMismatchesEnabled();
        SeparatedPseudoGenomeOutputBuilder* builder = this->createSeparatedPseudoGenomeOutputBuilder(outPgPrefix,
                isRevCompEnabled, areMismatchesEnabled);
        builder->setReadsSourceIterator(rlIt);
        builder->copyPseudoGenomeProperties(sPg);
        for(uint_reads_cnt_max i = 0; i < matchedReadsCount; i++) {
            uint_reads_cnt_max matchIdx = idxs[i];
            uint64_t currPos = builder->writeReadsFromIterator(readMatchPos[matchIdx]);
            DefaultReadsListEntry entry(currPos);
            entry.advanceEntryByPosition(readMatchPos[matchIdx], orgIndexesMapping->getReadOriginalIndex(matchIdx), readMatchRC[matchIdx]);
            this->updateEntry(entry, matchIdx);
            builder->writeExtraReadEntry(entry);
        }
        builder->writeReadsFromIterator();
        builder->build();
        delete(builder);
        closeEntryUpdating();
        cout << "... writing (" << outPgPrefix << ") output files completed in " << clock_millis() << " msec. " << endl << endl;
    }

    const vector<bool> DefaultReadsMatcher::getMatchedReadsBitmap() {
        vector<bool> res;
        res.reserve(readsCount);
        for(const uint64_t matchPos: readMatchPos)
            res.push_back(matchPos != DefaultReadsMatcher::NOT_MATCHED_POSITION);
        return res;
    }

    const vector<bool> mapReadsIntoPg(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                        uint_read_len_max matchPrefixLength, uint16_t readsExactMatchingChars, uint16_t minCharsPerMismatch,
                        char mismatchesMode, uint8_t minMismatches, bool dumpInfo, const string &pgDestFilePrefix,
                        IndexesMapping* orgIndexesMapping) {
        DefaultReadsMatcher* matcher;
        uint_read_len_max readLength = readsSet->maxReadLength();
        uint8_t maxMismatches = readLength / minCharsPerMismatch;
        if (readsExactMatchingChars > readLength)
            readsExactMatchingChars = readLength;
        uint8_t targetMismatches = readLength / readsExactMatchingChars - 1;
        if (targetMismatches == 0)
            switch (mismatchesMode) {
                case 'c': matcher = new CopMEMReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength,
                                                                 readsExactMatchingChars, 0);
                    break;
                default: matcher = new DefaultReadsExactMatcher(sPg, revComplPg, readsSet, matchPrefixLength);
            }
        else switch (mismatchesMode) {
                case 'd': matcher = new DefaultReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength,
                                                                  readsExactMatchingChars, maxMismatches, minMismatches);
                    break;
                case 'i': matcher = new InterleavedReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength,
                                                                      readsExactMatchingChars, maxMismatches, minMismatches);
                    break;
                case 'c': matcher = new CopMEMReadsApproxMatcher(sPg, revComplPg, readsSet, matchPrefixLength,
                                                                      readsExactMatchingChars, maxMismatches, minMismatches);
                    break;
                default:
                    fprintf(stderr, "Unknown mismatches mode: %c.\n", mismatchesMode);
                    exit(EXIT_FAILURE);

            }

        cout << "Target pseudogenome length: " << sPg->getPgSequence().length() << endl << endl;
        cout << "readsExactMatchingChars (minCharsPerMismatch): " << (int) readsExactMatchingChars <<
             "(" << (int) minCharsPerMismatch << ")" << endl;
        cout << "targetMismatches (maxMismatches): " << (int) targetMismatches << "(" << (int) maxMismatches << ")" << endl;
        matcher->matchConstantLengthReads();
        if (dumpInfo)
            matcher->writeMatchesInfo(pgDestFilePrefix);

        const vector<bool> res = matcher->getMatchedReadsBitmap();

        if (matchPrefixLength == DefaultReadsMatcher::DISABLED_PREFIX_MODE)
            matcher->writeIntoPseudoGenome(pgDestFilePrefix, orgIndexesMapping);

        delete(matcher);
        return res;
    }

}