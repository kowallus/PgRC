#include "ReadsMatchers.h"

#include "../readsset/PackedReadsSet.h"
#include "../readsset/persistance/ReadsSetPersistence.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../pseudogenome/persistence/SeparatedExtendedReadsList.h"

namespace PgTools {

    uint8_t countMismatches(const char *pattern, const char *text, uint64_t length, uint8_t maxMismatches) {
        uint8_t res = 0;
        const char *patEnd = pattern + length;
        while (pattern != patEnd) {
            if (*pattern++ != *text++) {
                if (res++ >= maxMismatches)
                    return UINT8_MAX;
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
    const uint32_t DefaultReadsMatcher::NOT_MATCHED_VALUE = UINT32_MAX;

    DefaultReadsMatcher::DefaultReadsMatcher(const string &pgFilePrefix, bool revComplPg, PackedReadsSet *readsSet,
                                             uint32_t matchPrefixLength) :
                                             pgFilePrefix(pgFilePrefix), revComplPg(revComplPg), readsSet(readsSet),
                                             matchPrefixLength(matchPrefixLength) {

        readLength = readsSet->readLength(0);
        matchingLength = readLength > matchPrefixLength ? matchPrefixLength : readLength;
        readsCount = readsSet->readsCount();
    }

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
        readMatchPos.insert(readMatchPos.begin(), readsCount, NOT_MATCHED_VALUE);
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

        string text = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(pgFilePrefix);
        this->matchConstantLengthReads(text.data(), text.length(), false);

        if (revComplPg) {
            text = PgSAHelpers::reverseComplement(text);
            this->matchConstantLengthReads(text.data(), text.length(), true);
        }

    }

    DefaultReadsExactMatcher::DefaultReadsExactMatcher(const string &pgFilePrefix, bool revComplPg,
                                                       PackedReadsSet *readsSet, uint32_t matchPrefixLength)
            : DefaultReadsMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength) {}

    void DefaultReadsExactMatcher::initMatching() {
        cout << "Feeding patterns...\n" << endl;
        this->hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(matchingLength);
        this->hashMatcher->addPackedPatterns(readsSet);
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        DefaultReadsMatcher::initMatching();
    }

    void DefaultReadsExactMatcher::matchConstantLengthReads(const char* txt, uint64_t length, bool revCompMode) {
        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;

        hashMatcher->iterateOver(txt, length);

        int i = 0;

        while (hashMatcher->moveNext()) {
            const uint64_t matchPosition = hashMatcher->getHashMatchTextPosition();
            const uint32_t matchReadIndex = hashMatcher->getHashMatchPatternIndex();

            bool exactMatch = readsSet->comparePackedReadWithPattern(matchReadIndex, txt + matchPosition) == 0;
            if (exactMatch) {
                if (readMatchPos[matchReadIndex] == NOT_MATCHED_VALUE) {
                    readMatchPos[matchReadIndex] = revCompMode?length-(matchPosition+matchingLength):matchPosition;
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

    AbstractReadsApproxMatcher::AbstractReadsApproxMatcher(const string &pgFilePrefix, bool revComplPg,
                                                           PackedReadsSet *readsSet, uint32_t matchPrefixLength,
                                                           uint8_t targetMismatches, uint8_t maxMismatches, uint8_t minMismatches)
            : DefaultReadsMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength), targetMismatches(targetMismatches),
            maxMismatches(maxMismatches), minMismatches(minMismatches){
        currentRead.resize(readsSet->getReadsSetProperties()->maxReadLength);
    }

    InterleavedReadsApproxMatcher::InterleavedReadsApproxMatcher(const string &pgFilePrefix, bool revComplPg,
                                                                 PackedReadsSet *readsSet, uint32_t matchPrefixLength,
                                                uint8_t targetMismatches, uint8_t maxMismatches, uint8_t minMismatches)
            : AbstractReadsApproxMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength, targetMismatches,
                    maxMismatches, minMismatches)  {}

    DefaultReadsApproxMatcher::DefaultReadsApproxMatcher(const string &pgFilePrefix, bool revComplPg,
                                                         PackedReadsSet *readsSet, uint32_t matchPrefixLength,
                                                         uint8_t targetMismatches, uint8_t maxMismatches, uint8_t minMismatches)
             :AbstractReadsApproxMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength, targetMismatches,
                     maxMismatches, minMismatches) {}




    void DefaultReadsApproxMatcher::initMatching() {
        cout << "Feeding patterns...\n" << endl;
        partLength = matchingLength / (targetMismatches + 1);
        hashMatcher = new DefaultConstantLengthPatternsOnTextHashMatcher(partLength);
        this->hashMatcher->addPackedPatterns(readsSet, targetMismatches + 1);
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;

        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, UINT8_MAX);
    }

    void DefaultReadsApproxMatcher::matchConstantLengthReads(const char* txt, uint64_t length, bool revCompMode) {

        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;
        hashMatcher->iterateOver(txt, length);

        int i = 0;
        while (hashMatcher->moveNext()) {
            const uint32_t matchPatternIndex = hashMatcher->getHashMatchPatternIndex();
            uint32_t matchReadIndex = matchPatternIndex / (targetMismatches + 1);
            if (readMismatchesCount[matchReadIndex] <= minMismatches)
                continue;
            uint64_t matchPosition = hashMatcher->getHashMatchTextPosition();
            const uint8_t positionShift = ((matchPatternIndex % (targetMismatches + 1)) * partLength);
            if (positionShift > matchPosition)
                continue;
            matchPosition -= positionShift;
            if (matchPosition + readLength > length)
                continue;
            if (readMatchPos[matchReadIndex] == (revCompMode?length-(matchPosition+matchingLength):matchPosition))
                continue;
            const uint8_t mismatchesCount = readsSet->countMismatchesVsPattern(matchReadIndex, txt + matchPosition,
                                                            matchingLength, maxMismatches);
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == UINT8_MAX)
                    matchedReadsCount++;
                else
                    multiMatchCount++;
                readMatchPos[matchReadIndex] = revCompMode?length-(matchPosition+matchingLength):matchPosition;
                if (revCompMode) readMatchRC[matchReadIndex] = true;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            } else if (mismatchesCount == UINT8_MAX)
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
        cout << "... approximate matching procedure completed in " << clock_millis() << " msec. " << endl;
        cout << "Matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;
    }

    void InterleavedReadsApproxMatcher::initMatching() {
        cout << "Feeding patterns...\n" << endl;
        partLength = matchingLength / (targetMismatches + 1);
        hashMatcher = new InterleavedConstantLengthPatternsOnTextHashMatcher(partLength, targetMismatches + 1);
        this->hashMatcher->addPackedPatterns(readsSet, targetMismatches + 1);
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;

        DefaultReadsMatcher::initMatching();
        readMismatchesCount.clear();
        readMismatchesCount.insert(readMismatchesCount.end(), readsCount, UINT8_MAX);
    }

    void InterleavedReadsApproxMatcher::matchConstantLengthReads(const char *txt, uint64_t length, bool revCompMode) {

        cout << "Matching" << (revCompMode?" in Pg reverse":"") << "...\n" << endl;
        hashMatcher->iterateOver(txt, length);

        int i = 0;
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
            if (matchPosition + readLength > length)
                continue;
            if (readMatchPos[matchReadIndex] == (revCompMode?length-(matchPosition+matchingLength):matchPosition))
                continue;
            const uint8_t mismatchesCount = readsSet->countMismatchesVsPattern(matchReadIndex, txt + matchPosition,
                                                            matchingLength, maxMismatches);
            if (mismatchesCount < readMismatchesCount[matchReadIndex]) {
                if (readMismatchesCount[matchReadIndex] == UINT8_MAX)
                    matchedReadsCount++;
                else
                    multiMatchCount++;
                readMatchPos[matchReadIndex] = revCompMode?length-(matchPosition+matchingLength):matchPosition;
                if (revCompMode) readMatchRC[matchReadIndex] = true;
                readMismatchesCount[matchReadIndex] = mismatchesCount;
            } else if (mismatchesCount == UINT8_MAX)
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
        cout << "... approximate matching procedure completed in " << clock_millis() << " msec. " << endl;
        cout << "Matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;
    }

    const string DefaultReadsMatcher::OFFSETS_SUFFIX = "_matched_offsets.txt";
    const string DefaultReadsMatcher::SUFFIXES_SUFFIX = "_matched_suffixes.txt";
    const string DefaultReadsMatcher::MISSED_READS_SUFFIX = "_missed.txt";

    void DefaultReadsExactMatcher::writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest) {
        clock_checkpoint();

        for (uint_reads_cnt_max i = 0; i < readMatchPos.size(); i++) {
            if (readMatchPos[i] == NOT_MATCHED_VALUE)
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

        string text = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(pgFilePrefix);
        vector<uint_reads_cnt_max> mismatchedReadsCount(maxMismatches + 1, 0);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            if (readMatchPos[i] == NOT_MATCHED_VALUE)
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

        for (uint8_t i = 0; i <= maxMismatches; i++)
            cout << "Matched " << mismatchedReadsCount[i] << " reads with " << (int) i << " mismatches." << endl;

        cout << "... writing info dump files completed in " << clock_millis() << " msec. " << endl;
    }

    const vector<uint_reads_cnt_max> DefaultReadsMatcher::getMatchedReadsIndexes() const {
        vector<uint_reads_cnt_max> matchedReads(matchedReadsCount);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            if (readMatchPos[i] != NOT_MATCHED_VALUE)
                matchedReads.push_back(i);
        }
        return matchedReads;
    }

    const vector<uint32_t> &DefaultReadsMatcher::getReadMatchPos() const {
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

    void AbstractReadsApproxMatcher::initEntryUpdating() {
        pg = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(pgFilePrefix);
    }

    void AbstractReadsApproxMatcher::updateEntry(DefaultReadsListEntry &entry, uint_reads_cnt_max matchIdx) {
        readsSet->getRead(matchIdx, (char_pg*) currentRead.data());
        if (readMatchRC[matchIdx])
            reverseComplementInPlace(currentRead);
        fillEntryWithMismatches(currentRead.data(), pg.data() + entry.pos, readMismatchesCount[matchIdx], entry);
    }

    void AbstractReadsApproxMatcher::closeEntryUpdating() {
        pg.clear();
    }

    void DefaultReadsMatcher::writeIntoPseudoGenome(const string &outPgPrefix, IndexesMapping* orgIndexesMapping) {
        clock_checkpoint();
        vector<uint_reads_cnt_max> idxs(matchedReadsCount);
        uint64_t counter = 0;
        for(uint_reads_cnt_max i = 0; i < readsCount; i++)
            if (readMatchPos[i] != NOT_MATCHED_VALUE)
                idxs[counter++] = i;

        std::sort(idxs.begin(), idxs.end(), [this](const uint_reads_cnt_max& idx1, const uint_reads_cnt_max& idx2) -> bool
        { return readMatchPos[idx1] < readMatchPos[idx2]; });

        initEntryUpdating();
        DefaultSeparatedExtendedReadsListIterator* rlIt = DefaultSeparatedExtendedReadsListIterator::getIterator(pgFilePrefix);
        SeparatedPseudoGenomeOutputBuilder* builder = this->createSeparatedPseudoGenomeOutputBuilder(outPgPrefix,
                rlIt->isRevCompEnabled(), rlIt->areMismatchesEnabled());
        builder->setReadsSourceIterator(rlIt);
        builder->copyPseudoGenomeHeader(pgFilePrefix);
        for(uint_reads_cnt_max i = 0; i < matchedReadsCount; i++) {
            uint_reads_cnt_max matchIdx = idxs[i];
            uint64_t currPos = builder->writeReadsFromIterator(readMatchPos[matchIdx]);
            DefaultReadsListEntry entry(currPos);
            entry.advanceEntryByPosition(readMatchPos[matchIdx], orgIndexesMapping->getReadOriginalIndex(matchIdx), readMatchRC[matchIdx]);
            this->updateEntry(entry, matchIdx);
            builder->writeExtraReadEntry(entry);
        }
        builder->writeReadsFromIterator();
        delete(rlIt);
        builder->build();
        delete(builder);
        closeEntryUpdating();
        cout << "... writing (" << outPgPrefix << ") output files completed in " << clock_millis() << " msec. " << endl << endl;
    }

    void mapReadsIntoPg(const string &pgFilePrefix, bool revComplPg, PackedReadsSet *readsSet,
                        uint_read_len_max matchPrefixLength, uint8_t targetMismatches, uint8_t maxMismatches,
                        char mismatchesMode, uint8_t minMismatches, bool dumpInfo, const string &pgDestFilePrefix,
                        IndexesMapping* orgIndexesMapping, bool divisionComplement,
                        const string &outDivisionFile) {
        DefaultReadsMatcher* matcher;
        cout << "targetMismatches (maxMismatches): " << (int) targetMismatches << "(" << (int) maxMismatches << ")" << endl;
        if (targetMismatches == 0)
            matcher = new DefaultReadsExactMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength);
        else switch (mismatchesMode) {
                case 'd': matcher = new DefaultReadsApproxMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength,
                                                                      targetMismatches, maxMismatches, minMismatches);
                    break;
                case 'i': matcher = new InterleavedReadsApproxMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength,
                                                            targetMismatches, maxMismatches, minMismatches);
                    break;
                default:
                    fprintf(stderr, "Unknown mismatches mode: %c.\n", mismatchesMode);
                    exit(EXIT_FAILURE);

            }
        matcher->matchConstantLengthReads();
        if (dumpInfo)
            matcher->writeMatchesInfo(pgDestFilePrefix);

        const vector<uint32_t> &readsMatchPos = matcher->getReadMatchPos();
        const vector<uint8_t> &readsMismatches = matcher->getReadMismatches();

        ReadsSetPersistence::writeOutputDivision(orgIndexesMapping, readsMatchPos,
                                                 DefaultReadsMatcher::NOT_MATCHED_VALUE, outDivisionFile, divisionComplement);

        if (matchPrefixLength == DefaultReadsMatcher::DISABLED_PREFIX_MODE)
            matcher->writeIntoPseudoGenome(pgDestFilePrefix, orgIndexesMapping);

        delete(matcher);
    }

}