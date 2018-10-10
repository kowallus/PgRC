#include "DefaultReadsMatcher.h"

#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include "../readsset/PackedReadsSet.h"
#include "../readsset/persistance/ReadsSetPersistence.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

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

    const uint_read_len_max DefaultReadsMatcher::DISABLED_PREFIX_MODE = (uint_read_len_max) -1;
    const uint32_t DefaultReadsMatcher::NOT_MATCHED_VALUE = UINT32_MAX;

    DefaultReadsMatcher::DefaultReadsMatcher(const string &pgFilePrefix, bool revComplPg, PackedReadsSet *readsSet,
                                             uint32_t matchPrefixLength, uint8_t maxMismatches) :
                                             pgFilePrefix(pgFilePrefix), revComplPg(revComplPg), readsSet(readsSet),
                                             matchPrefixLength(matchPrefixLength), maxMismatches(maxMismatches) {
        text = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(pgFilePrefix);
        if (revComplPg)
            text = text + "XXXXXX" + PgSAHelpers::reverseComplement(text);

        readLength = readsSet->readLength(0);
        matchingLength = readLength > matchPrefixLength ? matchPrefixLength : readLength;
        readsCount = readsSet->readsCount();
    }

    void DefaultReadsMatcher::matchConstantLengthReads() {
        if (maxMismatches)
            approxMatchConstantLengthReads();
        else
            exactMatchConstantLengthReads();

    }

    void DefaultReadsMatcher::exactMatchConstantLengthReads() {
        clock_checkpoint();

        cout << "Feeding patterns...\n" << endl;
        ConstantLengthPatternsOnTextHashMatcher hashMatcher(matchingLength);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            const string &read = readsSet->getRead(i);
            hashMatcher.addPattern(read.data(), i);
        }
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        cout << "Matching...\n" << endl;
        hashMatcher.iterateOver(text.data(), text.length());

        readMatchPos.clear();
        readMatchPos.insert(readMatchPos.begin(), readsCount, NOT_MATCHED_VALUE);
        matchedReadsCount = 0;

        int i = 0;
        uint64_t multiMatchCount = 0;
        uint64_t falseMatchCount = 0;
        while (hashMatcher.moveNext()) {
            const uint64_t matchPosition = hashMatcher.getHashMatchTextPosition();
            const uint32_t matchReadIndex = hashMatcher.getHashMatchPatternIndex();
            const string &matchedRead = readsSet->getRead(matchReadIndex);

            bool exactMatch = strncmp(text.data() + matchPosition, matchedRead.data(), matchingLength) == 0;
            if (exactMatch) {
                if (readMatchPos[matchReadIndex] == NOT_MATCHED_VALUE) {
                    readMatchPos[matchReadIndex] = matchPosition;
                    matchedReadsCount++;
                } else
                    multiMatchCount++;
            } else
                falseMatchCount++;
            if (i++ < 1) {
                cout << "Matched: " << matchReadIndex << "; "
                     << matchPosition << "; " << exactMatch << endl;
                cout << matchedRead << endl;
                const basic_string<char, char_traits<char>, allocator<char>> &pgPart = text.substr(matchPosition,
                                                                                                   matchingLength);
                cout << pgPart << endl;
            }
        }
        cout << "... exact matching procedure completed in " << clock_millis() << " msec. " << endl;
        cout << "Exact matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;
    }

    void DefaultReadsMatcher::approxMatchConstantLengthReads() {

        clock_checkpoint();

        cout << "Feeding patterns...\n" << endl;
        const uint_read_len_max partLength = matchingLength / (maxMismatches + 1);
        ConstantLengthPatternsOnTextHashMatcher hashMatcher(partLength);

        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            const string &read = readsSet->getRead(i);
            const char *readPtr = read.data();
            for (uint8_t j = 0; j <= maxMismatches; j++, readPtr += partLength)
                hashMatcher.addPattern(readPtr, i * (maxMismatches + 1) + j);
        }
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        cout << "Matching...\n" << endl;
        hashMatcher.iterateOver(text.data(), text.length());

        readMatchPos.clear();
        readMatchPos.insert(readMatchPos.end(), readsCount, NOT_MATCHED_VALUE);
        readMismatches.clear();
        readMismatches.insert(readMismatches.end(), readsCount, UINT8_MAX);
        matchedReadsCount = 0;

        int i = 0;
        uint64_t falseMatchCount = 0;
        uint64_t multiMatchCount = 0;
        while (hashMatcher.moveNext()) {
            const uint32_t matchPatternIndex = hashMatcher.getHashMatchPatternIndex();
            uint32_t matchReadIndex = matchPatternIndex / (maxMismatches + 1);
            if (readMismatches[matchReadIndex] <= minMismatches)
                continue;
            const string &matchedRead = readsSet->getRead(matchReadIndex);
            uint64_t matchPosition = hashMatcher.getHashMatchTextPosition();
            const uint8_t positionShift = ((matchPatternIndex % (maxMismatches + 1)) * partLength);
            if (positionShift > matchPosition)
                continue;
            matchPosition -= positionShift;
            if (matchPosition + readLength > text.length())
                continue;
            if (readMatchPos[matchReadIndex] == matchPosition)
                continue;
            const uint8_t mismatchesCount = countMismatches(matchedRead.data(), text.data() + matchPosition,
                                                            matchingLength, maxMismatches);
            if (mismatchesCount < readMismatches[matchReadIndex]) {
                if (readMismatches[matchReadIndex] == UINT8_MAX)
                    matchedReadsCount++;
                else
                    multiMatchCount++;
                readMatchPos[matchReadIndex] = matchPosition;
                readMismatches[matchReadIndex] = mismatchesCount;
            } else if (mismatchesCount == UINT8_MAX)
                falseMatchCount++;
            else
                multiMatchCount++;
            if (mismatchesCount < UINT8_MAX && i++ < 2) {
                cout << "Matched: " << matchReadIndex << " (" << matchPatternIndex << "); "
                     << matchPosition << "; " << (int) mismatchesCount << endl;
                cout << matchedRead << endl;
                const basic_string<char, char_traits<char>, allocator<char>> &pgPart = text.substr(matchPosition,
                                                                                                   matchingLength);
                cout << pgPart << endl;
            }
        }
        cout << "... approximate matching procedure completed in " << clock_millis() << " msec. " << endl;
        cout << "Matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;
    }

    const string DefaultReadsMatcher::OFFSETS_SUFFIX = "_matched_offsets.txt";
    const string DefaultReadsMatcher::SUFFIXES_SUFFIX = "_matched_suffixes.txt";
    const string DefaultReadsMatcher::MISSED_READS_SUFFIX = "_missed.txt";

    void DefaultReadsMatcher::writeExactMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest) {
        clock_checkpoint();

        cout << "Writing output files...\n" << endl;
        for (uint_reads_cnt_max i = 0; i < readMatchPos.size(); i++) {
            if (readMatchPos[i] == NOT_MATCHED_VALUE)
                missedPatternsDest << readsSet->getRead(i) << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << "\n";
                if (matchingLength < readLength)
                    suffixesDest << readsSet->getRead(i).substr(matchingLength);
            }
        }

        cout << "... writing output files completed in  " << clock_millis() << " msec. " << endl;
    }

    void DefaultReadsMatcher::writeApproxMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest) {
        clock_checkpoint();

        cout << "Writing output files...\n" << endl;
        vector<uint_reads_cnt_max> mismatchedReadsCount(maxMismatches + 1, 0);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            const string &read = readsSet->getRead(i);
            if (readMatchPos[i] == NOT_MATCHED_VALUE)
                missedPatternsDest << read << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i];
                reportMismatches(read.data(), text.data() + readMatchPos[i], matchingLength, offsetsDest);
                if (matchingLength < readLength)
                    suffixesDest << readsSet->getRead(i).substr(matchingLength);
                offsetsDest << "\n";
                mismatchedReadsCount[readMismatches[i]]++;
            }
        }

        for (uint8_t i = 0; i <= maxMismatches; i++)
            cout << "Matched " << mismatchedReadsCount[i] << " reads with " << (int) i << " mismatches." << endl;

        cout << "... writing output files completed in  " << clock_millis() << " msec. " << endl;
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

        if (maxMismatches == 0)
            writeExactMatchesInfo(offsetsDest, missedReadsDest, suffixesDest);
        else
            writeApproxMatchesInfo(offsetsDest, missedReadsDest, suffixesDest);

        offsetsDest.close();
        missedReadsDest.close();

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

}