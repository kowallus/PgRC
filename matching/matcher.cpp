#include "matcher.h"

#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include "../readsset/PackedReadsSet.h"

namespace PgTools {

    void exactMatchConstantLengthReads(string text, string readsFile, ofstream &offsetsDest,
                                       uint32_t matchPrefixLength, ofstream &missedPatternsDest,
                                       ofstream &suffixesDest) {
        clock_checkpoint();
        cout << "Reading reads set\n";
        PackedReadsSet *readsSet = PackedReadsSet::readReadsSet(readsFile);
        readsSet->printout();
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        cout << "Feeding patterns...\n" << endl;
        const uint_read_len_max readLength = readsSet->readLength(0);
        const uint_read_len_max matchingLength = readLength > matchPrefixLength ? matchPrefixLength : readLength;
        ConstantLengthPatternsOnTextHashMatcher hashMatcher(matchingLength);
        const uint_reads_cnt_max readsCount = readsSet->readsCount();
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            const string &read = readsSet->getRead(i);
            hashMatcher.addPattern(read.data(), i);
        }
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        cout << "Matching...\n" << endl;
        hashMatcher.iterateOver(text.data(), text.length());

        vector<uint32_t> readMatchPos(readsCount, UINT32_MAX);

        int i = 0;
        uint32_t matchCount = 0;
        uint64_t multiMatchCount = 0;
        uint64_t falseMatchCount = 0;
        while (hashMatcher.moveNext()) {
            const uint64_t matchPosition = hashMatcher.getHashMatchTextPosition();
            const uint32_t matchReadIndex = hashMatcher.getHashMatchPatternIndex();
            const string &matchedRead = readsSet->getRead(matchReadIndex);

            bool exactMatch = strncmp(text.data() + matchPosition, matchedRead.data(), matchingLength) == 0;
            if (exactMatch) {
                if (readMatchPos[matchReadIndex] == UINT32_MAX) {
                    readMatchPos[matchReadIndex] = matchPosition;
                    matchCount++;
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
        cout << "... finished matching in  " << clock_millis() << " msec. " << endl;
        cout << "Exact matched " << matchCount << " reads (" << (readsCount - matchCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;

        cout << "Writing output files...\n" << endl;
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            if (readMatchPos[i] == UINT32_MAX)
                missedPatternsDest << readsSet->getRead(i) << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << "\n";
                if (matchingLength < readLength)
                    suffixesDest << readsSet->getRead(i).substr(matchingLength);
            }
        }

        cout << "... matching and writing output files completed in  " << clock_millis() << " msec. " << endl;

        delete (readsSet);
    }

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

    void
    reportMismatches(const char *read, const char *pgPart, const uint_read_len_max length, ofstream &mismatchesDest) {
        uint64_t pos = 0;
        do {
            if (read[pos] != pgPart[pos])
                mismatchesDest << "\t" << pos << "\t" << read[pos];
        } while (++pos < length);
    }

    void approxMatchConstantLengthReads(string text, string readsFile, ofstream &offsetsDest, uint8_t maxMismatches,
                                        uint32_t matchPrefixLength, ofstream &missedPatternsDest,
                                        ofstream &suffixesDest) {
        uint8_t min_mismatches = 0;

        clock_checkpoint();
        cout << "Reading reads set\n";
        PackedReadsSet *readsSet = PackedReadsSet::readReadsSet(readsFile);
        readsSet->printout();
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        cout << "Feeding patterns...\n" << endl;
        const uint_read_len_max readLength = readsSet->readLength(0);
        const uint_read_len_max matchingLength = readLength > matchPrefixLength ? matchPrefixLength : readLength;
        const uint_read_len_max partLength = matchingLength / (maxMismatches + 1);
        ConstantLengthPatternsOnTextHashMatcher hashMatcher(partLength);
        const uint_reads_cnt_max readsCount = readsSet->readsCount();
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            const string &read = readsSet->getRead(i);
            const char *readPtr = read.data();
            for (uint8_t j = 0; j <= maxMismatches; j++, readPtr += partLength)
                hashMatcher.addPattern(readPtr, i * (maxMismatches + 1) + j);
        }
        cout << "... checkpoint " << clock_millis() << " msec. " << endl;
        cout << "Matching...\n" << endl;
        hashMatcher.iterateOver(text.data(), text.length());

        vector<uint32_t> readMatchPos(readsCount, UINT32_MAX);
        vector<uint8_t> readMismatches(readsCount, UINT8_MAX);

        int i = 0;
        uint32_t matchedReadsCount = 0;
        uint64_t falseMatchCount = 0;
        uint64_t multiMatchCount = 0;
        while (hashMatcher.moveNext()) {
            const uint32_t matchPatternIndex = hashMatcher.getHashMatchPatternIndex();
            uint32_t matchReadIndex = matchPatternIndex / (maxMismatches + 1);
            if (readMismatches[matchReadIndex] <= min_mismatches)
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
        cout << "... finished matching in  " << clock_millis() << " msec. " << endl;
        cout << "Matched " << matchedReadsCount << " reads (" << (readsCount - matchedReadsCount)
             << " left; " << multiMatchCount << " multi-matches). False matches reported: " << falseMatchCount << "."
             << endl;

        cout << "Writing output files...\n" << endl;
        vector<uint_reads_cnt_max> mismatchedReadsCount(maxMismatches + 1, 0);
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            const string &read = readsSet->getRead(i);
            if (readMatchPos[i] == UINT32_MAX)
                missedPatternsDest << read << "\n";
            else {
                offsetsDest << i << "\t" << readMatchPos[i] << "\n";
                reportMismatches(read.data(), text.data() + readMatchPos[i], matchingLength, offsetsDest);
                if (matchingLength < readLength)
                    suffixesDest << readsSet->getRead(i).substr(matchingLength);
                mismatchedReadsCount[readMismatches[i]]++;
            }
        }

        for (uint8_t i = 0; i <= maxMismatches; i++)
            cout << "Matched " << mismatchedReadsCount[i] << " reads with " << (int) i << " mismatches." << endl;

        cout << "... matching and writing output files completed in  " << clock_millis() << " msec. " << endl;

        delete (readsSet);
    }

}