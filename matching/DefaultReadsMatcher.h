#ifndef PGTOOLS_DEFAULTREADSMATCHER_H
#define PGTOOLS_DEFAULTREADSMATCHER_H

#include <iostream>
#include <fstream>
#include <string>


#include "../readsset/PackedReadsSet.h"

using namespace std;

namespace PgTools {

    uint8_t
    countMismatches(const char *pattern, const char *text, uint64_t length, uint8_t maxMismatches = UINT8_MAX);

    class DefaultReadsMatcher {
    private:
        string pgFilePrefix;
        bool revComplPg;
        PackedReadsSet *readsSet;
        uint32_t matchPrefixLength;
        uint8_t maxMismatches;
        uint8_t minMismatches = 0;

        string text;
        vector<uint32_t> readMatchPos;

        vector<uint8_t> readMismatches;
        uint_reads_cnt_max matchedReadsCount = 0;

        uint_reads_cnt_max readsCount;
        uint_read_len_max readLength;
        uint_read_len_max matchingLength;

        void exactMatchConstantLengthReads();

        void approxMatchConstantLengthReads();

        void writeExactMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest);
        void writeApproxMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest);

    public:
        static const uint_read_len_max DISABLED_PREFIX_MODE;
        static const uint32_t NOT_MATCHED_VALUE;

        DefaultReadsMatcher(const string &pgFilePrefix, bool revComplPg, PackedReadsSet *readsSet,
                            uint32_t matchPrefixLength, uint8_t maxMismatches);

        void matchConstantLengthReads();

        static const string OFFSETS_SUFFIX;
        static const string SUFFIXES_SUFFIX;
        static const string MISSED_READS_SUFFIX;

        void writeMatchesInfo(const string &outPrefix);

        const vector<uint_reads_cnt_max> getMatchedReadsIndexes() const;
        const vector<uint32_t> &getReadMatchPos() const;

    };

}

#endif //PGTOOLS_DEFAULTREADSMATCHER_H
