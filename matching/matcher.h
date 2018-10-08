#ifndef PGMATCHER_MATCHER_H
#define PGMATCHER_MATCHER_H

#include <iostream>
#include <fstream>
#include <string>

#include "../readsset/iterator/ReadsSetIterator.h"

using namespace std;

namespace PgTools {

    void
    exactMatchConstantLengthReads(const string &text, ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator,
            ofstream &offsetsDest, uint32_t matchPrefixLength, ofstream &missedPatternsDest, ofstream &suffixesDest);

    uint8_t countMismatches(const char *pattern, const char *text, uint64_t length, uint8_t maxMismatches = UINT8_MAX);

    void approxMatchConstantLengthReads(const string &text, ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator,
            ofstream &offsetsDest, uint8_t maxMismatches, uint32_t matchPrefixLength, ofstream &missedPatternsDest,
                                        ofstream &suffixesDest);

}

#endif //PGMATCHER_MATCHER_H
