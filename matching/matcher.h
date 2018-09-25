#ifndef PGMATCHER_MATCHER_H
#define PGMATCHER_MATCHER_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

namespace PgTools {

    void
    exactMatchConstantLengthReads(string text, string readsFile, ofstream &offsetsDest, uint32_t matchPrefixLength,
                                  ofstream &missedPatternsDest, ofstream &suffixesDest);

    uint8_t countMismatches(const char *pattern, const char *text, uint64_t length, uint8_t maxMismatches = UINT8_MAX);

    void approxMatchConstantLengthReads(string text, string readsFile, ofstream &offsetsDest, uint8_t maxMismatches,
                                        uint32_t matchPrefixLength, ofstream &missedPatternsDest,
                                        ofstream &suffixesDest);

}

#endif //PGMATCHER_MATCHER_H
