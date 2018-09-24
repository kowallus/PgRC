#ifndef PGMATCHER_MATCHER_H
#define PGMATCHER_MATCHER_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void exactMatchConstantLengthPatterns(string text, string readsFile, ofstream& offsetsDest, uint32_t matchPrefixLength,
        ofstream& missedPatternsDest);

uint8_t countMismatches(const char* pattern, const char* text, uint64_t length, uint8_t max_mismatches = UINT8_MAX);

void approxMatchConstantLengthPatterns(string text, string readsFile, ofstream& offsetsDest, uint8_t max_mismatches,
        uint32_t matchPrefixLength, ofstream& missedPatternsDest);

#endif //PGMATCHER_MATCHER_H
