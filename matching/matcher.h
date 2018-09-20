#ifndef PGMATCHER_MATCHER_H
#define PGMATCHER_MATCHER_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void exactMatchConstantLengthPatterns(string text, string readsFile, ofstream& offsetsDest,
                                      ofstream& missedPatternsDest);

void approxMatchConstantLengthPatterns(string text, string readsFile, ofstream& offsetsDest, uint8_t max_mismatches,
                                      ofstream& missedPatternsDest);

#endif //PGMATCHER_MATCHER_H
