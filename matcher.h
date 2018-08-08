#ifndef PGMATCHER_MATCHER_H
#define PGMATCHER_MATCHER_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void exactMatchConstantLengthPatterns(string text, ifstream& patternsSrc, ofstream& offsetsDest,
                                      ofstream& missedPatternsDest);

#endif //PGMATCHER_MATCHER_H
