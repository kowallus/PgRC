#ifndef PGTOOLS_CONSTANTLENGTHPATTERNSONTEXTHASHMATCHER_H
#define PGTOOLS_CONSTANTLENGTHPATTERNSONTEXTHASHMATCHER_H

#include <unordered_map>
#include "../rollinghash/cyclichash.h"

using namespace std;

class ConstantLengthPatternsOnTextHashMatcher {
private:
    unordered_multimap<uint32_t, uint32_t> hashToIndexMap;
    const uint32_t patternLength;

    CyclicHash<uint32_t> hf;

    const  char* txt = 0;
    uint64_t txtSize = 0;

    //iterator fields
    int64_t txtPos = -1;
    unordered_multimap<uint32_t, uint32_t>::iterator indexIter, indexIterEnd;

public:
    ConstantLengthPatternsOnTextHashMatcher(uint32_t patternLength);

    virtual ~ConstantLengthPatternsOnTextHashMatcher();

    void addPattern(const char* pattern, uint32_t idx);

    //iterator routines
    inline void iterateOver(const char* txt, uint64_t length);
    inline bool moveNext();
    uint32_t getHashMatchPatternIndex();
    uint64_t getHashMatchTextPosition();
};


void ConstantLengthPatternsOnTextHashMatcher::iterateOver(const char *txt, uint64_t length) {
    this->txt = txt;
    this->txtSize = length;
    hf.reset();
    for(uint32_t i = 0; i < txtSize && i < patternLength; i++)
        hf.eat(this->txt[i]);
    this->txtPos = -1;
    indexIter = hashToIndexMap.end();
    indexIterEnd = hashToIndexMap.end();
}

bool ConstantLengthPatternsOnTextHashMatcher::moveNext() {
    if (indexIter != indexIterEnd) {
        indexIter++;
        if (indexIter != indexIterEnd)
            return true;
    }
    while(++this->txtPos <= txtSize - patternLength) {
        auto indexIterRange = hashToIndexMap.equal_range(hf.hashvalue);
        hf.update(this->txt[this->txtPos], this->txt[this->txtPos + patternLength]);
        indexIter = indexIterRange.first;
        indexIterEnd = indexIterRange.second;
        if (indexIter != indexIterEnd)
            return true;
    }
    return false;
}

#endif //PGTOOLS_CONSTANTLENGTHPATTERNSONTEXTHASHMATCHER_H
