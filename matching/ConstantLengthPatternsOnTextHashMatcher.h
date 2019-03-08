#ifndef PGTOOLS_CONSTANTLENGTHPATTERNSONTEXTHASHMATCHER_H
#define PGTOOLS_CONSTANTLENGTHPATTERNSONTEXTHASHMATCHER_H

#include <unordered_map>
#include <vector>
#include "../rollinghash/cyclichash.h"
#include "../readsset/ReadsSetInterface.h"

using namespace std;

class DefaultConstantLengthPatternsOnTextHashMatcher {
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
    DefaultConstantLengthPatternsOnTextHashMatcher(uint32_t patternLength);

    virtual ~DefaultConstantLengthPatternsOnTextHashMatcher();

    void addPattern(const char* pattern, uint32_t idx);
    void addReadsSetOfPatterns(ConstantLengthReadsSetInterface *readsSet, uint8_t partsCount = 1,
                               vector<bool> matchedReadsBitmap = {});

    //iterator routines
    inline void iterateOver(const char* txt, uint64_t length);
    inline bool moveNext();
    uint32_t getHashMatchPatternIndex();
    uint64_t getHashMatchTextPosition();
};


void DefaultConstantLengthPatternsOnTextHashMatcher::iterateOver(const char *txt, uint64_t length) {
    this->txt = txt;
    this->txtSize = length;
    hf.reset();
    for(uint32_t i = 0; i < txtSize && i < patternLength; i++)
        hf.eat(this->txt[i]);
    this->txtPos = -1;
    indexIter = hashToIndexMap.end();
    indexIterEnd = hashToIndexMap.end();
}

bool DefaultConstantLengthPatternsOnTextHashMatcher::moveNext() {
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

class InterleavedConstantLengthPatternsOnTextHashMatcher {
private:
    unordered_multimap<uint32_t, uint32_t> hashToIndexMap;
    const uint32_t patternLength;
    const uint8_t patternParts;
    const uint32_t patternSpan;

    std::vector<CyclicHash<uint32_t>> hf;
    uint8_t currentHF = 0;

    const  char* txt = 0;
    uint64_t txtSize = 0;

    //iterator fields
    int64_t txtPos = -1;
    unordered_multimap<uint32_t, uint32_t>::iterator indexIter, indexIterEnd;

public:
    InterleavedConstantLengthPatternsOnTextHashMatcher(uint32_t patternLength, const uint8_t patternParts);

    virtual ~InterleavedConstantLengthPatternsOnTextHashMatcher();

    void addPattern(const char* pattern, uint32_t idx);

    //iterator routines
    inline void iterateOver(const char* txt, uint64_t length);
    inline bool moveNext();
    uint32_t getHashMatchPatternIndex();
    uint64_t getHashMatchTextPosition();

    void addPackedPatterns(ConstantLengthReadsSetInterface*readsSet, int partsCount,
                           vector<bool> matchedReadsBitmap = {});
};


void InterleavedConstantLengthPatternsOnTextHashMatcher::iterateOver(const char *txt, uint64_t length) {
    this->txt = txt;
    this->txtSize = length;
    for(uint8_t h = 0; h < patternParts; h++) {
        hf[h].reset();
        const uint32_t patternGuard = patternSpan + h;
        for (uint32_t i = h; i < txtSize && i < patternGuard; i += patternParts)
            hf[h].eat(this->txt[i]);
    }
    this->txtPos = -1;
    indexIter = hashToIndexMap.end();
    indexIterEnd = hashToIndexMap.end();
    currentHF = 0;
}

bool InterleavedConstantLengthPatternsOnTextHashMatcher::moveNext() {
    if (indexIter != indexIterEnd) {
        indexIter++;
        if (indexIter != indexIterEnd)
            return true;
    }
    while(++this->txtPos <= txtSize - patternSpan) {
        auto indexIterRange = hashToIndexMap.equal_range(hf[currentHF].hashvalue);
        hf[currentHF++].update(this->txt[this->txtPos], this->txt[this->txtPos + patternSpan]);
        if (currentHF == patternParts)
            currentHF = 0;
        indexIter = indexIterRange.first;
        indexIterEnd = indexIterRange.second;
        if (indexIter != indexIterEnd)
            return true;
    }
    return false;
}

#endif //PGTOOLS_CONSTANTLENGTHPATTERNSONTEXTHASHMATCHER_H
