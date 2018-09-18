#include "ConstantLengthPatternsOnTextHashMatcher.h"

ConstantLengthPatternsOnTextHashMatcher::ConstantLengthPatternsOnTextHashMatcher(uint32_t patternLength)
        : patternLength(patternLength), hf(patternLength, 32) {}

ConstantLengthPatternsOnTextHashMatcher::~ConstantLengthPatternsOnTextHashMatcher() {

}

void ConstantLengthPatternsOnTextHashMatcher::addPattern(const char *pattern, uint32_t idx) {
    if (this->txt != 0) {
        cerr << "Adding patterns not permitted during interation";
        exit(EXIT_FAILURE);
    }
    hf.reset();
    for(uint32_t i = 0; i < patternLength; i++)
        hf.eat(pattern[i]);
    hashToIndexMap.insert(std::pair<uint32_t, uint32_t>(hf.hashvalue, idx));
}

uint32_t ConstantLengthPatternsOnTextHashMatcher::getHashMatchPatternIndex() {
    return indexIter->second;
}

uint64_t ConstantLengthPatternsOnTextHashMatcher::getHashMatchTextPosition() {
    return this->txtPos - patternLength;
}
