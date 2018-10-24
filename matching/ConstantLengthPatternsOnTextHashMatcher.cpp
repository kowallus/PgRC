#include "ConstantLengthPatternsOnTextHashMatcher.h"

DefaultConstantLengthPatternsOnTextHashMatcher::DefaultConstantLengthPatternsOnTextHashMatcher(uint32_t patternLength)
        : patternLength(patternLength), hf(patternLength, 32) {
}

DefaultConstantLengthPatternsOnTextHashMatcher::~DefaultConstantLengthPatternsOnTextHashMatcher() {

}

void DefaultConstantLengthPatternsOnTextHashMatcher::addPattern(const char *pattern, uint32_t idx) {
    if (this->txt != 0) {
        cerr << "Adding patterns not permitted during iteration";
        exit(EXIT_FAILURE);
    }
    hf.reset();
    for(uint32_t i = 0; i < patternLength; i++)
        hf.eat(pattern[i]);
    hashToIndexMap.insert(std::pair<uint32_t, uint32_t>(hf.hashvalue, idx));
}

uint32_t DefaultConstantLengthPatternsOnTextHashMatcher::getHashMatchPatternIndex() {
    return indexIter->second;
}

uint64_t DefaultConstantLengthPatternsOnTextHashMatcher::getHashMatchTextPosition() {
    return this->txtPos;
}

InterleavedConstantLengthPatternsOnTextHashMatcher::InterleavedConstantLengthPatternsOnTextHashMatcher(
        uint32_t patternLength, const uint8_t patternParts)
        : patternLength(patternLength), patternParts(patternParts), patternSpan(patternLength*patternParts) {
    hf.push_back(CyclicHash<uint32_t>(patternLength, 32));
    for(uint8_t i = 1; i < patternParts; i++)
        hf.push_back(hf[0]);
}

InterleavedConstantLengthPatternsOnTextHashMatcher::~InterleavedConstantLengthPatternsOnTextHashMatcher() {

}

void InterleavedConstantLengthPatternsOnTextHashMatcher::addPattern(const char *pattern, uint32_t idx) {
    if (this->txt != 0) {
        cerr << "Adding patterns not permitted during iteration";
        exit(EXIT_FAILURE);
    }
    hf[0].reset();
    for(uint32_t i = 0; i < patternSpan; i += patternParts)
        hf[0].eat(pattern[i]);
    hashToIndexMap.insert(std::pair<uint32_t, uint32_t>(hf[0].hashvalue, idx));
}

uint32_t InterleavedConstantLengthPatternsOnTextHashMatcher::getHashMatchPatternIndex() {
    return indexIter->second;
}

uint64_t InterleavedConstantLengthPatternsOnTextHashMatcher::getHashMatchTextPosition() {
    return this->txtPos;
}
