#include "SymbolsPackingFacility.h"

#include <cassert>

namespace PgIndex {

    static const vector<char> binaryCodes { 0, 1 };
    static const vector<char> trenaryCodes { 0, 1, 2 };
    static const vector<char> quaternaryCodes { 0, 1, 2, 3 };

    static const vector<char> acgtSymbols { 'A', 'C', 'G', 'T' };
    static const vector<char> acgtnSymbols { 'A', 'C', 'G', 'N', 'T' };

    SymbolsPackingFacility SymbolsPackingFacility::ACGTPacker(acgtSymbols, true);
    SymbolsPackingFacility SymbolsPackingFacility::ACGTNPacker(acgtnSymbols, true);
    
    SymbolsPackingFacility::SymbolsPackingFacility(const vector<char> symbolsList, bool isGloballyManaged):
         symbolsCount(symbolsList.size()),
         globallyManaged(isGloballyManaged),
         symbolsPerElement(SymbolsPackingFacility::maxSymbolsPerElement(symbolsCount)) {
        uint_max combinationCount = powuint(symbolsCount, symbolsPerElement);
        maxValue = combinationCount - 1;
        if (maxValue > (int) (uint8_t) - 1)
            cout << "ERROR in symbols packaging: max value for type: " << (int) (uint8_t) - 1 << " while max " << " \n";

        reverse = new char_pg*[combinationCount];
        reverseFlat = new char_pg[combinationCount * symbolsPerElement];

        clear = new uint8_t*[combinationCount];
        clearFlat = new uint8_t[combinationCount * symbolsPerElement];

        std::copy(symbolsList.begin(), symbolsList.end(), std::begin(this->symbolsList));
        memset(symbolOrder, -1, UCHAR_MAX);
        for (uint_symbols_cnt i = 0; i < symbolsCount; i++)
            symbolOrder[(unsigned char) symbolsList[(unsigned char) i]] = i;

        buildReversePackAndClearIndexes();
    }

    SymbolsPackingFacility::SymbolsPackingFacility(ReadsSetProperties* readsSetProperties, uchar symbolsPerElement)
    : symbolsCount(readsSetProperties->symbolsCount),
    symbolsPerElement(symbolsPerElement) {
        uint_max combinationCount = powuint(symbolsCount, symbolsPerElement);
        maxValue = combinationCount - 1;
        if (maxValue > (int) (uint8_t) - 1)
            cout << "ERROR in symbols packaging: max value for type: " << (int) (uint8_t) - 1 << " while max " << " \n";

        reverse = new char_pg*[combinationCount];
        reverseFlat = new char_pg[combinationCount * symbolsPerElement];

        clear = new uint8_t*[combinationCount];
        clearFlat = new uint8_t[combinationCount * symbolsPerElement];

        std::copy(std::begin(readsSetProperties->symbolsList), std::end(readsSetProperties->symbolsList), std::begin(symbolsList));
        std::copy(std::begin(readsSetProperties->symbolOrder), std::end(readsSetProperties->symbolOrder), std::begin(symbolOrder));

        buildReversePackAndClearIndexes();
    }

    SymbolsPackingFacility::~SymbolsPackingFacility() {
        delete[]reverse;
        delete[]reverseFlat;
        delete[]clear;
        delete[]clearFlat;
        delete[]packLUT0;
        delete[]packLUT1;
    }

    void SymbolsPackingFacility::buildReversePackAndClearIndexes() {
        char_pg* rPtr = reverseFlat;
        uint8_t* cPtr = clearFlat;
        for (uint_max i = 0; i <= maxValue; i++) {
            reverse[i] = rPtr;
            rPtr += symbolsPerElement;
            clear[i] = cPtr;
            cPtr += symbolsPerElement;
        }

        uint8_t* currentClear = new uint8_t[symbolsPerElement]();
        uint_symbols_cnt* sequence = new uint_symbols_cnt[symbolsPerElement]();
        for (uint_max i = 0; i <= maxValue; i++) {
            for (uchar j = 0; j < symbolsPerElement; j++) {
                reverse[i][j] = symbolsList[sequence[j]];
                clear[i][j] = currentClear[j];
            }

            uchar j = symbolsPerElement - 1;
            while ((++sequence[j] == symbolsCount) && ((int) j > 0)) {
                sequence[j] = 0;
                currentClear[j] = i + 1;
                j--;
            }
        }
        delete[]currentClear;
        delete[]sequence;

        packLUT0 = new uchar[PACK_LUT_SIZE]();
        packLUT1 = new uchar[PACK_LUT_SIZE]();
        symbolsPerLUT1 = symbolsPerElement - SYMBOLS_PER_LUT_0;
        for (uint_max i = 0; i <= maxValue; i++) {
            uint16_t temp = 0;
            if (reverse[i][SYMBOLS_PER_LUT_0] == symbolsList[0]
                && (symbolsPerLUT1 == 1 || reverse[i][SYMBOLS_PER_LUT_0 + 1] == symbolsList[0])) {
                memcpy(&temp, reverse[i], SYMBOLS_PER_LUT_0);
                packLUT0[temp & PACK_MASK] = (uint8_t) i;
            }
            if (reverse[i][0] == symbolsList[0] && reverse[i][1] == symbolsList[0]) {
                memcpy(&temp, reverse[i] + SYMBOLS_PER_LUT_0, symbolsPerLUT1);
                packLUT1[temp & PACK_MASK] = (uint8_t) i;
            }
        }
    }

    uint8_t SymbolsPackingFacility::clearSuffix(const uint8_t value, uchar prefixLength) {
        return clear[value][prefixLength];
    }

    uint8_t SymbolsPackingFacility::getMaxValue() {
        return maxValue;
    }

    bool SymbolsPackingFacility::isCompatible(uchar symbolsPerElement, uchar symbolsCount) {
        return powuint(symbolsCount, symbolsPerElement) - 1 <= (uint8_t) - 1;
    }

    uchar SymbolsPackingFacility::maxSymbolsPerElement(uchar symbolsCount) {
        for (int i = 0; i < UCHAR_MAX; i++)
            if (!SymbolsPackingFacility::isCompatible(i + 1, symbolsCount))
                return i;
        return UCHAR_MAX;
    }

    uint8_t SymbolsPackingFacility::packPrefixSymbols(const char_pg* symbols, const uint_max length) {
        uint8_t value = 0;
        for (uchar j = 0; j < length; j++) {
            validateSymbol((uchar) symbols[j]);
            value = value * symbolsCount + symbolOrder[(uchar) symbols[j]];
        }

        return value;
    }

    uint_max SymbolsPackingFacility::packSequence(const char_pg* source, const uint_max length, uint8_t* dest) {
        uint_max i = 0;

        const char_pg* guard = source + length - symbolsPerElement;

        while (source <= guard) {
            dest[i++] = packSymbols(source);
            source += symbolsPerElement;
        }

        uint_max left = guard + symbolsPerElement - source;
        if (left > 0)
            dest[i++] = packSuffixSymbols(source, left);

        return i;
    }

    string SymbolsPackingFacility::packSequence(const char_pg *source, const uint_max length) {
        size_t packedLength = (length + symbolsPerElement - 1) / symbolsPerElement * sizeof(uint8_t);
        string tmp;
        tmp.resize(packedLength);
        this->packSequence(source, length, (uint8_t*) tmp.data());
        return tmp;
    }

    uint8_t SymbolsPackingFacility::packSuffixSymbols(const char_pg* symbols, const uint_max length) {
        uint8_t value = 0;
        for (uchar j = 0; j < symbolsPerElement; j++) {
            value *= symbolsCount;
            if (j < length) {
                validateSymbol((uchar) symbols[j]);
                value += symbolOrder[(uchar) symbols[j]];
            }
        }
        return value;
    }

    uint8_t SymbolsPackingFacility::packSymbols(const char_pg* symbols) {
        uint16_t temp0 = 0, temp1 = 0;
        memcpy(&temp0, symbols, SYMBOLS_PER_LUT_0);
        memcpy(&temp1, symbols + SYMBOLS_PER_LUT_0, symbolsPerLUT1);
        return packLUT0[temp0 & PACK_MASK] + packLUT1[temp1 & PACK_MASK];
    }

    const string SymbolsPackingFacility::reverseSequence(const uint8_t* sequence, const uint_max pos, const uint_max length) {
        string res;
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);
        res.append(reverse[sequence[i++]], reminder, this->symbolsPerElement - reminder);
        while (res.size() < length)
            res.append(reverse[sequence[i++]], symbolsPerElement);
        res.resize(length);
        return res;
    }

    void SymbolsPackingFacility::reverseSequence(const uint8_t* sequence, const uint_max pos, const uint_max length, string &res) {
        res.clear();
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);
        uint8_t value = sequence[i++];
        for (uchar j = reminder; j < symbolsPerElement; j++) {
            res.push_back(reverse[value][j]);
            if (res.size() < length)
                return;
        }
        res.append(reverse[sequence[i++]], reminder, this->symbolsPerElement - reminder);
        while (res.size() < length - symbolsPerElement)
            res.append(reverse[sequence[i++]], symbolsPerElement);
        value = sequence[i];
        for (uchar j = 0; res.size() < length; j++)
            res.push_back(reverse[value][j]);
    }

    void SymbolsPackingFacility::reverseSequence(const uint8_t* sequence, const uint_max pos, const uint_max length, char_pg* destPtr) {
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);

        char_pg* ptr = destPtr;
        const char_pg* endPtr = destPtr + length;
        uint8_t value = sequence[i++];
        for (uchar j = reminder; j < symbolsPerElement; j++) {
            *ptr++ = reverse[value][j];
            if (ptr == endPtr)
                return;
        }
        while(true) {
            value = sequence[i++];
            for (uchar j = 0; j < symbolsPerElement; j++) {
                *ptr++ = reverse[value][j];
                if (ptr == endPtr)
                    return;
            }
        }
   }

    char_pg SymbolsPackingFacility::reverseValue(uint8_t value, uchar position) {
        return reverse[value][position];
    }

    char_pg SymbolsPackingFacility::reverseValue(uint8_t* sequence, uint_read_len_max pos) {
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger<uint_max>(pos, this->symbolsPerElement, i);

        uint8_t value = sequence[i];
        return reverse[value][reminder];
    }

    string SymbolsPackingFacility::reverseValue(uint8_t value) {
        string res;
        res.resize(symbolsPerElement);
        for (uchar j = 0; j < symbolsPerElement; j++)
            res[j] = reverse[value][j];
        return res;
    }

    int SymbolsPackingFacility::compareSequences(uint8_t* lSeq, uint8_t* rSeq, const uint_max length) {
        uint_max i = length;
        while (i >= symbolsPerElement) {
            int cmp = (int) *lSeq++ - *rSeq++;
            if (cmp)
                return cmp;
            i -= symbolsPerElement;
        } 
        
        int j = 0;
        while (i--) {
            int cmp = (int) reverse[*lSeq][j] - reverse[*rSeq][j];
            if (cmp)
                return cmp;
            j++;
        }
        
        return 0;
    }

    int SymbolsPackingFacility::compareSequences(uint8_t* lSeq, uint8_t* rSeq, uint_max pos, uint_max length) {
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);

        lSeq += i;
        rSeq += i;
        for (uchar j = reminder; j < symbolsPerElement; j++) {
            int cmp = (int) reverse[*lSeq][j] - reverse[*rSeq][j];
            if (cmp)
                return cmp;
            if (--length == 0)
                return 0;
        }
        
        return compareSequences(lSeq + 1, rSeq + 1, length);
    }

    int SymbolsPackingFacility::compareSuffixWithPrefix(uint8_t* sufSeq, uint8_t* preSeq, uint_max sufPos, uint_max length) {
        uint_max i = divideBySmallInteger(sufPos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(sufPos, this->symbolsPerElement, i);
        
        sufSeq += i;
        if (reminder == 0)
            return compareSequences(sufSeq, preSeq, length);
        
        uchar sufIdx = reminder;
        uchar preIdx = 0;
        while (true) {
            int cmp = (int) reverse[*sufSeq][sufIdx] - reverse[*preSeq][preIdx];
            if (cmp)
                return cmp;
            if (--length == 0)
                return 0;
            if (++preIdx == symbolsPerElement) {
                preIdx = 0; preSeq++;
            } 
            if (++sufIdx == symbolsPerElement) {
                sufIdx = 0; sufSeq++;
            }
        }
    }

    int SymbolsPackingFacility::compareSequenceWithUnpacked(uint_ps_element_min *seq,
                                                                          const char *pattern,
                                                                          uint_read_len_max length) {
        uint_max i = length;
        while (i >= symbolsPerElement) {
            uint_ps_element_min pSeq = packSymbols(pattern);
            int cmp = (int) *seq++ - pSeq;
            if (cmp)
                return cmp;
            pattern += symbolsPerElement;
            i -= symbolsPerElement;
        }

        int j = 0;
        while (i--) {
            int cmp = (int) reverse[*seq][j] - *pattern++;
            if (cmp)
                return cmp;
            j++;
        }

        return 0;
    }

    uint8_t SymbolsPackingFacility::countSequenceMismatchesVsUnpacked(uint_ps_element_min *seq,
                                                                                    const char *pattern,
                                                                                    uint_read_len_max length,
                                                                                    uint8_t maxMismatches) {
        uint8_t res = 0;
        uint_max i = length;
        while (i >= symbolsPerElement) {
            uint_ps_element_min pSeq = packSymbols(pattern);
            if (*seq != pSeq) {
                for(uint8_t j = 0; j < symbolsPerElement; j++) {
                    if (reverse[*seq][j] != *pattern++)
                        if (res++ >= maxMismatches)
                            return UINT8_MAX;
                }
            } else
                pattern += symbolsPerElement;
            seq++;
            i -= symbolsPerElement;
        }

        int j = 0;
        while (i--) {
            if (reverse[*seq][j] != *pattern++) {
                if (res++ >= maxMismatches)
                    return UINT8_MAX;
            }
            j++;
        }

        return res;
    }

    void SymbolsPackingFacility::validateSymbol(uchar symbol) {
        if (symbolOrder[symbol] == -1) {
            fprintf(stdout, "Unexpected symbol '%c'. Packer supports: %s.\n", symbol, symbolsList);
            exit(EXIT_FAILURE);
        }
    }

    SymbolsPackingFacility *
    SymbolsPackingFacility::getInstance(ReadsSetProperties *properties, uchar symbolsPerElement) {
        if (symbolsPerElement == 3 &&
            strcmp(SymbolsPackingFacility::ACGTNPacker.symbolsList, properties->symbolsList) == 0)
            return &SymbolsPackingFacility::ACGTNPacker;
        if (symbolsPerElement == 4 &&
            strcmp(SymbolsPackingFacility::ACGTPacker.symbolsList, properties->symbolsList) == 0)
            return &SymbolsPackingFacility::ACGTPacker;
        return new SymbolsPackingFacility(properties, symbolsPerElement);
    }
}