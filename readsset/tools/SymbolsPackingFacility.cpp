#include "SymbolsPackingFacility.h"

#include <cassert>

namespace PgSAIndex {

    static const vector<char> binaryCodes { 0, 1 };
    static const vector<char> trenaryCodes { 0, 1, 2 };
    static const vector<char> quaternaryCodes { 0, 1, 2, 3 };

    static const vector<char> acgtSymbols { 'A', 'C', 'G', 'T' };
    static const vector<char> acgtnSymbols { 'A', 'C', 'G', 'N', 'T' };

    template<typename uint_element>
    SymbolsPackingFacility<uint_element> SymbolsPackingFacility<uint_element>::ACGTPacker(acgtSymbols, true);
    template<typename uint_element>
    SymbolsPackingFacility<uint_element> SymbolsPackingFacility<uint_element>::ACGTNPacker(acgtnSymbols, true);

    template<typename uint_element>
    SymbolsPackingFacility<uint_element>::SymbolsPackingFacility(const vector<char> symbolsList, bool isGloballyManaged):
         symbolsCount(symbolsList.size()),
         globallyManaged(isGloballyManaged),
         symbolsPerElement(SymbolsPackingFacility<uint_element>::maxSymbolsPerElement(symbolsCount)) {
        uint_max combinationCount = powuint(symbolsCount, symbolsPerElement);
        maxValue = combinationCount - 1;
        if (maxValue > (int) (uint_element) - 1)
            cout << "ERROR in symbols packaging: max value for type: " << (int) (uint_element) - 1 << " while max " << " \n";

        reverse = new char_pg*[combinationCount];
        reverseFlat = new char_pg[combinationCount * symbolsPerElement];

        clear = new uint_element*[combinationCount];
        clearFlat = new uint_element[combinationCount * symbolsPerElement];

        std::copy(symbolsList.begin(), symbolsList.end(), std::begin(this->symbolsList));
        memset(symbolOrder, -1, UCHAR_MAX);
        for (uint_symbols_cnt i = 0; i < symbolsCount; i++)
            symbolOrder[(unsigned char) symbolsList[(unsigned char) i]] = i;

        buildReversePackAndClearIndexes();
    }

    template<typename uint_element>
    SymbolsPackingFacility<uint_element>::SymbolsPackingFacility(ReadsSetProperties* readsSetProperties, uchar symbolsPerElement)
    : symbolsCount(readsSetProperties->symbolsCount),
    symbolsPerElement(symbolsPerElement) {
        uint_max combinationCount = powuint(symbolsCount, symbolsPerElement);
        maxValue = combinationCount - 1;
        if (maxValue > (int) (uint_element) - 1)
            cout << "ERROR in symbols packaging: max value for type: " << (int) (uint_element) - 1 << " while max " << " \n";

        reverse = new char_pg*[combinationCount];
        reverseFlat = new char_pg[combinationCount * symbolsPerElement];

        clear = new uint_element*[combinationCount];
        clearFlat = new uint_element[combinationCount * symbolsPerElement];

        std::copy(std::begin(readsSetProperties->symbolsList), std::end(readsSetProperties->symbolsList), std::begin(symbolsList));
        std::copy(std::begin(readsSetProperties->symbolOrder), std::end(readsSetProperties->symbolOrder), std::begin(symbolOrder));

        buildReversePackAndClearIndexes();
    }

    template<typename uint_element>
    SymbolsPackingFacility<uint_element>::~SymbolsPackingFacility() {
        delete[]reverse;
        delete[]reverseFlat;
        delete[]clear;
        delete[]clearFlat;
        delete[]packLUT;
    }

    template<typename uint_element>
    void SymbolsPackingFacility<uint_element>::buildReversePackAndClearIndexes() {
        char_pg* rPtr = reverseFlat;
        uint_element* cPtr = clearFlat;
        for (uint_max i = 0; i <= maxValue; i++) {
            reverse[i] = rPtr;
            rPtr += symbolsPerElement;
            clear[i] = cPtr;
            cPtr += symbolsPerElement;
        }

        uint_element* currentClear = new uint_element[symbolsPerElement]();
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

        packLUT = new uint_element[PACK_LUT_SIZE];
        for (uint_max i = 0; i <= maxValue; i++) {
            uint32_t temp = 0;
            memcpy(&temp, reverse[i], symbolsPerElement);
            packLUT[temp & PACK_MASK] = i;
        }
    }

    template<typename uint_element>
    uint_element SymbolsPackingFacility<uint_element>::clearSuffix(const uint_element value, uchar prefixLength) {
        return clear[value][prefixLength];
    }

    template<typename uint_element>
    uint_element SymbolsPackingFacility<uint_element>::getMaxValue() {
        return maxValue;
    }

    template<typename uint_element>
    bool SymbolsPackingFacility<uint_element>::isCompatibile(uchar symbolsPerElement, uchar symbolsCount) {
        return powuint(symbolsCount, symbolsPerElement) - 1 <= (uint_element) - 1;
    }

    template<typename uint_element>
    uchar SymbolsPackingFacility<uint_element>::maxSymbolsPerElement(uchar symbolsCount) {
        for (int i = 0; i < UCHAR_MAX; i++)
            if (!SymbolsPackingFacility<uint_element>::isCompatibile(i + 1, symbolsCount))
                return i;
        return UCHAR_MAX;
    }

    template<typename uint_element>
    uint_element SymbolsPackingFacility<uint_element>::packPrefixSymbols(const char_pg* symbols, const uint_max length) {
        uint_element value = 0;
        for (uchar j = 0; j < length; j++) {
            validateSymbol((uchar) symbols[j]);
            value = value * symbolsCount + symbolOrder[(uchar) symbols[j]];
        }

        return value;
    }

    template<typename uint_element>
    uint_max SymbolsPackingFacility<uint_element>::packSequence(const char_pg* source, const uint_max length, uint_element* dest) {
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

    template<typename uint_element>
    string SymbolsPackingFacility<uint_element>::packSequence(const char_pg *source, const uint_max length) {
        size_t packedLength = (length + symbolsPerElement - 1) / symbolsPerElement * sizeof(uint_element);
        string tmp;
        tmp.resize(packedLength);
        this->packSequence(source, length, (uint_element*) tmp.data());
        return tmp;
    }

    template<typename uint_element>
    uint_element SymbolsPackingFacility<uint_element>::packSuffixSymbols(const char_pg* symbols, const uint_max length) {
        uint_element value = 0;
        for (uchar j = 0; j < symbolsPerElement; j++) {
            value *= symbolsCount;
            if (j < length) {
                validateSymbol((uchar) symbols[j]);
                value += symbolOrder[(uchar) symbols[j]];
            }
        }
        return value;
    }

    template<typename uint_element>
    uint_element SymbolsPackingFacility<uint_element>::packSymbols(const char_pg* symbols) {
        uint32_t temp = 0;
        memcpy(&temp, symbols, symbolsPerElement);
        return packLUT[temp & PACK_MASK];
    }

    template<typename uint_element>
    inline const string SymbolsPackingFacility<uint_element>::reverseSequence(const uint_element* sequence, const uint_max pos, const uint_max length) {
        string res;
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);
        res.append(reverse[sequence[i++]], reminder, this->symbolsPerElement - reminder);
        while (res.size() < length)
            res.append(reverse[sequence[i++]], symbolsPerElement);
        res.resize(length);
        return res;
    }

    template<typename uint_element>
    inline void SymbolsPackingFacility<uint_element>::reverseSequence(const uint_element* sequence, const uint_max pos, const uint_max length, string &res) {
        res.clear();
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);
        res.append(reverse[sequence[i++]], reminder, this->symbolsPerElement - reminder);
        while (res.size() < length - symbolsPerElement)
            res.append(reverse[sequence[i++]], symbolsPerElement);
        uint_element value = sequence[i];
        for (uchar j = 0; res.size() < length; j++)
            res.push_back(reverse[value][j]);
    }

    template<typename uint_element>
    inline void SymbolsPackingFacility<uint_element>::reverseSequence(const uint_element* sequence, const uint_max pos, const uint_max length, char_pg* destPtr) {
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);

        char_pg* ptr = destPtr;
        const char_pg* endPtr = destPtr + length;
        uint_element value = sequence[i++];
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

    template<typename uint_element>
    inline char_pg SymbolsPackingFacility<uint_element>::reverseValue(uint_element value, uchar position) {
        return reverse[value][position];
    }

    template<typename uint_element>
    inline char_pg SymbolsPackingFacility<uint_element>::reverseValue(uint_element* sequence, uint_read_len_max pos) {
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger<uint_max>(pos, this->symbolsPerElement, i);

        uint_element value = sequence[i++];
        return reverse[value][reminder];
    }

    template<typename uint_element>
    inline string SymbolsPackingFacility<uint_element>::reverseValue(uint_element value) {
        string res;
        res.resize(symbolsPerElement);
        for (uchar j = 0; j < symbolsPerElement; j++)
            res[j] = reverse[value][j];
        return res;
    }

    template<typename uint_element>
    inline int SymbolsPackingFacility<uint_element>::compareSequences(uint_element* lSeq, uint_element* rSeq, const uint_max length) {
        uint_max i = length;
        while (i >= symbolsPerElement) {
            if (*lSeq > *rSeq)
                return 1;
            if (*lSeq++ < *rSeq++)
                return -1;
            i -= symbolsPerElement;
        } 
        
        int j = 0;
        while (i--) {
            if (reverse[*lSeq][j] > reverse[*rSeq][j])
                return 1;
            if (reverse[*lSeq][j] < reverse[*rSeq][j])
                return -1;
            j++;
        }
        
        return 0;
    }

    template<typename uint_element>
    inline int SymbolsPackingFacility<uint_element>::compareSequences(uint_element* lSeq, uint_element* rSeq, uint_max pos, uint_max length) {
        uint_max i = divideBySmallInteger(pos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);

        lSeq += i;
        rSeq += i;
        for (uchar j = reminder; j < symbolsPerElement; j++) {
            if (reverse[*lSeq][j] > reverse[*rSeq][j])
                return 1;
            if (reverse[*lSeq][j] < reverse[*rSeq][j])
                return -1;
            if (--length == 0)
                return 0;
        }
        
        return compareSequences(lSeq + 1, rSeq + 1, length);
    }

    template<typename uint_element>
    inline int SymbolsPackingFacility<uint_element>::compareSuffixWithPrefix(uint_element* sufSeq, uint_element* preSeq, uint_max sufPos, uint_max length) {
        uint_max i = divideBySmallInteger(sufPos, symbolsPerElement);
        uint_max reminder = moduloBySmallInteger(sufPos, this->symbolsPerElement, i);
        
        sufSeq += i;
        if (reminder == 0)
            return compareSequences(sufSeq, preSeq, length);
        
        uchar sufIdx = reminder;
        uchar preIdx = 0;
        while (true) {
            if (reverse[*sufSeq][sufIdx] > reverse[*preSeq][preIdx])
                return 1;
            if (reverse[*sufSeq][sufIdx] < reverse[*preSeq][preIdx])
                return -1;
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

    template<typename uint_element>
    int SymbolsPackingFacility<uint_element>::compareSequenceWithUnpacked(uint_ps_element_min *seq,
                                                                          const char *pattern,
                                                                          uint_read_len_max length) {
        uint_max i = length;
        while (i >= symbolsPerElement) {
            uint_ps_element_min pSeq = packSymbols(pattern);
            if (*seq > pSeq)
                return 1;
            if (*seq++ < pSeq)
                return -1;
            pattern += symbolsPerElement;
            i -= symbolsPerElement;
        }

        int j = 0;
        while (i--) {
            if (reverse[*seq][j] > *pattern)
                return 1;
            if (reverse[*seq][j] < *pattern++)
                return -1;
            j++;
        }

        return 0;
    }

    template<typename uint_element>
    uint8_t SymbolsPackingFacility<uint_element>::countSequenceMismatchesVsUnpacked(uint_ps_element_min *seq,
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

    template<typename uint_element>
    void SymbolsPackingFacility<uint_element>::validateSymbol(uchar symbol) {
        if (symbolOrder[symbol] == -1) {
            fprintf(stdout, "Unexpected symbol '%c'. Packer supports: %s.\n", symbol, symbolsList);
            exit(EXIT_FAILURE);
        }
    }

    template<typename uint_element>
    SymbolsPackingFacility<uint_element> *
    SymbolsPackingFacility<uint_element>::getInstance(ReadsSetProperties *properties, uchar symbolsPerElement) {
        if (symbolsPerElement == 3 &&
            strcmp(SymbolsPackingFacility<uint_element>::ACGTNPacker.symbolsList, properties->symbolsList))
            return &SymbolsPackingFacility<uint_element>::ACGTNPacker;
        if (symbolsPerElement == 4 &&
            strcmp(SymbolsPackingFacility<uint_element>::ACGTPacker.symbolsList, properties->symbolsList))
            return &SymbolsPackingFacility<uint_element>::ACGTPacker;
        return new SymbolsPackingFacility<uint_element>(properties, symbolsPerElement);
    }

    template class SymbolsPackingFacility<uint_ps_element_min>;
}