#ifndef SYMBOLSPACKINGFACILITY_H
#define	SYMBOLSPACKINGFACILITY_H

#include "../utils/helper.h"
#include "../pgrc/pg-config.h"
#include "../readsset/ReadsSetBase.h"

using namespace PgHelpers;

namespace PgIndex {

    class SymbolsPackingFacility {
        private:

            const static uint16_t PACK_LUT_SIZE = 1 << 11;
            const static uint16_t PACK_MASK = PACK_LUT_SIZE - 1;

            uint_max maxValue;
            const uint_symbols_cnt symbolsCount;
            const uchar symbolsPerElement;
            
            char symbolsList[UCHAR_MAX] = {};
            int symbolOrder[UCHAR_MAX] = {};

            uint8_t* packLUT0;
            uint8_t* packLUT1;
            const uchar SYMBOLS_PER_LUT_0 = 2;
            uchar symbolsPerLUT1;

            // for a given value reverse[value][pos] returns character at the given position
            char_pg** reverse;
            char_pg* reverseFlat;
            
            // for a given value clear[value][prefixLenght] returns value of first prefixLength symbols followed by 0 symbols.
            uint8_t** clear;
            uint8_t* clearFlat;

            const bool globallyManaged = false;

            void buildReversePackAndClearIndexes();

            inline void validateSymbol(uchar symbol);

        public:
            SymbolsPackingFacility(ReadsSetProperties* readsSetProperties, uchar symbolsPerElement);
            SymbolsPackingFacility(const vector<char> symbolsList, bool isGloballyManaged = false);
            
            ~SymbolsPackingFacility();

            bool isGloballyManaged() { return globallyManaged; };

            uint8_t getMaxValue();
            
            // value should be not greater then maxValue
            char_pg reverseValue(uint8_t value, uchar position);
            char_pg reverseValue(uint8_t* sequence, uint_read_len_max position);

            // value should be not greater then maxValue
            string reverseValue(uint8_t value);
            
            const string reverseSequence(const uint8_t* sequence, const uint_max pos, const uint_max length);
            void reverseSequence(const uint8_t* sequence, const uint_max pos, const uint_max length, string& res);
            void reverseSequence(const uint8_t* sequence, const uint_max pos, const uint_max length, char_pg* destPtr);
          
            // sequence should consist of at least symbolsPerElement symbols; result should be not greater then maxValue
            inline uint8_t packSymbols(const char_pg* symbols);
            
            // result should be not greater then maxValue
            uint8_t packSuffixSymbols(const char_pg* symbols, const uint_max length);
            
            // result should be not greater then maxValue
            uint8_t packPrefixSymbols(const char_pg* symbols, const uint_max length);
            
            uint_max packSequence(const char_pg* source, const uint_max length, uint8_t* dest);
            string packSequence(const char_pg* source, const uint_max length);

            uint8_t clearSuffix(const uint8_t value, uchar prefixLength);
            
            static bool isCompatible(uchar symbolsPerElement, uchar symbolsCount);
            
            static uchar maxSymbolsPerElement(uchar symbolsCount);
          
            int compareSequences(uint8_t* lSeq, uint8_t* rSeq, const uint_max length);
            int compareSequences(uint8_t* lSeq, uint8_t* rSeq, uint_max pos, uint_max length);
            int compareSuffixWithPrefix(uint8_t* sufSeq, uint8_t* preSeq, uint_max sufPos, uint_max length);

            int compareSequenceWithUnpacked(uint_ps_element_min *seq, const char *pattern, uint_read_len_max length);

            uint8_t countSequenceMismatchesVsUnpacked(uint_ps_element_min *seq, const char *pattern, uint_read_len_max length,
                                              uint8_t maxMismatches);

            static SymbolsPackingFacility ACGTPacker, ACGTNPacker;

            static SymbolsPackingFacility *getInstance(ReadsSetProperties *properties, uchar symbolsPerElement);
    };
    
}

#endif	/* SYMBOLSPACKINGFACILITY_H */

