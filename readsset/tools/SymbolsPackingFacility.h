/* 
 * File:   SymbolsPackingFacility.h
 * Author: coach
 *
 * Created on 12 luty 2014, 11:43
 */

#ifndef SYMBOLSPACKINGFACILITY_H
#define	SYMBOLSPACKINGFACILITY_H

#include "../../utils/helper.h"
#include "../../pgsaconfig.h"
#include "../ReadsSetBase.h"

using namespace PgSAHelpers;

namespace PgSAIndex {
    
    template<typename uint_element>
    class SymbolsPackingFacility {
        private:
            
            uint_max maxValue;
            const uint_symbols_cnt symbolsCount;
            const uchar symbolsPerElement;
            
            char symbolsList[UCHAR_MAX] = {};
            int symbolOrder[UCHAR_MAX] = {};

            // for a given value reverse[value][pos] returns character at the given position
            char_pg** reverse;
            char_pg* reverseFlat;
            
            // for a given value clear[value][prefixLenght] returns value of first prefixLength symbols followed by 0 symbols.
            uint_element** clear;
            uint_element* clearFlat;
            
            void buildReverseAndClearIndexes();

            inline void validateSymbol(uchar symbol);

        public:
            SymbolsPackingFacility(ReadsSetProperties* readsSetProperties, uchar symbolsPerElement);
            
            ~SymbolsPackingFacility();
            
            uint_element getMaxValue();
            
            // value should be not greater then maxValue
            char_pg reverseValue(uint_element value, uchar position);
            char_pg reverseValue(uint_element* sequence, uint_read_len_max position);

            // value should be not greater then maxValue
            string reverseValue(uint_element value);
            
            const string reverseSequence(const uint_element* sequence, const uint_max pos, const uint_max length);
            void reverseSequence(const uint_element* sequence, const uint_max pos, const uint_max length, string& res);
            void reverseSequence(const uint_element* sequence, const uint_max pos, const uint_max length, char_pg* destPtr);
          
            // sequence should consist of at least symbolsPerElement symbols; result should be not greater then maxValue
            uint_element packSymbols(const char_pg* symbols);
            
            // result should be not greater then maxValue
            uint_element packSuffixSymbols(const char_pg* symbols, const uint_max length);
            
            // result should be not greater then maxValue
            uint_element packPrefixSymbols(const char_pg* symbols, const uint_max length);
            
            uint_max packSequence(const char_pg* source, const uint_max length, uint_element* dest);

            uint_element clearSuffix(const uint_element value, uchar prefixLength);
            
            static bool isCompatibile(uchar symbolsPerElement, uchar symbolsCount);
            
            static uchar maxSymbolsPerElement(uchar symbolsCount);
          
            int compareSequences(uint_element* lSeq, uint_element* rSeq, const uint_max length);
            int compareSequences(uint_element* lSeq, uint_element* rSeq, uint_max pos, uint_max length);
            int compareSuffixWithPrefix(uint_element* sufSeq, uint_element* preSeq, uint_max sufPos, uint_max length);

        int compareSequenceWithUnpacked(uint_ps_element_min *seq, const char *pattern, uint_read_len_max length);

        uint8_t countSequenceMismatchesVsUnpacked(uint_ps_element_min *seq, const char *pattern, uint_read_len_max length,
                                              uint8_t maxMismatches);
    };
    
}

#endif	/* SYMBOLSPACKINGFACILITY_H */

