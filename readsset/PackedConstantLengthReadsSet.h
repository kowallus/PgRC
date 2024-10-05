#ifndef PACKEDREADSSET_H_INCLUDED
#define PACKEDREADSSET_H_INCLUDED

#include <algorithm>
#include "ReadsSetBase.h"
#include "iterator/ReadsSetIterator.h"
#include "ReadsSetInterface.h"
#include "../coders/SymbolsPackingFacility.h"

using namespace PgIndex;

namespace PgReadsSet {

    class PackedConstantLengthReadsSet: public ConstantLengthReadsSetInterface
    {
        private:
            vector<uint_ps_element_min> packedReads;
            uchar packedLength;

        public:

            SymbolsPackingFacility* sPacker;

            PackedConstantLengthReadsSet(uint_read_len_max readLength, const char* symbolsList, uint8_t symbolsCount);
            ~PackedConstantLengthReadsSet() override;

            void reserve(uint_reads_cnt_max readsCount);
            void resize(uint_reads_cnt_max readsCount);
            void addRead(const char* read, uint_read_len_max readLength);
            void copyRead(uint_reads_cnt_max srcIdx, uint_reads_cnt_max destIdx, uint_reads_cnt_max n = 1);
            void copyPackedRead(const uint_ps_element_min *packedSequence, uint_reads_cnt_max destIdx,
                    uint_reads_cnt_max n = 1);

            inline uint_read_len_max minReadLength() { return properties->minReadLength; };
            inline uint_read_len_max maxReadLength() override { return properties->maxReadLength; };
            inline uint_reads_cnt_max readsCount() override { return properties->readsCount; };

            inline bool isReadLengthConstant() override { return properties->constantReadLength; };

            inline const uint_ps_element_min* getPackedRead(uint_reads_cnt_max i) { return packedReads.data() + i * (size_t) packedLength;};
            inline const string getReadPrefix(uint_reads_cnt_max i, uint_read_len_max skipSuffix) { return sPacker->reverseSequence(packedReads.data() + i * (size_t) packedLength, 0, readLength(i) - skipSuffix);};
            inline void getReadSuffix(char *destPtr, uint_reads_cnt_max i, uint_read_len_max suffixPos) { sPacker->reverseSequence(packedReads.data() + i * (size_t) packedLength, suffixPos, readLength(i) - suffixPos, destPtr);};
            inline const string getRead(uint_reads_cnt_max i) override { return sPacker->reverseSequence(packedReads.data() + (size_t) packedLength * i, 0, readLength(i));};
            inline void getRead(uint_reads_cnt_max i, string &res) { sPacker->reverseSequence(packedReads.data() + (size_t) packedLength * i, 0, readLength(i), res);};
            inline void getRead(uint_reads_cnt_max i, char_pg* res) override { sPacker->reverseSequence(packedReads.data() + (size_t) packedLength * i, 0, readLength(i), res);};
            inline char getReadSymbol(uint_reads_cnt_max i, uint_read_len_max pos) override { return sPacker->reverseValue(packedReads.data() + i * (size_t) packedLength, pos); };
            inline uint_read_len_max readLength(uint_reads_cnt_max i) override { return maxReadLength(); };

            int comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx);
            int comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx, uint_read_len_max offset);
            int compareSuffixWithPrefix(uint_reads_cnt_max sufIdx, uint_reads_cnt_max preIdx, uint_read_len_max sufOffset);

            int compareReadWithPattern(const uint_reads_cnt_max i, const char *pattern) override;
            int compareReadWithPattern(const uint_reads_cnt_max i, const char *pattern, int length);
            uint8_t countMismatchesVsPattern(uint_reads_cnt_max i, const char *pattern, uint_read_len_max length, uint8_t maxMismatches) override;

            uint_read_len_max maxReadLengthVirtual() { return maxReadLength(); };
            uint_reads_cnt_max readsCountVirtual() { return readsCount(); };
            bool isReadLengthConstantVirtual() { return isReadLengthConstant(); };
            const string getReadVirtual(uint_reads_cnt_max i) { return getRead(i); };
            uint_read_len_max readLengthVirtual(uint_reads_cnt_max i) { return readLength(i); };

            void printout() {
                properties->printout();
            };

            template<class ReadsSourceIterator>
            static PackedConstantLengthReadsSet* loadReadsSet(ReadsSourceIterator* readsIterator, ReadsSetProperties* properties = nullptr);

    };
}

#endif // PACKEDREADSSET_H_INCLUDED
