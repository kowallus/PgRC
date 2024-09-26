#ifndef READSSETINTERFACE_H_INCLUDED
#define READSSETINTERFACE_H_INCLUDED

#include "../utils/helper.h"
#include "../pgrc/pg-config.h"
#include "ReadsSetBase.h"

using namespace PgReadsSet;

namespace PgReadsSet {

    template < typename uint_read_len, typename uint_reads_cnt >
    class ReadsSetInterface
    {
        public:

            virtual ~ReadsSetInterface() {};

            virtual bool isReadLengthConstant() = 0;
            virtual uint_read_len maxReadLength() = 0;

            virtual uint_reads_cnt readsCount() = 0;

            virtual const string getRead(uint_reads_cnt) = 0;
            virtual uint_read_len readLength(uint_reads_cnt) = 0;
    };

    class ConstantLengthReadsSetInterface: public ReadsSetInterface<uint_read_len_max, uint_reads_cnt_max>,
            public ReadsSetBase
    {
    public:

        virtual ~ConstantLengthReadsSetInterface() {};

        virtual char getReadSymbol(uint_reads_cnt_max index, uint_read_len_max pos) = 0;

        virtual int compareReadWithPattern(const uint_reads_cnt_max index, const char *pattern) = 0;

        virtual uint8_t countMismatchesVsPattern(uint_reads_cnt_max index, const char *pattern, const uint_read_len_max patLength,
                                           uint8_t mismatchesLimit) = 0;
        virtual void getRead(uint_reads_cnt_max index, char* res) = 0;
        virtual const string getRead(uint_reads_cnt_max index) = 0;
    };

    class SumOfConstantLengthReadsSets: public ConstantLengthReadsSetInterface
    {
    private:
        ConstantLengthReadsSetInterface* clrs1 = 0;
        ConstantLengthReadsSetInterface* clrs2 = 0;
        uint_reads_cnt_max idxBeg2 = 0;

        inline ConstantLengthReadsSetInterface *getCrls(uint_reads_cnt_max i) const
            { return (i < idxBeg2 ? clrs1 : clrs2); }
        inline uint_reads_cnt_max getClrsIdx(uint_reads_cnt_max i) const { return i < idxBeg2 ? i : i - idxBeg2; }

    public:
        SumOfConstantLengthReadsSets(ConstantLengthReadsSetInterface *clrs1, ConstantLengthReadsSetInterface *clrs2)
                : clrs1(clrs1), clrs2(clrs2), idxBeg2(clrs1->readsCount()) {}

        bool isReadLengthConstant() override { return true; }
        uint_read_len_max maxReadLength() override { return clrs1->maxReadLength(); }
        uint_reads_cnt_max readsCount() override { return clrs1->readsCount() + clrs2->readsCount(); }

        const string getRead(uint_reads_cnt_max i) override {
            return getCrls(i)->getRead(getClrsIdx(i));
        }

        void getRead(uint_reads_cnt_max i, char *res) override {
            getCrls(i)->getRead(getClrsIdx(i), res);
        }

        char getReadSymbol(uint_reads_cnt_max i, uint_read_len_max pos) override {
            return getCrls(i)->getReadSymbol(getClrsIdx(i), pos);
        }

        uint_read_len_max readLength(uint_reads_cnt_max i) override {
            return getCrls(i)->readLength(getClrsIdx(i));
        }

        int compareReadWithPattern(const uint_reads_cnt_max i, const char *pattern) override {
            return getCrls(i)->compareReadWithPattern(getClrsIdx(i), pattern);
        }

        uint8_t
        countMismatchesVsPattern(uint_reads_cnt_max i, const char *pattern, const uint_read_len_max patLength,
                                 uint8_t mismatchesLimit) override {
            return getCrls(i)->countMismatchesVsPattern(getClrsIdx(i), pattern, patLength, mismatchesLimit);
        }
    };

    typedef ReadsSetInterface< uint_read_len_std, uint_reads_cnt_std > StandardReadsSetInterface;
}

#endif // READSSETINTERFACE_H_INCLUDED
