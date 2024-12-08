#ifndef PGTOOLS_DIVISIONREADSSETDECORATOR_H
#define PGTOOLS_DIVISIONREADSSETDECORATOR_H

#include "ReadsSetIterator.h"

namespace PgTools {

    template < typename uint_read_len >
    class QualityDividingReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate< uint_read_len >* coreIterator;
        int64_t allCounter = -1;
        double error_level;
        int suffix_pos;
        bool suffix_simplified_mode;

    public:

        QualityDividingReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level,
                                        bool suffix_simplified_mode = true, double read_length = 0);

        ~QualityDividingReadsSetIterator() override;

        uint_reads_cnt_max getReadOriginalIndex();

        bool moveNext() override;
        string& getRead() override;
        string& getQualityInfo() override;
        uint_read_len getReadLength() override;
        void rewind() override;

        bool isQualityHigh();
        bool containsN();

        IndexesMapping* retainVisitedIndexesMapping() override;
    };

    template < typename uint_read_len >
    class DividedReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t allCounter = -1;
        uint64_t currentDivIdx;
        std::istream* divSource;
        bool visitComplement;
        bool ignoreNReads;
        bool ignoreNoNReads;
        bool plainTextReadMode = false;

        vector<uint_reads_cnt_max> indexesMapping;
        bool isIgnored();
    public:
        DividedReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, std::istream* divSource,
                bool visitComplement = false, bool ignoreNReads = false, bool ignoreNoNReads = false);

        bool moveNext() override;
        string& getRead() override;
        string& getQualityInfo() override;
        uint_read_len getReadLength() override;
        void rewind() override;

        IndexesMapping* retainVisitedIndexesMapping() override;
    };

}


#endif //PGTOOLS_DIVISIONREADSSETDECORATOR_H
