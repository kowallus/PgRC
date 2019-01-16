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
        bool filterNReads;
        bool visitGoodReads;
        bool isQualityGood();
        bool containsN();
        vector<uint_reads_cnt_max> indexesMapping;

    public:

        QualityDividingReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level,
                                        bool filterNreads = false, bool visitGoodReads = true);

        virtual ~QualityDividingReadsSetIterator();

        uint_reads_cnt_max getReadOriginalIndex();

        bool moveNext();
        string getRead();
        string getQualityInfo();
        uint_read_len getReadLength();
        void rewind();

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

        bool moveNext();
        string getRead();
        string getQualityInfo();
        uint_read_len getReadLength();
        void rewind();

        IndexesMapping* retainVisitedIndexesMapping() override;
    };

}


#endif //PGTOOLS_DIVISIONREADSSETDECORATOR_H
