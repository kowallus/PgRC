#ifndef PGTOOLS_DIVISIONREADSSETDECORATOR_H
#define PGTOOLS_DIVISIONREADSSETDECORATOR_H

#include "ReadsSetIterator.h"

namespace PgTools {

    template < typename uint_read_len >
    class QualityDividingReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate< uint_read_len >* coreIterator;
        int64_t counter = -1;
        double error_level;
        bool filterNReads;
        bool visitGoodReads;
        bool isQualityGood();
        bool containsN();

    public:

        QualityDividingReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level,
                                        bool filterNreads = false, bool visitGoodReads = true);

        virtual ~QualityDividingReadsSetIterator();

        uint64_t getReadOriginalIndex();

        bool moveNextVirtual();
        string getReadVirtual();
        string getQualityInfoVirtual();
        uint_read_len getReadLengthVirtual();
        void rewindVirtual();
    };

    template < typename uint_read_len >
    class DividedReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t counter = -1;
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

        bool moveNextVirtual();
        string getReadVirtual();
        string getQualityInfoVirtual();
        uint_read_len getReadLengthVirtual();
        void rewindVirtual();

        const vector<uint_reads_cnt_max> getVisitedIndexesMapping() override;
    };

}


#endif //PGTOOLS_DIVISIONREADSSETDECORATOR_H
