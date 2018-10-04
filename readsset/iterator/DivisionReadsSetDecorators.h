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
        bool visitGoodReads = true;
        bool isQualityGood();

    public:

        QualityDividingReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level,
                                        bool visitGoodReads = true);

        virtual ~QualityDividingReadsSetIterator();

        uint64_t getReadIndex();

        bool moveNextVirtual();
        string getReadVirtual();
        string getQualityInfoVirtual();
        uint_read_len getReadLengthVirtual();
        void rewindVirtual();
    };

    template < typename uint_read_len >
    class DividedReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate< uint_read_len > coreIterator;
        int counter = 0;
        double error_level;

    public:
        bool moveNextVirtual();
        string getReadVirtual();
        string getQualityInfoVirtual();
        uint_read_len getReadLengthVirtual();
        void rewindVirtual();
    };

}


#endif //PGTOOLS_DIVISIONREADSSETDECORATOR_H
