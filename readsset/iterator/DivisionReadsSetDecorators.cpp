#include "DivisionReadsSetDecorators.h"

#include "../../helper.h"

using namespace PgSAHelpers;

namespace PgTools {

    template<typename uint_read_len>
    QualityDividingReadsSetIterator<uint_read_len>::QualityDividingReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level, bool visitGoodReads)
            :coreIterator(coreIterator), error_level(error_level), visitGoodReads(visitGoodReads) {}

    template<typename uint_read_len>
    QualityDividingReadsSetIterator<uint_read_len>::~QualityDividingReadsSetIterator() {}

    template<typename uint_read_len>
    bool QualityDividingReadsSetIterator<uint_read_len>::moveNextVirtual() {
        while (coreIterator->moveNextVirtual()) {
            counter++;
            if (isQualityGood() == visitGoodReads)
                return true;
        }
        return false;
    }

    template<typename uint_read_len>
    bool QualityDividingReadsSetIterator<uint_read_len>::isQualityGood() {
        return 1 - qualityScore2correctProb(getQualityInfoVirtual()) < error_level;
    }

    template<typename uint_read_len>
    string QualityDividingReadsSetIterator<uint_read_len>::getReadVirtual() {
        return coreIterator->getReadVirtual();
    }

    template<typename uint_read_len>
    string QualityDividingReadsSetIterator<uint_read_len>::getQualityInfoVirtual() {
        return coreIterator->getQualityInfoVirtual();
    }

    template<typename uint_read_len>
    uint_read_len QualityDividingReadsSetIterator<uint_read_len>::getReadLengthVirtual() {
        return coreIterator->getReadLengthVirtual();
    }

    template<typename uint_read_len>
    void QualityDividingReadsSetIterator<uint_read_len>::rewindVirtual() {
        counter = -1;
        return coreIterator->rewindVirtual();
    }

    template<typename uint_read_len>
    uint64_t QualityDividingReadsSetIterator<uint_read_len>::getReadIndex() {
        return counter;
    }

    template class QualityDividingReadsSetIterator<uint_read_len_min>;
    template class QualityDividingReadsSetIterator<uint_read_len_std>;
}
