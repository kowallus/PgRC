#include "DivisionReadsSetDecorators.h"

#include "../../utils/helper.h"

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
        counter++;
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
        coreIterator->rewindVirtual();
    }

    template<typename uint_read_len>
    uint64_t QualityDividingReadsSetIterator<uint_read_len>::getReadOriginalIndex() {
        return counter;
    }

    template<typename uint_read_len>
    DividedReadsSetIterator<uint_read_len>::DividedReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, std::istream *divSource, bool visitComplement)
            :coreIterator(coreIterator), divSource(divSource), visitComplement(visitComplement) {
        readValue(*divSource, currentDivIdx);
    }

    template<typename uint_read_len>
    bool DividedReadsSetIterator<uint_read_len>::moveNextVirtual() {
        while (coreIterator->moveNextVirtual()) {
            counter++;
            if (visitComplement) {
                if (counter != currentDivIdx)
                    return true;
                else
                    readValue(*divSource, currentDivIdx);
            } else {
                if (counter == currentDivIdx) {
                    readValue(*divSource, currentDivIdx);
                    return true;
                }
            }
        }
        counter++;
        return false;
    }

    template<typename uint_read_len>
    string DividedReadsSetIterator<uint_read_len>::getReadVirtual() {
        return coreIterator->getReadVirtual();
    }

    template<typename uint_read_len>
    string DividedReadsSetIterator<uint_read_len>::getQualityInfoVirtual() {
        return coreIterator->getQualityInfoVirtual();
    }

    template<typename uint_read_len>
    uint_read_len DividedReadsSetIterator<uint_read_len>::getReadLengthVirtual() {
        return coreIterator->getReadLengthVirtual();
    }

    template<typename uint_read_len>
    void DividedReadsSetIterator<uint_read_len>::rewindVirtual() {
        counter = -1;
        divSource->clear();
        divSource->seekg(0);
        readValue(*divSource, currentDivIdx);
        coreIterator->rewindVirtual();
    }

    template class QualityDividingReadsSetIterator<uint_read_len_min>;
    template class QualityDividingReadsSetIterator<uint_read_len_std>;

    template class DividedReadsSetIterator<uint_read_len_min>;
    template class DividedReadsSetIterator<uint_read_len_std>;
}
