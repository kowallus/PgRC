#include "DivisionReadsSetDecorators.h"

#include "../../utils/helper.h"

using namespace PgSAHelpers;

namespace PgTools {

    template<typename uint_read_len>
    QualityDividingReadsSetIterator<uint_read_len>::QualityDividingReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level,
                bool filterNReads, bool visitGoodReads)
            :coreIterator(coreIterator), error_level(error_level),
             filterNReads(filterNReads), visitGoodReads(visitGoodReads) {}

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
        return (1 - qualityScore2correctProb(getQualityInfoVirtual()) < error_level) &&
                (!filterNReads || coreIterator->getReadVirtual().find('N') == string::npos);
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
        return (uint64_t) counter;
    }

    template<typename uint_read_len>
    bool QualityDividingReadsSetIterator<uint_read_len>::containsN() {
        return coreIterator->getReadVirtual().find('N') != string::npos;
    }

    template<typename uint_read_len>
    DividedReadsSetIterator<uint_read_len>::DividedReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, std::istream *divSource, bool visitComplement,
            bool ignoreNReads, bool ignoreNoNReads)
            :coreIterator(coreIterator), divSource(divSource), visitComplement(visitComplement),
            ignoreNReads(ignoreNReads), ignoreNoNReads(ignoreNoNReads) {
        plainTextReadMode = readReadMode(*divSource);
        readValue(*divSource, currentDivIdx, plainTextReadMode);
    }

    template<typename uint_read_len>
    bool DividedReadsSetIterator<uint_read_len>::moveNextVirtual() {
        if (counter >= 0)
            indexesMapping.push_back(counter);
        while (coreIterator->moveNextVirtual()) {
            counter++;
            if (visitComplement) {
                if ((counter != currentDivIdx) && !isIgnored())
                    return true;

                readValue(*divSource, currentDivIdx, plainTextReadMode);
            } else {
                if (counter == currentDivIdx) {
                    readValue(*divSource, currentDivIdx, plainTextReadMode);
                    if (!isIgnored())
                        return true;
                }
            }
        }
        indexesMapping.push_back(++counter);
        return false;
    }

    template<typename uint_read_len>
    const vector<uint_reads_cnt_max> DividedReadsSetIterator<uint_read_len>::getVisitedIndexesMapping() {
        return indexesMapping;
    }

    template<typename uint_read_len>
    bool DividedReadsSetIterator<uint_read_len>::isIgnored() {
        return (ignoreNReads && coreIterator->getReadVirtual().find('N') != string::npos)
            || (ignoreNoNReads && coreIterator->getReadVirtual().find('N') == string::npos);
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
        indexesMapping.clear();
        divSource->clear();
        divSource->seekg(0);
        readReadMode(*divSource);
        readValue(*divSource, currentDivIdx, plainTextReadMode);
        coreIterator->rewindVirtual();
    }

    template class QualityDividingReadsSetIterator<uint_read_len_min>;
    template class QualityDividingReadsSetIterator<uint_read_len_std>;

    template class DividedReadsSetIterator<uint_read_len_min>;
    template class DividedReadsSetIterator<uint_read_len_std>;
}
