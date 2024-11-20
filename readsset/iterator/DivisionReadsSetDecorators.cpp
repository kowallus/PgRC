#include "DivisionReadsSetDecorators.h"

#include "../../utils/helper.h"

using namespace PgHelpers;

namespace PgTools {

    template<typename uint_read_len>
    QualityDividingReadsSetIterator<uint_read_len>::QualityDividingReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, double error_level, bool suffix_simplified_mode, double read_length)
            :coreIterator(coreIterator), error_level(error_level), suffix_simplified_mode(suffix_simplified_mode) {
                suffix_pos = read_length * (1 - error_level);
            }

    template<typename uint_read_len>
    QualityDividingReadsSetIterator<uint_read_len>::~QualityDividingReadsSetIterator() {}

    template<typename uint_read_len>
    bool QualityDividingReadsSetIterator<uint_read_len>::moveNext() {
        while (coreIterator->moveNext()) {
            allCounter++;
            return true;
        }
        allCounter++;
        return false;
    }

    template<typename uint_read_len>
    bool QualityDividingReadsSetIterator<uint_read_len>::isQualityHigh() {
        string& quality = getQualityInfo();
        if (suffix_simplified_mode) {
            return quality[suffix_pos] > '#';
        } else {
            double q = qualityScore2correctProbArithAvg(quality);
            return (1 - q <= error_level);
        }
    }

    template<typename uint_read_len>
    string& QualityDividingReadsSetIterator<uint_read_len>::getRead() {
        return coreIterator->getRead();
    }

    template<typename uint_read_len>
    string& QualityDividingReadsSetIterator<uint_read_len>::getQualityInfo() {
        return coreIterator->getQualityInfo();
    }

    template<typename uint_read_len>
    uint_read_len QualityDividingReadsSetIterator<uint_read_len>::getReadLength() {
        return coreIterator->getReadLength();
    }

    template<typename uint_read_len>
    void QualityDividingReadsSetIterator<uint_read_len>::rewind() {
        allCounter = -1;
        coreIterator->rewind();
    }

    template<typename uint_read_len>
    uint_reads_cnt_max QualityDividingReadsSetIterator<uint_read_len>::getReadOriginalIndex() {
        return allCounter;
    }

    template<typename uint_read_len>
    bool QualityDividingReadsSetIterator<uint_read_len>::containsN() {
        return coreIterator->getRead().find('N') != string::npos;
    }

    template<typename uint_read_len>
    IndexesMapping *QualityDividingReadsSetIterator<uint_read_len>::retainVisitedIndexesMapping() {
        return new DirectMapping(allCounter);
    }

    template<typename uint_read_len>
    DividedReadsSetIterator<uint_read_len>::DividedReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator, std::istream *divSource, bool visitComplement,
            bool ignoreNReads, bool ignoreNoNReads)
            :coreIterator(coreIterator), divSource(divSource), visitComplement(visitComplement),
            ignoreNReads(ignoreNReads), ignoreNoNReads(ignoreNoNReads) {
        plainTextReadMode = confirmTextReadMode(*divSource);
        uint_reads_cnt_max readsCount = 0;
        readValue(*divSource, readsCount, plainTextReadMode);
        readValue(*divSource, currentDivIdx, plainTextReadMode);
    }

    template<typename uint_read_len>
    bool DividedReadsSetIterator<uint_read_len>::moveNext() {
        while (coreIterator->moveNext()) {
            allCounter++;
            if (visitComplement) {
                if ((allCounter != currentDivIdx) && !isIgnored()) {
                    indexesMapping.push_back(allCounter);
                    return true;
                }

                readValue(*divSource, currentDivIdx, plainTextReadMode);
            } else {
                if (allCounter == currentDivIdx) {
                    readValue(*divSource, currentDivIdx, plainTextReadMode);
                    if (!isIgnored()) {
                        indexesMapping.push_back(allCounter);
                        return true;
                    }
                }
            }
        }
        allCounter++;
        return false;
    }

    template<typename uint_read_len>
    IndexesMapping* DividedReadsSetIterator<uint_read_len>::retainVisitedIndexesMapping() {
        return new VectorMapping(std::move(indexesMapping), allCounter);
    }

    template<typename uint_read_len>
    bool DividedReadsSetIterator<uint_read_len>::isIgnored() {
        return (ignoreNReads && coreIterator->getRead().find('N') != string::npos)
            || (ignoreNoNReads && coreIterator->getRead().find('N') == string::npos);
    }

    template<typename uint_read_len>
    string& DividedReadsSetIterator<uint_read_len>::getRead() {
        return coreIterator->getRead();
    }

    template<typename uint_read_len>
    string& DividedReadsSetIterator<uint_read_len>::getQualityInfo() {
        return coreIterator->getQualityInfo();
    }

    template<typename uint_read_len>
    uint_read_len DividedReadsSetIterator<uint_read_len>::getReadLength() {
        return coreIterator->getReadLength();
    }

    template<typename uint_read_len>
    void DividedReadsSetIterator<uint_read_len>::rewind() {
        allCounter = -1;
        indexesMapping.clear();
        divSource->clear();
        divSource->seekg(0);
        confirmTextReadMode(*divSource);
        uint_reads_cnt_max readsCount = 0;
        readValue(*divSource, readsCount, plainTextReadMode);
        readValue(*divSource, currentDivIdx, plainTextReadMode);
        coreIterator->rewind();
    }

    template class QualityDividingReadsSetIterator<uint_read_len_min>;
    template class QualityDividingReadsSetIterator<uint_read_len_std>;

    template class DividedReadsSetIterator<uint_read_len_min>;
    template class DividedReadsSetIterator<uint_read_len_std>;
}
