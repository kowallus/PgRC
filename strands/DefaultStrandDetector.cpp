#include "DefaultStrandDetector.h"

#include "../pseudogenome/readslist/ReadsListTypes.h"

//#define DETECTOR_VALIDATION

namespace PgTools {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::DefaultStrandDetector(
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList)
            : AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>(readsList) {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::~DefaultStrandDetector() {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_pg_len DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getReadPosition(
            uint_reads_cnt rlIdx) {
        return this->readsList->getReadPosition(rlIdx);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len
    DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getFollowingReadOverlap(
            uint_reads_cnt rlIdx) {
        return this->readsList->getReadLength(rlIdx) - (getReadPosition(rlIdx + 1) - getReadPosition(rlIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_reads_cnt DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getOriginalIndex(
            uint_reads_cnt rlIdx) {
        return this->readsList->getReadOriginalIndex(rlIdx);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::matchReadStrands(uint_read_len overlap_threshold) {
#ifdef DETECTOR_VALIDATION
        vector<uint_reads_cnt> destIdxs(readsCount + 1);
        vector<uint_reads_cnt> srcIdxs(readsCount + 1);
#endif
        uint_reads_cnt readsCount = this->readsList->getReadsCount();
        vector<uint_reads_cnt> rlIdxs(readsCount);
        for(uint_reads_cnt i = 0; i < readsCount - 1; i++)
            rlIdxs[i] = i;

        std::sort(rlIdxs.begin(), rlIdxs.end(), [this](const uint_reads_cnt& rlIdx1, const uint_reads_cnt& rlIdx2) -> bool
        { return getFollowingReadOverlap(rlIdx1) > getFollowingReadOverlap(rlIdx2); });

        for(uint_reads_cnt i = 0; i < readsCount - 1; i++) {
            uint_reads_cnt rlIdx = rlIdxs[i];
            if (getFollowingReadOverlap(rlIdx) < overlap_threshold) {
                cout << "Number of reads in readsList with at least " << (int) overlap_threshold << " symbols overlap: " << i << endl;
                break;
            }

            this->pairReads(getOriginalIndex(rlIdx) + 1, getOriginalIndex(rlIdx + 1) + 1, true);
#ifdef DETECTOR_VALIDATION
            destIdxs[getOriginalIndex(rlIdx + 1) + 1] = getOriginalIndex(rlIdx) + 1;
            srcIdxs[getOriginalIndex(rlIdx) + 1] = getOriginalIndex(rlIdx + 1) + 1;
#endif
        }
    }

    template class DefaultStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class DefaultStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
}