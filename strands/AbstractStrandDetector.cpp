#include "AbstractStrandDetector.h"

#include "../pseudogenome/readslist/ReadsListTypes.h"

namespace PgTools {

    using namespace PgSAHelpers;

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::AbstractStrandDetector(
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList)
            : readsList(readsList) {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::~AbstractStrandDetector() {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::init() {
        headRead.resize(readsList->getReadsCount() + 1, 0);
        headStrand.resize(readsList->getReadsCount() + 1, true);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::dispose() {
        headRead.clear();
        headStrand.clear();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    vector<int8_t> AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::detectStrands(
            uint_read_len_max overlap_threshold, bool paired_reads, bool concatenated_readssrc, int8_t groups_limit) {

        init();

        if (paired_reads)
            matchPairedReads(concatenated_readssrc);

        matchReadStrands();

        vector<int8_t> res = prepareClassification(groups_limit);
        quick_stats();

        dispose();

        return res;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool
    AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::pairReads(uint_reads_cnt destIdx,
                                                                                                  uint_reads_cnt srcIdx,
                                                                                                  bool sameStrand) {
        uint_reads_cnt headIdx = headRead[destIdx]?getHead(destIdx):destIdx;
        sameStrand = ((headRead[destIdx]?headStrand[destIdx]:true) == sameStrand);

        uint_reads_cnt joiningIdx = headRead[srcIdx]?getHead(srcIdx):srcIdx;
        sameStrand = ((headRead[srcIdx]?headStrand[srcIdx]:true) == sameStrand);

        if (headIdx == joiningIdx) {
            if (!sameStrand) {
                matchingContradictionCount++;
                return false;
            }
            return true;
        }

        headRead[joiningIdx] = headIdx;
        headStrand[joiningIdx] = sameStrand;

        return true;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_reads_cnt
    AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getHead(uint_reads_cnt idx) {
        uint_reads_cnt headIdx = headRead[idx];
        if (headRead[idx] == 0)
            return idx;

        bool sameStrand = headStrand[headIdx];
        headRead[idx] = getHead(headIdx);
        headStrand[idx] = (headStrand[idx] == (headStrand[headIdx] == sameStrand));

        return headRead[idx];
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::isInHeadStrand(
            uint_reads_cnt idx) {
        getHead(idx);
        return headStrand[idx];
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool
    AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::isHeadOf(uint_reads_cnt head,
                                                                                                 uint_reads_cnt idx) {
        return getHead(idx) == head;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::matchPairedReads(
            bool concatanated_readssrc) {
        uint_reads_cnt readsCount = readsList->getReadsCount();
        if (readsCount % 2 == 1) {
            cerr << "Unpaired input. Expected even reads count (not " << readsCount << ")";
            exit(EXIT_FAILURE);
        }
        if (concatanated_readssrc) {
            for(uint_reads_cnt i = 1; i <= readsCount / 2; i++)
                pairReads(i, i + readsCount/2, false);
        } else {
            for (uint_reads_cnt i = 1; i <= readsCount; i += 2)
                pairReads(i, i + 1, false);
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    vector<int8_t>
    AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::prepareClassification(
            int8_t groups_limit) {
        clock_checkpoint();
        uint_reads_cnt readsCount = readsList->getReadsCount();
        vector<int8_t> result(readsCount);

        for(int i = 1; i <= readsCount; i++) {
            uint_reads_cnt headIdx = getHead(i);
            readsCountPerHead[headIdx]++;
        }

        heads.clear();
        heads.reserve(readsCountPerHead.size());
        for(auto kv : readsCountPerHead)
            heads.push_back(kv.first);

        std::sort(heads.begin(), heads.end(), [this](const uint_reads_cnt& headIdx1, const uint_reads_cnt& headIdx2) -> bool
            { return readsCountPerHead[headIdx1] > readsCountPerHead[headIdx2]; });

        cout << "Reads/strands classification assembled in " << clock_millis() << " msec\n\n";
        return result;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::quick_stats() {
        cout << "Matched read groups count: " << readsCountPerHead.size() << endl;
        const uint_reads_cnt GROUPS_TO_SHOW = 5;
        cout << "Largest (" << GROUPS_TO_SHOW << ") matched read groups count: " << readsCountPerHead[heads[0]];
        uint_reads_cnt i = 0;
        while (++i < readsCountPerHead.size() && i < GROUPS_TO_SHOW) {
            cout << ", " << readsCountPerHead[heads[i]];
        };
        cout << endl;
        cout << "Contradicted read strands during matching: " << matchingContradictionCount << endl;
    }

    template class AbstractStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class AbstractStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class AbstractStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class AbstractStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
}
