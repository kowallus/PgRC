#ifndef PGTOOLS_ABSTRACTSTRANDDETECTOR_H
#define PGTOOLS_ABSTRACTSTRANDDETECTOR_H

#include <unordered_map>

#include "StrandDetectorBase.h"
#include "../pseudogenome/readslist/ReadsListInterface.h"
#include "../helper.h"

using namespace PgSAIndex;

namespace PgTools {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class AbstractStrandDetector : public StrandDetectorBase {
    private:
        void init();
        void dispose();

        void quick_stats();

        std::unordered_map<uint_reads_cnt, uint_reads_cnt> readsCountPerHead;
        std::vector<uint_reads_cnt> heads;

        vector<int8_t> prepareClassification(int8_t groups_limit);

        void matchPairedReads(bool concatanated_readssrc);

    protected:
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList;

        // auxiliary structures
        vector<uint_reads_cnt> headRead;
        vector<bool> headStrand;
        vector<uint_reads_cnt> matchingContradictionsCount;

        bool pairReads(uint_reads_cnt destIdx, uint_reads_cnt srcIdx, bool sameStrand);

        uint_reads_cnt getHead(uint_reads_cnt idx);
        bool isInHeadStrand(uint_reads_cnt idx);
        bool isHeadOf(uint_reads_cnt head, uint_reads_cnt idx);

        virtual void matchReadStrands(uint_read_len overlap_threshold) = 0;

        // statistics
        uint_reads_cnt matchingContradictionTotal = 0;

    public:
        AbstractStrandDetector(
                ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList);

        virtual ~AbstractStrandDetector();

        vector<int8_t> detectStrands(uint_read_len_max overlap_threshold, bool paired_reads, bool concatenated_readssrc, int8_t groups_limit) override;
    };

}

#endif //PGTOOLS_ABSTRACTSTRANDDETECTOR_H
