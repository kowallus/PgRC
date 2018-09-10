#ifndef PGTOOLS_DEFAULTSTRANDDETECTOR_H
#define PGTOOLS_DEFAULTSTRANDDETECTOR_H

#include "StrandDetectorBase.h"
#include "../pseudogenome/readslist/ReadsListInterface.h"

using namespace PgSAIndex;

namespace PgTools {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class DefaultStrandDetector: public StrandDetectorBase {
    private:
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList = 0;
    public:
        DefaultStrandDetector(ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList);
        ~DefaultStrandDetector();

        vector<int8_t> detectStrands(int8_t groups_limit, bool paired_reads, bool concatanated_readssrc) override;
    };

}

#endif //PGTOOLS_DEFAULTSTRANDDETECTOR_H
