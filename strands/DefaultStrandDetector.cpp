#include "DefaultStrandDetector.h"

#include "../pseudogenome/readslist/ReadsListTypes.h"
#include "../pseudogenome/readslist/ReadsListTypes.h"

namespace PgTools {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::DefaultStrandDetector(
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList)
            : readsList(readsList){
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::~DefaultStrandDetector() {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    vector<int8_t> DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::detectStrands(
            int8_t groups_limit, bool paired_reads, bool concatanated_readssrc) {
        return vector<int8_t>();
    }

    template class DefaultStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class DefaultStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
}