#include "DefaultStrandDetector.h"

#include "../pseudogenome/readslist/ReadsListTypes.h"

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
    void DefaultStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::matchReadStrands() {

    }

    template class DefaultStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultStrandDetector<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class DefaultStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultStrandDetector<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
}