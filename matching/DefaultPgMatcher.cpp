#include "DefaultPgMatcher.h"

#include "../pseudogenome/DefaultPseudoGenome.h"
#include "../pseudogenome/PackedPseudoGenome.h"

namespace PgTools {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass>
    DefaultPgMatcher<uint_read_len, uint_reads_cnt, uint_pg_len, PseudoGenomeClass>::DefaultPgMatcher(
            PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len, PseudoGenomeClass> *pg)
            :pg(pg), readsList(pg->template getReadsList<typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>()) {}

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass>
    DefaultPgMatcher<uint_read_len, uint_reads_cnt, uint_pg_len, PseudoGenomeClass>::~DefaultPgMatcher() {

    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass>
    void DefaultPgMatcher<uint_read_len, uint_reads_cnt, uint_pg_len, PseudoGenomeClass>::exactMatchPg(string text,
                                                                                                       ofstream &offsetsDest,
                                                                                                       uint32_t minMatchLength,
                                                                                                       bool textFromSamePg) {

    }

    template class DefaultPgMatcher<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename PgSAIndex::DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>>;
    template class DefaultPgMatcher<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename PgSAIndex::DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>>;
    template class DefaultPgMatcher<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename PgSAIndex::DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>>;
    template class DefaultPgMatcher<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename PgSAIndex::DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>>;

    template class DefaultPgMatcher<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min>>;
    template class DefaultPgMatcher<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min>>;
    template class DefaultPgMatcher<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min>>;
    template class DefaultPgMatcher<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min>>;

    template class DefaultPgMatcher<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std>>;
    template class DefaultPgMatcher<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std>>;
    template class DefaultPgMatcher<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_std>>;
    template class DefaultPgMatcher<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename PgSAIndex::PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_std>>;

}

