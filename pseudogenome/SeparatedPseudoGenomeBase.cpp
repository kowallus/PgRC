#include "SeparatedPseudoGenomeBase.h"

namespace PgTools {

    const string SeparatedPseudoGenomeBase::PG_FILES_EXTENSION = ".pg";

    const string SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX = "_map_off" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX = "_map_len" + PG_FILES_EXTENSION;

    const string SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX = "" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX = "_prop" + PG_FILES_EXTENSION;

    const string SeparatedPseudoGenomeBase::READSLIST_POSITIONS_FILE_SUFFIX = "_rl_pos" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = "_rl_idx" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_REVERSECOMPL_FILE_SUFFIX = "_rl_rc" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX = "_rl_mis_cnt" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX =
            "_rl_mis_sym" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX =
            "_rl_mis_pos" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_OFFSETS_FILE_SUFFIX = "_rl_off" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX =
            "_rl_mis_roff" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX =
            "_rl_pr_idx" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX =
            "_rl_pr_off" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX =
            "_rl_pr_sf" + PG_FILES_EXTENSION;

    SeparatedPseudoGenomeBase::SeparatedPseudoGenomeBase(PgSAIndex::uint_pg_len_max length,
                                                         ReadsSetProperties *properties)
            : PseudoGenomeBase(length, properties) {}

    SeparatedPseudoGenomeBase::SeparatedPseudoGenomeBase(PgSAIndex::uint_pg_len_max length, istream &src)
            : PseudoGenomeBase(length, src) {}

    void SeparatedPseudoGenomeBase::getPseudoGenomeProperties(const string &pseudoGenomePrefix,
                                                              PseudoGenomeHeader *&pgh,
                                                              ReadsSetProperties *&rsProp,
                                                              bool &plainTextReadMode) {
        ifstream pgPropSrc(pseudoGenomePrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX,
                           ios_base::in | ios_base::binary);
        if (pgPropSrc.fail()) {
            fprintf(stderr, "Cannot read pseudogenome properties (%s does not open).\n",
                    pseudoGenomePrefix.c_str());
            exit(EXIT_FAILURE);
        }
        pgh = new PseudoGenomeHeader(pgPropSrc);
        rsProp = new ReadsSetProperties(pgPropSrc);
        plainTextReadMode = PgSAHelpers::confirmTextReadMode(pgPropSrc);
        pgPropSrc.close();
    }

}