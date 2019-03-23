#ifndef PGTOOLS_SEPARATEDEXTENDEDPSEUDOGENOMEBASE_H
#define PGTOOLS_SEPARATEDEXTENDEDPSEUDOGENOMEBASE_H

#include "PseudoGenomeBase.h"

namespace PgTools {

    using namespace PgSAIndex;

    class SeparatedPseudoGenomeBase : public PseudoGenomeBase {
    public:

        SeparatedPseudoGenomeBase(uint_pg_len_max length, ReadsSetProperties *properties);

        SeparatedPseudoGenomeBase(uint_pg_len_max length, istream &src);

        virtual ~SeparatedPseudoGenomeBase() {};



        const static string PG_FILES_EXTENSION;

        const static string PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX;
        const static string PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX;

        const static string PSEUDOGENOME_FILE_SUFFIX;
        const static string PSEUDOGENOME_PROPERTIES_SUFFIX;

        const static string READSLIST_POSITIONS_FILE_SUFFIX;
        const static string READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX;
        const static string READSLIST_REVERSECOMPL_FILE_SUFFIX;
        const static string READSLIST_MISMATCHES_COUNT_FILE_SUFFIX;
        const static string READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX;
        const static string READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX;
        const static string READSLIST_OFFSETS_FILE_SUFFIX;
        const static string READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX;
        const static string READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX;
        const static string READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX;
        const static string READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX;

        static void getPseudoGenomeProperties(const string &pseudoGenomePrefix, PseudoGenomeHeader *&pgh,
                                              ReadsSetProperties *&rsProp, bool &plainTextReadMode);

    };
}

#include "PseudoGenomeBase.h"

#endif //PGTOOLS_SEPARATEDEXTENDEDPSEUDOGENOMEBASE_H
