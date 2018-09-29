#ifndef PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
#define PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H

#include "../DefaultPseudoGenome.h"
#include "../PackedPseudoGenome.h"
#include "../../helper.h"

namespace PgTools {

    class SeparatedPseudoGenomePersistence {
    public:
        static void writePseudoGenome(PseudoGenomeBase* pgb, string pseudoGenomePrefix);

        const static string PSEUDOGENOME_FILE_SUFFIX;
        const static string READSLIST_OFFSETS_FILE_SUFFIX;
        const static string READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX;
        const static string PSEUDOGENOME_PROPERTIES_SUFFIX;

    };

    class SeparatedReadsListWriterBase {
    public:
        virtual ~SeparatedReadsListWriterBase() {};

        virtual void writeReadsList(string pseudoGenomePrefix) = 0;
    };

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class SeparatedReadsListWriter: public SeparatedReadsListWriterBase {
    private:
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList;
    public:
        SeparatedReadsListWriter(
                ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList) : readsList(
                readsList) {}

    public:
        void writeReadsList(string pseudoGenomePrefix) {
            std::ofstream destRlOffsets(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX, std::ios::out | std::ios::binary);
            std::ofstream destRlIndexes(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX, std::ios::out | std::ios::binary);
            for (uint_reads_cnt i = 0; i < readsList->getReadsCount(); i++) {
                PgSAHelpers::writeValue<uint_pg_len_max>(destRlOffsets, readsList->getReadPosition(i));
                PgSAHelpers::writeValue<uint_reads_cnt_std>(destRlIndexes, readsList->getReadOriginalIndex(i));
            }
            destRlOffsets.close();
            destRlIndexes.close();
        };
    };

}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
