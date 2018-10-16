#ifndef PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
#define PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H

#include "../DefaultPseudoGenome.h"
#include "../PackedPseudoGenome.h"
#include "../../utils/helper.h"
#include "../../readsset/persistance/ReadsSetPersistence.h"

namespace PgTools {

    class SeparatedPseudoGenomePersistence {
    private:

        static bool acceptTemporaryPseudoGenomeElement(const string &pseudoGenomePrefix, const string& fileSuffix);

    public:
        static void writePseudoGenome(PseudoGenomeBase* pgb, const string &pseudoGenomePrefix, string divisionFile = "", bool divisionComplement = false);

        static std::ifstream getPseudoGenomeSrc(const string &pseudoGenomePrefix);
        static string getPseudoGenome(const string &pseudoGenomePrefix);

        static std::ifstream getPseudoGenomeElementSrc(const string &pseudoGenomePrefix, const string& fileSuffix);
        static std::ofstream getPseudoGenomeElementDest(const string &pseudoGenomePrefix, const string &fileSuffix,
                                                        bool temporary = false);
        static void acceptTemporaryPseudoGenomeElements(const string &pseudoGenomePrefix);

        const static string PSEUDOGENOME_FILE_SUFFIX;
        const static string READSLIST_OFFSETS_FILE_SUFFIX;
        const static string READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX;
        const static string PSEUDOGENOME_PROPERTIES_SUFFIX;
        const static string READSLIST_REVERSECOMPL_FILE_SUFFIX;
        const static string READSLIST_MISMATCHESCOUNT_FILE_SUFFIX;
        const static string READSLIST_MISMATCHEDSYMBOLS_FILE_SUFFIX;
        const static string READSLIST_MISMATCHESPOS_FILE_SUFFIX;

        const static string TEMPORARY_FILE_SUFFIX;

    };

    class SeparatedReadsListWriterBase {
    public:
        virtual ~SeparatedReadsListWriterBase() {};

        virtual void writeReadsList(const string &pseudoGenomePrefix, string divisionFile = "", bool divisionComplement = false) = 0;
        virtual void writeReadsList(std::ofstream* destRlOffsets, std::ofstream* destRlIndexes, string divisionFile = "", bool divisionComplement = false) = 0;
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

        void writeReadsList(std::ofstream *destRlOffsets, std::ofstream *destRlIndexes, string divisionFile,
                            bool divisionComplement) override {
            // TODO: Write whole arrays in binary mode (for performance)
            const uint_reads_cnt_max readsCount = readsList->getReadsCount();
            const vector<uint_reads_cnt_max> orgIndexesMapping = ReadsSetPersistence::getReadsOriginalIndexes(divisionFile, divisionComplement, readsCount);
            for (uint_reads_cnt i = 0; i < readsList->getReadsCount(); i++) {
                PgSAHelpers::writeValue<uint_pg_len_max>(*destRlOffsets, readsList->getReadPosition(i));
                PgSAHelpers::writeValue<uint_reads_cnt_std>(*destRlIndexes, orgIndexesMapping[readsList->getReadOriginalIndex(i)]);
            }
        }

        void writeReadsList(const string &pseudoGenomePrefix, string divisionFile, bool divisionComplement) {
            std::ofstream destRlOffsets(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX, std::ios::out | std::ios::binary);
            std::ofstream destRlIndexes(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX, std::ios::out | std::ios::binary);
            writeReadsList(&destRlOffsets, &destRlIndexes, divisionFile, divisionComplement);
            destRlOffsets.close();
            destRlIndexes.close();
        };
    };

    class SeparatedPseudoGenomeOutputBuilder {
    private:
        const string &pseudoGenomePrefix;
        ofstream* pgDest = 0;
        ofstream* pgPropDest = 0;
        ofstream* rlOffDest = 0;
        ofstream* rlOrgIdxDest = 0;
        ofstream* rlRevCompDest = 0;
        ofstream* rlMisCntDest = 0;
        ofstream* rlMisSymDest = 0;
        ofstream* rlMisPosDest = 0;

        std::ofstream* getSingletonDest(ofstream* &dest, const string &fileSuffix);

        void freeDest(ofstream* &dest);
        void freeDests();
    public:

        SeparatedPseudoGenomeOutputBuilder(const string &pseudoGenomePrefix);

        std::ofstream* getPseudoGenomeElementDest(const string &fileSuffix);

        void build();
    };



}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
