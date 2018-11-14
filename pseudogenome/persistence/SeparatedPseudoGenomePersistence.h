#ifndef PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
#define PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H

#include "../DefaultPseudoGenome.h"
#include "../PackedPseudoGenome.h"
#include "../../utils/helper.h"
#include "../../readsset/persistance/ReadsSetPersistence.h"
#include "../readslist/iterator/ReadsListIteratorExtendedWrapper.h"

namespace PgTools {

    class SeparatedPseudoGenomePersistence {
    private:

        static bool acceptTemporaryPseudoGenomeElement(const string &pseudoGenomePrefix, const string& fileSuffix);

        static void appendIndexesFromPg(string pgFilePrefix, vector<uint_reads_cnt_std> &idxs);
        static void writePairMapping(basic_string<char> &pgFilePrefix, vector<uint_reads_cnt_std> orgIdxs);

    public:
        static void writePseudoGenome(PseudoGenomeBase* pgb, const string &pseudoGenomePrefix,
                string divisionFile = "", bool divisionComplement = false, bool revComplPairFile = false);

        static std::ifstream getPseudoGenomeSrc(const string &pseudoGenomePrefix);
        static string getPseudoGenome(const string &pseudoGenomePrefix);

        static void getPseudoGenomePropertes(const string &pseudoGenomePrefix, PseudoGenomeHeader* &pgh, bool& plainTextReadMode);

        static std::ifstream getPseudoGenomeElementSrc(const string &pseudoGenomePrefix, const string& fileSuffix);
        static std::ofstream getPseudoGenomeElementDest(const string &pseudoGenomePrefix, const string &fileSuffix,
                                                        bool temporary = false);
        static void acceptTemporaryPseudoGenomeElements(const string &pseudoGenomePrefix);

        const static string PSEUDOGENOME_FILE_SUFFIX;
        const static string READSLIST_POSITIONS_FILE_SUFFIX;
        const static string READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX;
        const static string PSEUDOGENOME_PROPERTIES_SUFFIX;
        const static string READSLIST_REVERSECOMPL_FILE_SUFFIX;
        const static string READSLIST_MISMATCHES_COUNT_FILE_SUFFIX;
        const static string READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX;
        const static string READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX;

        const static string READSLIST_OFFSETS_FILE_SUFFIX;
        const static string READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX;

        const static string READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX;
        const static string READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX;

        const static string TEMPORARY_FILE_SUFFIX;

        static bool enableReadPositionRepresentation;
        static bool enableRevOffsetMismatchesRepresentation;

        static void dumpPgPairs(vector<string> pgFilePrefixes);

    };

    class SeparatedPseudoGenomeOutputBuilder {
    private:
        const string &pseudoGenomePrefix;
        ofstream* pgDest = 0;
        ofstream* pgPropDest = 0;
        ofstream* rlPosDest = 0;
        ofstream* rlOrgIdxDest = 0;
        ofstream* rlRevCompDest = 0;
        ofstream* rlMisCntDest = 0;
        ofstream* rlMisSymDest = 0;
        ofstream* rlMisPosDest = 0;

        ofstream* rlOffDest = 0;
        ofstream* rlMisRevOffDest = 0;

        bool disableRevComp = false;
        bool disableMismatches = false;

        void initDest(ofstream* &dest, const string &fileSuffix);
        void initReadsListDests();

        DefaultReadsListIteratorInterface *rlIt = 0;
        bool iterationPaused = false;

        uint_pg_len_max lastWrittenPos = 0;
        void writeReadEntry(const DefaultReadsListEntry &rlEntry);
        uint_reads_cnt_max readsCounter = 0;

        PseudoGenomeHeader* pgh = 0;

        void freeDest(ofstream* &dest);
        void freeDests();
    public:

        SeparatedPseudoGenomeOutputBuilder(const string &pseudoGenomePrefix, bool disableRevComp = false, bool disableMismatches = false);

        void copyPseudoGenomeHeader(const string &pseudoGenomePrefix);

        void setReadsSourceIterator(DefaultReadsListIteratorInterface *rlIt);
        uint_pg_len_max writeReadsFromIterator(uint_pg_len_max stopPos = (uint_pg_len_max) -1);
        void writeExtraReadEntry(const DefaultReadsListEntry &rlEntry);

        void writePseudoGenome(PseudoGenomeBase* pgb, string divisionFile, bool divisionComplement, bool revComplPairFile = "false");

        void build();
    };



}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
