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

    public:
        static void writePseudoGenome(PseudoGenomeBase* pgb, const string &pseudoGenomePrefix, string divisionFile = "", bool divisionComplement = false);

        static std::ifstream getPseudoGenomeSrc(const string &pseudoGenomePrefix);
        static string getPseudoGenome(const string &pseudoGenomePrefix);

        static std::ifstream getPseudoGenomeElementSrc(const string &pseudoGenomePrefix, const string& fileSuffix);
        static std::ofstream getPseudoGenomeElementDest(const string &pseudoGenomePrefix, const string &fileSuffix,
                                                        bool temporary = false);
        static void acceptTemporaryPseudoGenomeElements(const string &pseudoGenomePrefix);

        const static string PSEUDOGENOME_FILE_SUFFIX;
        const static string READSLIST_POSITIONS_FILE_SUFFIX;
        const static string READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX;
        const static string PSEUDOGENOME_PROPERTIES_SUFFIX;
        const static string READSLIST_REVERSECOMPL_FILE_SUFFIX;
        const static string READSLIST_MISMATCHESCOUNT_FILE_SUFFIX;
        const static string READSLIST_MISMATCHEDSYMBOLS_FILE_SUFFIX;
        const static string READSLIST_MISMATCHESOFFSETS_FILE_SUFFIX;

        const static string TEMPORARY_FILE_SUFFIX;

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
        ofstream* rlMisOffDest = 0;

        bool disableRevComp = false;
        bool disableMismatches = false;

        void initDest(ofstream* &dest, const string &fileSuffix);
        void initReadsListDests();

        uint_reads_cnt_max readsCounter = 0;

    private:
        PseudoGenomeHeader* pgh = 0;

        void freeDest(ofstream* &dest);
        void freeDests();
    public:

        SeparatedPseudoGenomeOutputBuilder(const string &pseudoGenomePrefix, bool disableRevComp = false, bool disableMismatches = false);

        void copyPseudoGenomeHeader(const string &pseudoGenomePrefix);

        void writeReadEntry(const DefaultReadsListEntry &rlEntry);

        void writeReads(DefaultReadsListIteratorInterface* rlIt, uint_pg_len_max stopPos = (uint_pg_len_max) -1);

        void writePseudoGenome(PseudoGenomeBase* pgb, string divisionFile, bool divisionComplement);

        void build();
    };



}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
