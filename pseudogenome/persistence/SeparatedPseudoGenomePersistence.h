#ifndef PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
#define PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H

#include <memory>
#include "../DefaultPseudoGenome.h"
#include "../PackedPseudoGenome.h"
#include "../SeparatedPseudoGenome.h"
#include "../../utils/helper.h"
#include "../../readsset/persistance/ReadsSetPersistence.h"
#include "../readslist/iterator/ReadsListIteratorExtendedWrapper.h"
#include "../readslist/SeparatedExtendedReadsList.h"
#include "../../coders/CodersLib.h"

namespace PgTools {

    class SeparatedPseudoGenomePersistence {
    private:

        static bool acceptTemporaryPseudoGenomeElement(const string &pseudoGenomePrefix, const string& fileSuffix, bool alwaysRemoveExisting);

        static void appendIndexesFromPg(string pgFilePrefix, vector<uint_reads_cnt_std> &idxs);
        static void writePairMapping(basic_string<char> &pgFilePrefix, vector<uint_reads_cnt_std> orgIdxs);

    public:
        static void writePseudoGenome(PseudoGenomeBase* pgb, const string &pseudoGenomePrefix,
                IndexesMapping* orgIndexesMapping, bool revComplPairFile = false);
        static void writeSeparatedPseudoGenome(SeparatedPseudoGenome *sPg, const string &pseudoGenomePrefix,
            bool skipPgSequence = false);
        static void compressSeparatedPseudoGenomeReadsList(SeparatedPseudoGenome *sPg, ostream* pgrcOut,
                uint8_t coder_level, bool ignoreOffDest, bool skipPgSequence = true);
        static SeparatedPseudoGenome* loadSeparatedPseudoGenome(const string &pgPrefix, bool skipReadsList = false);

        static std::ifstream getPseudoGenomeSrc(const string &pseudoGenomePrefix);
        static string loadPseudoGenomeSequence(const string &pseudoGenomePrefix);

        static std::ifstream getPseudoGenomeElementSrc(const string &pseudoGenomePrefix, const string& fileSuffix);
        static std::ofstream getPseudoGenomeElementDest(const string &pseudoGenomePrefix, const string &fileSuffix,
                                                        bool temporary = false);
        static void acceptTemporaryPseudoGenomeElements(const string &pseudoGenomePrefix, bool clearReadsListDescriptionIfOverride);

        const static string TEMPORARY_FILE_SUFFIX;

        static bool enableReadPositionRepresentation;
        static bool enableRevOffsetMismatchesRepresentation;

        static void dumpPgPairs(vector<string> pgFilePrefixes);
        static void compressReadsOrder(ostream &pgrcOut, const vector<uint_reads_cnt_std>& orgIdxs, uint8_t coder_level,
                bool completeOrderInfo = false, bool ignorePairOrderInformation = false, bool singleFileMode = true);
        static void decompressReadsOrder(istream &pgrcIn, vector<uint_reads_cnt_std>& rlIdxOrder,
                                       bool completeOrderInfo = false, bool ignorePairOrderInformation = false, bool singleFileMode = true);

        static void writePseudoGenomeSequence(string &pgSequence, string pgPrefix);

        template <typename uint_pg_len>
        static void compressReadsPgPositions(ostream &pgrcOut, vector<uint_pg_len_max> orgIdx2PgPos,
                uint_pg_len_max joinedPgLength, uint8_t coder_level, bool singleFileMode, bool deltaPairEncodingEnabled = true);
        template <typename uint_pg_len>
        static void decompressReadsPgPositions(istream &pgrcIn, vector<uint_pg_len> &pgPos, PgRCParams* params);
    };

    class SeparatedPseudoGenomeOutputBuilder {
    private:
        const string pseudoGenomePrefix;

        bool onTheFlyMode() { return !pseudoGenomePrefix.empty(); }

        ostream* pgDest = nullptr;
        ostream* pgPropDest = nullptr;
        ostream* rlPosDest = nullptr;
        ostream* rlOrgIdxDest = nullptr;
        ostream* rlRevCompDest = nullptr;
        ostream* rlMisCntDest = nullptr;
        ostream* rlMisSymDest = nullptr;
        ostream* rlMisPosDest = nullptr;

        ostream* rlOffDest = nullptr;
        ostream* rlMisRevOffDest = nullptr;

        bool disableRevComp = false;
        bool disableMismatches = false;

        void initDest(ostream *&dest, const string &fileSuffix);
        void initReadsListDests();

        DefaultReadsListIteratorInterface *rlIt = nullptr;
        bool iterationPaused = false;

        uint_pg_len_max lastWrittenPos = 0;
        uint_reads_cnt_max readsCounter = 0;

        PseudoGenomeHeader* pgh = nullptr;
        ReadsSetProperties* rsProp = nullptr;

        void prebuildAssert(bool requireOnTheFlyMode);
        void buildProps();

        void freeDest(ostream* &dest);
        void freeDests();

        string toStringAndSeparateZeros(ostream* dest, string& zeroFlags);
        string toString(ostream* dest);
        void compressRlMisRevOffDest(ostream &pgrcOut, uint8_t coder_level, bool separateFirstOffsetMode = false, bool transposeMode = false);
        void destToFile(ostream *dest, const string &fileName);
    public:

        SeparatedPseudoGenomeOutputBuilder(bool disableRevComp = false, bool disableMismatches = false);
        SeparatedPseudoGenomeOutputBuilder(const string pseudoGenomePrefix, bool disableRevComp = false, bool disableMismatches = false);

        virtual ~SeparatedPseudoGenomeOutputBuilder();

        void copyPseudoGenomeProperties(const string &pseudoGenomePrefix);
        void copyPseudoGenomeProperties(SeparatedPseudoGenome* sPg);

        void writeReadEntry(const DefaultReadsListEntry &rlEntry);

        void setReadsSourceIterator(DefaultReadsListIteratorInterface *rlIt);
        uint_pg_len_max writeReadsFromIterator(uint_pg_len_max stopPos = (uint_pg_len_max) -1);
        void writeExtraReadEntry(const DefaultReadsListEntry &rlEntry);

        void writePseudoGenome(PseudoGenomeBase* pgb, IndexesMapping* orgIndexesMapping, bool revComplPairFile = "false");
        void feedSeparatedPseudoGenome(SeparatedPseudoGenome *sPg, bool skipPgSequence = false);
        void appendPseudoGenome(const string &pg);

        void build();
        void build(const string &pgPrefix);

        void compressedBuild(ostream &pgrcOut, uint8_t coder_level, bool ignoreOffDest = false);

        void updateOriginalIndexesIn(SeparatedPseudoGenome *sPg);

        string reorderingSymbolsExclusiveMismatchEncoding(string &rlMisSym);
    };



}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOMEPERSISTENCE_H
