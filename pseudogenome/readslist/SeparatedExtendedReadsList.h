#ifndef PGTOOLS_SEPARATEDEXTENDEDREADSLIST_H
#define PGTOOLS_SEPARATEDEXTENDEDREADSLIST_H

#include "iterator/ExtendedReadsListIteratorInterface.h"
#include "../SeparatedPseudoGenomeBase.h"
#include "../../coders/SymbolsPackingFacility.h"
#include "../DefaultPseudoGenome.h"
#include "../../pgrc/pgrc-params.h"

namespace PgTools {

    template<int maxMismatches>
    class SeparatedExtendedReadsListIterator;

    typedef SeparatedExtendedReadsListIterator<UINT8_MAX> DefaultSeparatedExtendedReadsListIterator;
    typedef SeparatedExtendedReadsListIterator<0> SimpleSeparatedReadsListIterator;

    class ExtendedReadsListWithConstantAccessOption : public DefaultReadsListIteratorInterface {

        DefaultReadsListEntry entry;
        int64_t current = -1;
        uint64_t curMisCumCount = 0;

    public:

        uint_read_len_max readLength;
        uint_reads_cnt_std readsCount;

        vector<uint_read_len_min> off;
        vector<uint_reads_cnt_std> orgIdx;
        vector<uint8_t> revComp;
        vector<uint_read_len_min> misCnt;
        vector<uint8_t> misSymCode;
        vector<uint_read_len_min> misOff;

        // constant access features
        vector<uint_pg_len_max> pos;
        vector<uint_reads_cnt_max> misCumCount;

        ExtendedReadsListWithConstantAccessOption(uint_read_len_max readLength) : readLength(readLength) {}

        virtual ~ExtendedReadsListWithConstantAccessOption() {};

        void enableConstantAccess(bool disableIterationMode = false, bool skipPositions = false);

        bool isConstantAccessEnalbed();

        inline uint8_t getMisCount(uint_reads_cnt_max rlIdx) {
            return misCumCount[rlIdx + 1] - misCumCount[rlIdx];
        }

        inline uint8_t getMisSymCode(uint_reads_cnt_max rlIdx, uint8_t misIdx) {
            return misSymCode[misCumCount[rlIdx] + misIdx];
        }

        inline uint8_t getMisOff(uint_reads_cnt_max rlIdx, uint8_t misIdx) {
            return misOff[misCumCount[rlIdx] + misIdx];
        }

        void copyMismatchesToEntry(uint_reads_cnt_max idx, DefaultReadsListEntry &entry) {
            if (misCumCount.size() == 0)
                return;
            uint8_t mismatchesCount = getMisCount(idx);
            for (uint8_t i = 0; i < mismatchesCount; i++)
                entry.addMismatch(getMisSymCode(idx, i), getMisOff(idx, i));
        }

        // iterator routines

        bool moveNext() override;
        void rewind() override;
        DefaultReadsListEntry &peekReadEntry() { return entry; };

        // reads list factories

        static ExtendedReadsListWithConstantAccessOption *loadConstantAccessExtendedReadsList(const string &pseudoGenomePrefix,
                                                                                    uint_pg_len_max pgLengthPosGuard = 0,
                                                                                    bool skipMismatches = false);

        static ExtendedReadsListWithConstantAccessOption *loadConstantAccessExtendedReadsList(istream &pgrcIn,
                PseudoGenomeHeader *pgh, ReadsSetProperties *rsProp, const string validationPgPrefix,
                PgRCParams* params, bool disableRevCompl = false, bool disableMismatches = false,
                bool separateFirstOffsetMode = false);

        static ExtendedReadsListWithConstantAccessOption *loadConstantAccessExtendedReadsList(
                DefaultSeparatedExtendedReadsListIterator &rl,
                uint_pg_len_max pgLengthPosGuard = 0, bool skipMismatches = false);

        bool isRevCompEnabled();

        bool areMismatchesEnabled();

        bool getRevComp(uint_reads_cnt_std idx);
    };

    template<int maxMismatches>
    class SeparatedExtendedReadsListIterator
            : public ExtendedReadsListIteratorInterface<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> {
    private:
        const string pseudoGenomePrefix;

        bool fromFileMode() { return !pseudoGenomePrefix.empty(); }

        istream *rlPosSrc = 0;
        istream *rlOrgIdxSrc = 0;
        istream *rlRevCompSrc = 0;
        istream *rlMisCntSrc = 0;
        istream *rlMisSymSrc = 0;
        istream *rlMisOffSrc = 0;

        istream *rlOffSrc = 0;
        istream *rlMisRevOffSrc = 0;

        bool ownProps = true;
        PseudoGenomeHeader *pgh = 0;
        ReadsSetProperties *rsProp = 0;
        bool plainTextReadMode = false;

        int64_t current = -1;

        ReadsListEntry<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> entry;

        void initSrc(istream *&src, const string &fileSuffix);

        void initSrcs();

        void freeSrc(istream *&src);

        void freeSrcs();

        SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix);


    public:
        ~SeparatedExtendedReadsListIterator() override;

        bool moveNext() override;

        void rewind() override;

        ReadsListEntry<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &
        peekReadEntry() override;

        bool isRevCompEnabled();

        bool areMismatchesEnabled();

        static SeparatedExtendedReadsListIterator<maxMismatches> *getIterator(const string &pseudoGenomePrefix);


        friend ExtendedReadsListWithConstantAccessOption *ExtendedReadsListWithConstantAccessOption::
        loadConstantAccessExtendedReadsList(const string &pseudoGenomePrefix,
                                            uint_pg_len_max pgLengthPosGuard, bool skipMismatches);

        friend ExtendedReadsListWithConstantAccessOption *ExtendedReadsListWithConstantAccessOption::
        loadConstantAccessExtendedReadsList(DefaultSeparatedExtendedReadsListIterator &rl,
                                            uint_pg_len_max pgLengthPosGuard, bool skipMismatches);

        friend ExtendedReadsListWithConstantAccessOption *ExtendedReadsListWithConstantAccessOption::
        loadConstantAccessExtendedReadsList(istream &pgrcIn,
                                            PseudoGenomeHeader *pgh, ReadsSetProperties *rsProp,
                                            const string validationPgPrefix, PgRCParams* params, bool disableRevCompl,
                                            bool disableMismatches, bool separateFirstOffsetMode);
    };

}

#endif //PGTOOLS_SEPARATEDEXTENDEDREADSLIST_H
