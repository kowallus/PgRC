#ifndef PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H
#define PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H

#include "iterator/ExtendedReadsListIteratorInterface.h"
#include "../persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    class ConstantAccessExtendedReadsList {
    public:

        uint_read_len_max readLength;

        vector<uint_pg_len_max> pos;
        vector<uint_reads_cnt_std> orgIdx;
        vector<bool> revComp;
        vector<uint_reads_cnt_max> misCumCount;
        vector<uint8_t> misSymCode;
        vector<uint8_t> misOff;
        bool misRevOffMode = false;

        ConstantAccessExtendedReadsList(uint_read_len_max readLength) : readLength(readLength) {}

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
            for(uint8_t i = 0; i < mismatchesCount; i++)
                entry.addMismatch(getMisSymCode(idx, i), getMisOff(idx, i));
            if (misRevOffMode)
                convertMisOffsets2RevOffsets(entry.mismatchOffset, entry.mismatchesCount, readLength);
        }
    };

    template <int maxMismatches>
    class SeparatedExtendedReadsList : public ExtendedReadsListIteratorInterface<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> {
    private:
        const string &pseudoGenomePrefix;
        ifstream *rlPosSrc = 0;
        ifstream *rlOrgIdxSrc = 0;
        ifstream *rlRevCompSrc = 0;
        ifstream *rlMisCntSrc = 0;
        ifstream *rlMisSymSrc = 0;
        ifstream *rlMisOffSrc = 0;

        ifstream* rlOffSrc = 0;
        ifstream* rlMisRevOffSrc = 0;

        PseudoGenomeHeader* pgh = 0;
        bool plainTextReadMode = false;

        int64_t current = -1;

        ReadsListEntry<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> entry;

        void initSrc(ifstream *&src, const string &fileSuffix);
        void initSrcs();
        void freeSrc(ifstream *&src);
        void freeSrcs();

        SeparatedExtendedReadsList(const string &pseudoGenomePrefix);

    public:
        ~SeparatedExtendedReadsList() override;

        bool moveNext() override;

        ReadsListEntry<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max>&
        peekReadEntry() override;

        bool isRevCompEnabled();

        bool areMismatchesEnabled();

        static SeparatedExtendedReadsList<maxMismatches>* getIterator(const string &pseudoGenomePrefix);
        static ConstantAccessExtendedReadsList* loadConstantAccessExtendedReadsList(const string &pseudoGenomePrefix,
                                                                                    uint_pg_len_max pgLengthPosGuard = 0);
    };

    typedef SeparatedExtendedReadsList<UINT8_MAX> DefaultSeparatedExtendedReadsList;
    typedef SeparatedExtendedReadsList<UINT8_MAX> DefaultSeparatedExtendedReadsListIterator;
    typedef SeparatedExtendedReadsList<0> SimpleSeparatedReadsListIterator;

}

#endif //PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H
