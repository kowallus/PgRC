#ifndef PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H
#define PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H

#include "../readslist/iterator/ExtendedReadsListIteratorInterface.h"
#include "SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    template <int maxMismatches>
    class SeparatedExtendedReadsListIterator : public ExtendedReadsListIteratorInterface<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> {
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

    public:
        SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix);

        ~SeparatedExtendedReadsListIterator() override;

        bool moveNext() override;

        ReadsListEntry<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max>&
        peekReadEntry() override;

        bool isRevCompEnabled();

        bool areMismatchesEnabled();
    };

    typedef SeparatedExtendedReadsListIterator<UINT8_MAX> DefaultSeparatedExtendedReadsListIterator;
    typedef SeparatedExtendedReadsListIterator<0> SimpleSeparatedReadsListIterator;

}

#endif //PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H
