#ifndef PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H
#define PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H

#include "../readslist/iterator/ExtendedReadsListIteratorInterface.h"
#include "SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    class SeparatedExtendedReadsListIterator : public PgTools::DefaultReadsListIteratorInterface {
    private:
        const string &pseudoGenomePrefix;
        ifstream *rlPosSrc = 0;
        ifstream *rlOrgIdxSrc = 0;
        ifstream *rlRevCompSrc = 0;
        ifstream *rlMisCntSrc = 0;
        ifstream *rlMisSymSrc = 0;
        ifstream *rlMisOffSrc = 0;

        PgTools::DefaultReadsListEntry entry;

        void initSrc(ifstream *&src, const string &fileSuffix);

        void initSrcs();

        void freeSrc(ifstream *&src);

        void freeSrcs();

    public:
        SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix);

        ~SeparatedExtendedReadsListIterator() override;

        bool moveNext() override;

        const PgTools::ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, PgSAIndex::uint_pg_len_max> &
        peekReadEntry() override;
    };
}

#endif //PGTOOLS_SEPARATEDEXTENDEDREADSLISTITERATOR_H
