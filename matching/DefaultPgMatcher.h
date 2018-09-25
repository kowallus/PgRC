#ifndef PGTOOLS_DEFAULTPGMATCHER_H
#define PGTOOLS_DEFAULTPGMATCHER_H

#include "../pseudogenome/PseudoGenomeInterface.h"
#include "../pseudogenome/readslist/ReadsListTypes.h"
#include "PgMatcherBase.h"

namespace PgTools {

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass>
    class DefaultPgMatcher: public PgMatcherBase {
    private:
        PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len, PseudoGenomeClass>* pg;
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type> *readsList;
    public:
        DefaultPgMatcher(
                PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len, PseudoGenomeClass> *pg);

        void exactMatchPg(string text, ofstream &offsetsDest, uint32_t minMatchLength, bool textFromSamePg) override;

        virtual ~DefaultPgMatcher();
    };

}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
