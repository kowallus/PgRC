#ifndef PGTOOLS_PGRC_DATA_H
#define PGTOOLS_PGRC_DATA_H

#include "pg-config.h"
#include "../readsset/DividedPCLReadsSets.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    class PgRCData {
    public:

        DividedPCLReadsSets *divReadsSets = 0;
        SeparatedPseudoGenome *hqPg = 0;
        SeparatedPseudoGenome *lqPg = 0;
        SeparatedPseudoGenome *nPg = 0;

        vector<uint_reads_cnt_std> rlIdxOrder;
        vector<uint_pg_len_max> orgIdx2PgPos;
        vector<uint_pg_len_std> orgIdx2StdPgPos;

        void disposeChainData() {
            if (divReadsSets) {
                delete (divReadsSets);
                divReadsSets = 0;
            }
            if (hqPg) {
                delete (hqPg);
                hqPg = 0;
            }
            if (lqPg) {
                delete (lqPg);
                lqPg = 0;
            }
            if (nPg) {
                delete (nPg);
                nPg = 0;
            }
        }

    };

}

#endif //PGTOOLS_PGRC_DATA_H
