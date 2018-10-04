#ifndef PGTOOLS_READSSETPERSISTENCE_H
#define PGTOOLS_READSSETPERSISTENCE_H

#include "../DefaultReadsSet.h"
#include "../PackedReadsSet.h"
#include "../iterator/ReadsSetIterator.h"

namespace PgSAReadsSet {

    class ReadsSetPersistence {

    public:
        static ReadsSourceIteratorTemplate<uint_read_len_max>* createReadsIterator(string filename, string pairfile = "");
        static FASTQReadsSourceIterator<uint_read_len_max>* createFastQReadsIterator(string filename, string pairfile = "");

    };

}


#endif //PGTOOLS_READSSETPERSISTENCE_H
