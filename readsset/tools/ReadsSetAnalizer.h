#ifndef PGTOOLS_READSSETANALIZER_H
#define PGTOOLS_READSSETANALIZER_H

#include "../ReadsSetBase.h"

namespace PgSAReadsSet {
    class ReadsSetAnalizer {
    public:

        template<class ReadsSourceIterator>
        static ReadsSetProperties* analizeReadsSet(ReadsSourceIterator* readsIterator);
    };
}

#endif //PGTOOLS_READSSETANALIZER_H
