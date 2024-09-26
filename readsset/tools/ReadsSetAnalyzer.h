#ifndef PGTOOLS_READSSETANALIZER_H
#define PGTOOLS_READSSETANALIZER_H

#include "../ReadsSetBase.h"

namespace PgReadsSet {
    class ReadsSetAnalyzer {
    public:

        template<class ReadsSourceIterator>
        static ReadsSetProperties* analyzeReadsSet(ReadsSourceIterator *readsIterator);
    };
}

#endif //PGTOOLS_READSSETANALIZER_H
