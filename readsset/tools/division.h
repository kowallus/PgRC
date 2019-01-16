#ifndef PGTOOLS_DIVISION_H
#define PGTOOLS_DIVISION_H

#include "../iterator/ReadsSetIterator.h"

namespace PgTools {

    void divideReads(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, string outputFile,
            double error_limit, bool filterAllNReads = false);

    //DividedPackedReadsSet* generateDividedReadsSet(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, double error_limit, bool filterAllNReads = false);


}

#endif //PGTOOLS_DIVISION_H
