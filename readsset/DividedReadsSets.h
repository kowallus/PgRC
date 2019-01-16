#ifndef PGTOOLS_DIVIDEDREADSSETS_H
#define PGTOOLS_DIVIDEDREADSSETS_H

#include "PackedConstantLengthReadsSet.h"

namespace PgSAReadsSet {

    class DividedReadsSets {
    private:
        PackedConstantLengthReadsSet* hqReadsSet = 0;
        PackedConstantLengthReadsSet* lqReadsSet = 0;
        bool nReadsLQ;
        PackedConstantLengthReadsSet* nReadsSet = 0;
        bool separateNReadsSet;


    public:
        DividedReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, double error_limit,
                         bool separateNReadsSet = false, bool nReadsLQ = false);

        PackedConstantLengthReadsSet *getHqReadsSet() const;

        PackedConstantLengthReadsSet *getLqReadsSet() const;

        bool areNReadsLQ() const;

        PackedConstantLengthReadsSet *getNReadsSet() const;

        bool isSeparateNReadsSet() const;


    };

}


#endif //PGTOOLS_DIVIDEDREADSSETS_H
