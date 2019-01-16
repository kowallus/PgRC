#include "DividedReadsSets.h"

namespace PgSAReadsSet {
    DividedReadsSets::DividedReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator,
                                       double error_limit, bool separateNReadsSet, bool nReadsLQ) :
                                       nReadsLQ(nReadsLQ), separateNReadsSet(separateNReadsSet) {

    }

    PackedConstantLengthReadsSet *DividedReadsSets::getHqReadsSet() const {
        return hqReadsSet;
    }

    PackedConstantLengthReadsSet *DividedReadsSets::getLqReadsSet() const {
        return lqReadsSet;
    }

    bool DividedReadsSets::areNReadsLQ() const {
        return nReadsLQ;
    }

    PackedConstantLengthReadsSet *DividedReadsSets::getNReadsSet() const {
        return nReadsSet;
    }

    bool DividedReadsSets::isSeparateNReadsSet() const {
        return separateNReadsSet;
    }
}
