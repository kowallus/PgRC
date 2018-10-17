#include "ExtendedReadsListIterators.h"

PgTools::SeparatedExtendedReadsListIterator::SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix)
: pseudoGenomePrefix(pseudoGenomePrefix) {
    initSrcs();
}

PgTools::SeparatedExtendedReadsListIterator::~SeparatedExtendedReadsListIterator() {
    freeSrcs();
}

void PgTools::SeparatedExtendedReadsListIterator::initSrc(ifstream *&src) {

}

void PgTools::SeparatedExtendedReadsListIterator::initSrcs() {

}

void PgTools::SeparatedExtendedReadsListIterator::freeSrc(ifstream *&src) {
    if (src) {
        src->close();
        delete(src);
        src = 0;
    }
}

void PgTools::SeparatedExtendedReadsListIterator::freeSrcs() {

}

bool PgTools::SeparatedExtendedReadsListIterator::moveNext() {
    return false;
}

const PgTools::ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &PgTools::SeparatedExtendedReadsListIterator::peekReadEntry() {
    return entry;
}

