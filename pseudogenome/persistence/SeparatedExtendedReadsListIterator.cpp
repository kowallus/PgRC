#include "SeparatedExtendedReadsListIterator.h"

namespace PgTools {

    SeparatedExtendedReadsListIterator::SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix) {
        initSrcs();
    }

    SeparatedExtendedReadsListIterator::~SeparatedExtendedReadsListIterator() {
        freeSrcs();
    }

    void SeparatedExtendedReadsListIterator::initSrc(ifstream *&src, const string &fileSuffix) {
        src = new ifstream(pseudoGenomePrefix + fileSuffix, ios_base::in | ios_base::binary);
        if (src->fail()) {
            delete (src);
            src = 0;
        }
    }

    void SeparatedExtendedReadsListIterator::initSrcs() {
        initSrc(rlPosSrc, SeparatedPseudoGenomePersistence::READSLIST_POSITIONS_FILE_SUFFIX);
    }

    void SeparatedExtendedReadsListIterator::freeSrc(ifstream *&src) {
        if (src) {
            src->close();
            delete (src);
            src = 0;
        }
    }

    void SeparatedExtendedReadsListIterator::freeSrcs() {

    }

    bool SeparatedExtendedReadsListIterator::moveNext() {
        return false;
    }

    const PgTools::ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &
    SeparatedExtendedReadsListIterator::peekReadEntry() {
        return entry;
    }
}