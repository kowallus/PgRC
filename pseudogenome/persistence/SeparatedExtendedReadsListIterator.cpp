#include "SeparatedExtendedReadsListIterator.h"

namespace PgTools {

    SeparatedExtendedReadsListIterator::SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix) {
        ifstream pgPropSrc(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX,
                           ios_base::in | ios_base::binary);
        if (pgPropSrc.fail()) {
            fprintf(stderr, "Cannot read pseudogenome properties (%s does not open).\n",
                    pseudoGenomePrefix.c_str());
            exit(EXIT_FAILURE);
        }
        pgh = new PseudoGenomeHeader(pgPropSrc);
        pgPropSrc.close();
        initSrcs();
    }

    SeparatedExtendedReadsListIterator::~SeparatedExtendedReadsListIterator() {
        delete(pgh);
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
        initSrc(rlOrgIdxSrc, SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        initSrc(rlRevCompSrc, SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        initSrc(rlMisCntSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHESCOUNT_FILE_SUFFIX);
        initSrc(rlMisSymSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHEDSYMBOLS_FILE_SUFFIX);
        initSrc(rlMisOffSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHESOFFSETS_FILE_SUFFIX);
    }

    void SeparatedExtendedReadsListIterator::freeSrc(ifstream *&src) {
        if (src) {
            src->close();
            delete (src);
            src = 0;
        }
    }

    void SeparatedExtendedReadsListIterator::freeSrcs() {
        freeSrc(rlPosSrc);
        freeSrc(rlOrgIdxSrc);
        freeSrc(rlRevCompSrc);
        freeSrc(rlMisCntSrc);
        freeSrc(rlMisSymSrc);
        freeSrc(rlMisOffSrc);
    }

    bool SeparatedExtendedReadsListIterator::moveNext() {
        if (++current < pgh->getReadsCount()) {
            uint_pg_len_max pos;
            uint_reads_cnt_std idx;
            uint8_t revComp = 0;
            PgSAHelpers::readValue<uint_pg_len_max>(*rlPosSrc, pos);
            PgSAHelpers::readValue<uint_reads_cnt_std>(*rlOrgIdxSrc, idx);
            if (rlRevCompSrc)
                PgSAHelpers::readValue<uint8_t>(*rlRevCompSrc, revComp);
            entry.advanceEntryByPosition(pos, idx, revComp==1);
            if (rlMisCntSrc) {
                uint8_t mismatchesCount;
                PgSAHelpers::readValue<uint8_t>(*rlMisCntSrc, mismatchesCount);
                for(uint8_t i = 0; i < mismatchesCount; i++) {
                    uint8_t mismatchCode;
                    uint_read_len_max mismatchOffset;
                    PgSAHelpers::readValue<uint8_t>(*rlMisSymSrc, mismatchCode);
                    PgSAHelpers::readValue<uint_read_len_max>(*rlMisOffSrc, mismatchOffset);
                    entry.addMismatch(mismatchOffset, mismatchCode);
                }
            }
            return true;
        }
        return false;
    }

    const PgTools::ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &
    SeparatedExtendedReadsListIterator::peekReadEntry() {
        return entry;
    }
}