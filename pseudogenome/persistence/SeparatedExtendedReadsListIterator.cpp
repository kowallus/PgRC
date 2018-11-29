#include "SeparatedExtendedReadsListIterator.h"

namespace PgTools {

    template <int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>::SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix) {
        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pseudoGenomePrefix, pgh, plainTextReadMode);
        if (PgSAReadsSet::isReadLengthMin(pgh->getMaxReadLength()))
            PgSAHelpers::bytePerReadLengthMode = true;
        initSrcs();
    }

    template <int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>::~SeparatedExtendedReadsListIterator() {
        delete(pgh);
        freeSrcs();
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::initSrc(ifstream *&src, const string &fileSuffix) {
        src = new ifstream(pseudoGenomePrefix + fileSuffix, ios_base::in | ios_base::binary);
        if (src->fail()) {
            delete (src);
            src = 0;
        }
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::initSrcs() {
        initSrc(rlPosSrc, SeparatedPseudoGenomePersistence::READSLIST_POSITIONS_FILE_SUFFIX);
        initSrc(rlOrgIdxSrc, SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        initSrc(rlRevCompSrc, SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        initSrc(rlMisCntSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
        if (maxMismatches == 0 && rlMisCntSrc) {
            fprintf(stderr, "WARNING: mismatches unsupported in current routine working on %s Pg\n", pseudoGenomePrefix.c_str());
        }
        initSrc(rlMisSymSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
        initSrc(rlMisOffSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);

        initSrc(rlOffSrc, SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX);
        initSrc(rlMisRevOffSrc, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::freeSrc(ifstream *&src) {
        if (src) {
            src->close();
            delete (src);
            src = 0;
        }
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::freeSrcs() {
        freeSrc(rlPosSrc);
        freeSrc(rlOrgIdxSrc);
        freeSrc(rlRevCompSrc);
        freeSrc(rlMisCntSrc);
        freeSrc(rlMisSymSrc);
        freeSrc(rlMisOffSrc);

        freeSrc(rlOffSrc);
        freeSrc(rlMisRevOffSrc);
    }

    template <int maxMismatches>
    bool SeparatedExtendedReadsListIterator<maxMismatches>::moveNext() {
        if (++current < pgh->getReadsCount()) {
            uint_reads_cnt_std idx = 0;
            uint8_t revComp = 0;
            PgSAHelpers::readValue<uint_reads_cnt_std>(*rlOrgIdxSrc, idx, plainTextReadMode);
            if (rlRevCompSrc)
                PgSAHelpers::readValue<uint8_t>(*rlRevCompSrc, revComp, plainTextReadMode);
            if (rlOffSrc) {
                uint_read_len_max offset = 0;
                PgSAHelpers::readReadLengthValue(*rlOffSrc, offset, plainTextReadMode);
                entry.advanceEntryByOffset(offset, idx, revComp == 1);
            } else {
                uint_pg_len_max pos = 0;
                PgSAHelpers::readValue<uint_pg_len_max>(*rlPosSrc, pos, plainTextReadMode);
                entry.advanceEntryByPosition(pos, idx, revComp == 1);
            }
            if (maxMismatches != 0 && rlMisCntSrc) {
                uint8_t mismatchesCount = 0;
                PgSAHelpers::readValue<uint8_t>(*rlMisCntSrc, mismatchesCount, plainTextReadMode);
                for(uint8_t i = 0; i < mismatchesCount; i++) {
                    uint8_t mismatchCode = 0;
                    uint_read_len_max mismatchOffset = 0;
                    PgSAHelpers::readValue<uint8_t>(*rlMisSymSrc, mismatchCode, plainTextReadMode);
                    if (rlMisOffSrc)
                        PgSAHelpers::readReadLengthValue(*rlMisOffSrc, mismatchOffset, plainTextReadMode);
                    else
                        PgSAHelpers::readReadLengthValue(*rlMisRevOffSrc, mismatchOffset, plainTextReadMode);
                    entry.addMismatch(mismatchCode, mismatchOffset);
                }
                if (!rlMisOffSrc)
                    convertMisOffsets2RevOffsets(entry.mismatchOffset, entry.mismatchesCount, pgh->getMaxReadLength());
            }
            return true;
        }
        return false;
    }

    template <int maxMismatches>
    PgTools::ReadsListEntry<maxMismatches, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &
    SeparatedExtendedReadsListIterator<maxMismatches>::peekReadEntry() {
        return entry;
    }

    template <int maxMismatches>
    bool SeparatedExtendedReadsListIterator<maxMismatches>::isRevCompEnabled() {
        return rlRevCompSrc;
    }

    template <int maxMismatches>
    bool SeparatedExtendedReadsListIterator<maxMismatches>::areMismatchesEnabled() {
        return rlMisCntSrc;
    }

    template class SeparatedExtendedReadsListIterator<UINT8_MAX>;
    template class SeparatedExtendedReadsListIterator<0>;
}