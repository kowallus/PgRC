#include "SeparatedExtendedReadsList.h"

#include "../../utils/LzmaLib.h"
#include "../SeparatedPseudoGenomeBase.h"

namespace PgTools {

    template <int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>::SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix) {
        ownProps = true;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pseudoGenomePrefix, pgh, rsProp, plainTextReadMode);
        if (PgSAReadsSet::isReadLengthMin(pgh->getMaxReadLength()))
            PgSAHelpers::bytePerReadLengthMode = true;
        initSrcs();
    }

    template<int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>::SeparatedExtendedReadsListIterator(istream& pgrcIn,
            PseudoGenomeHeader* pgh, ReadsSetProperties* rsProp, const string pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix), pgh(pgh), rsProp(rsProp),
            plainTextReadMode(false) {
        ownProps = false;
        if (PgSAReadsSet::isReadLengthMin(pgh->getMaxReadLength()))
            PgSAHelpers::bytePerReadLengthMode = true;
        decompressSrcs(pgrcIn);
    }

    template<int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>* SeparatedExtendedReadsListIterator<maxMismatches>::getIterator(const string &pseudoGenomePrefix) {
        return new SeparatedExtendedReadsListIterator(pseudoGenomePrefix);
    }

    template <int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>::~SeparatedExtendedReadsListIterator() {
        if (ownProps) {
            delete (pgh);
            delete (rsProp);
        }
        freeSrcs();
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::initSrc(istream *&src, const string &fileSuffix) {
        if (!pseudoGenomePrefix.empty()) {
            src = new ifstream(pseudoGenomePrefix + fileSuffix, ios_base::in | ios_base::binary);
            if (src->fail()) {
                delete (src);
                src = 0;
            }
        }
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::initSrcs() {
        initSrc(rlPosSrc, SeparatedPseudoGenomeBase::READSLIST_POSITIONS_FILE_SUFFIX);
        initSrc(rlOrgIdxSrc, SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        initSrc(rlRevCompSrc, SeparatedPseudoGenomeBase::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        initSrc(rlMisCntSrc, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
        if (maxMismatches == 0 && rlMisCntSrc) {
            fprintf(stderr, "WARNING: mismatches unsupported in current routine working on %s Pg\n", pseudoGenomePrefix.c_str());
        }
        initSrc(rlMisSymSrc, SeparatedPseudoGenomeBase::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
        initSrc(rlMisOffSrc, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);

        initSrc(rlOffSrc, SeparatedPseudoGenomeBase::READSLIST_OFFSETS_FILE_SUFFIX);
        initSrc(rlMisRevOffSrc, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::decompressSrc(istream *&src, istream &pgrcIn,
                                                                          SymbolsPackingFacility<uint8_t> *symPacker) {
        uint64_t resLength = 0;
        if (symPacker)
            PgSAHelpers::readValue<uint64_t>(pgrcIn, resLength, false);
        string tmp;
        readCompressed(pgrcIn, tmp);
        if (symPacker)
            tmp = symPacker->reverseSequence((uint8_t*) tmp.data(), 0, resLength);
        src = new istringstream(tmp);
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::decompressMisRevOffSrc(istream &pgrcIn, bool transposeMode) {
        uint8_t mismatchesCountSrcsLimit = 0;
        PgSAHelpers::readValue<uint8_t>(pgrcIn, mismatchesCountSrcsLimit, false);
        vector<uint8_t> misCnt2SrcIdx(UINT8_MAX, mismatchesCountSrcsLimit);
        if (mismatchesCountSrcsLimit == 1) {
            decompressSrc(rlMisRevOffSrc, pgrcIn);
            return;
        }
        for(uint8_t m = 1; m < mismatchesCountSrcsLimit; m++)
            PgSAHelpers::readValue<uint8_t>(pgrcIn, misCnt2SrcIdx[m], false);

        istream* srcs[UINT8_MAX];
        for(uint8_t m = 1; m <= mismatchesCountSrcsLimit; m++) {
            cout << (int) m << ": ";
            decompressSrc(srcs[m], pgrcIn);
        }

        if (transposeMode) {
            for (uint8_t s = 1; s < mismatchesCountSrcsLimit; s++) {
                if (misCnt2SrcIdx[s] == misCnt2SrcIdx[s - 1] || misCnt2SrcIdx[s] == misCnt2SrcIdx[s + 1])
                    continue;
                string matrix = ((istringstream *) srcs[s])->str();
                uint64_t readsCount = matrix.size() / s / (bytePerReadLengthMode ? 1 : 2);
                if (bytePerReadLengthMode)
                    ((istringstream *) srcs[s])->str(transpose<uint8_t>(matrix, s, readsCount));
                else
                    ((istringstream *) srcs[s])->str(transpose<uint16_t>(matrix, s, readsCount));
            }
        }

        ostringstream misRevOffDest;
        uint8_t misCnt = 0;
        uint16_t revOff = 0;
        for(uint_reads_cnt_max i = 0; i < this->rsProp->readsCount; i++) {
            PgSAHelpers::readValue<uint8_t>(*rlMisCntSrc, misCnt, false);
            for(uint8_t m = 0; m < misCnt; m++) {
                PgSAHelpers::readReadLengthValue(*srcs[misCnt2SrcIdx[misCnt]], revOff, false);
                PgSAHelpers::writeReadLengthValue(misRevOffDest, revOff);
            }
        }
        rlMisCntSrc->seekg(0, ios::beg);
        rlMisRevOffSrc = new istringstream(misRevOffDest.str());
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::decompressSrcs(istream &pgrcIn) {
        decompressSrc(rlOffSrc, pgrcIn);
        decompressSrc(rlRevCompSrc, pgrcIn);//, &SymbolsPackingFacility<uint8_t>::BinaryPacker);
        decompressSrc(rlMisCntSrc, pgrcIn);
        if (maxMismatches == 0 && rlMisCntSrc) {
            fprintf(stderr, "WARNING: mismatches unsupported in current routine working on %s Pg\n", pseudoGenomePrefix.c_str());
        }
        decompressSrc(rlMisSymSrc, pgrcIn);//, &SymbolsPackingFacility<uint8_t>::QuaternaryPacker);
        decompressMisRevOffSrc(pgrcIn);

        initSrc(rlOrgIdxSrc, SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
    }

    template <int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::freeSrc(istream *&src) {
        if (src) {
/*            if (fromFileMode())
                ((ifstream*) src)->close();*/
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
                    convertMisRevOffsets2Offsets(entry.mismatchOffset, entry.mismatchesCount, pgh->getMaxReadLength());
            }
            return true;
        }
        return false;
    }

    template<int maxMismatches>
    void SeparatedExtendedReadsListIterator<maxMismatches>::rewind() {
        fprintf(stderr, "Error: Rewinding SeparatedExtendedReadsListIterator unimplemented.");
        exit(EXIT_FAILURE);
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

    ExtendedReadsListWithConstantAccessOption *ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(
            const string &pseudoGenomePrefix, uint_pg_len_max pgLengthPosGuard, bool skipMismatches) {
        DefaultSeparatedExtendedReadsListIterator rl(pseudoGenomePrefix);
        return loadConstantAccessExtendedReadsList(rl, pgLengthPosGuard, skipMismatches);
    }

    ExtendedReadsListWithConstantAccessOption *ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(
            DefaultSeparatedExtendedReadsListIterator &rl, uint_pg_len_max pgLengthPosGuard, bool skipMismatches) {
        ExtendedReadsListWithConstantAccessOption *res = new ExtendedReadsListWithConstantAccessOption(rl.pgh->getMaxReadLength());
        if (rl.plainTextReadMode) {
            fprintf(stderr, "Unsupported text plain read mode in creating ExtendedReadsListWithConstantAccessOption for %s\n\n",
                    rl.pseudoGenomePrefix.c_str());
            exit(EXIT_FAILURE);
        }

        const uint_reads_cnt_max readsCount = rl.pgh->getReadsCount();
        if (rl.rlOrgIdxSrc) {
            res->orgIdx.resize(readsCount);
            PgSAHelpers::readArray(*(rl.rlOrgIdxSrc), res->orgIdx.data(), sizeof(uint_reads_cnt_std) * readsCount);
        } else
            res->orgIdx.resize(readsCount, 0);
        if (rl.rlRevCompSrc) {
            res->revComp.resize(readsCount);
            PgSAHelpers::readArray(*(rl.rlRevCompSrc), res->revComp.data(), sizeof(uint8_t) * readsCount);
        }
        res->pos.reserve(readsCount + 1);
        if (rl.rlOffSrc) {
            uint_pg_len_max pos = 0;
            uint_read_len_max offset = 0;
            for(uint_reads_cnt_max i = 0; i < readsCount; i++) {
                PgSAHelpers::readReadLengthValue(*(rl.rlOffSrc), offset, rl.plainTextReadMode);
                pos = pos + offset;
                res->pos.push_back(pos);
            }
        } else {
            res->pos.resize(readsCount);
            PgSAHelpers::readArray(*(rl.rlPosSrc), res->pos.data(), sizeof(uint_pg_len_max) * readsCount);
        }
        if (pgLengthPosGuard)
            res->pos.push_back(pgLengthPosGuard);
        if (!skipMismatches && rl.rlMisCntSrc) {
            uint_reads_cnt_max cumCount = 0;
            uint8_t mismatchesCount = 0;
            res->misCumCount.reserve(readsCount + 1);
            res->misCumCount.push_back(0);
            for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                PgSAHelpers::readValue<uint8_t>(*(rl.rlMisCntSrc), mismatchesCount, rl.plainTextReadMode);
                cumCount += mismatchesCount;
                res->misCumCount.push_back(cumCount);
            }
            res->misSymCode.resize(cumCount);
            PgSAHelpers::readArray(*(rl.rlMisSymSrc), res->misSymCode.data(), sizeof(uint8_t) * cumCount);
            res->misOff.resize(cumCount);
            bool misRevOffMode = rl.rlMisOffSrc == 0;
            PgSAHelpers::readArray(misRevOffMode?*(rl.rlMisRevOffSrc):*(rl.rlMisOffSrc), res->misOff.data(),
                    sizeof(uint8_t) * cumCount);
            if (misRevOffMode) {
                for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                    PgSAHelpers::convertMisRevOffsets2Offsets<uint8_t>(res->misOff.data() + res->misCumCount[i],
                            res->getMisCount(i), res->readLength);
                }
            }
        }
        cout << "Loaded Pg reads list containing " << readsCount << " reads." << endl;
        return res;
    }

    ExtendedReadsListWithConstantAccessOption* ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(
            istream& pgrcIn, PseudoGenomeHeader* pgh, ReadsSetProperties* rsProp, const string validationPgPrefix,
            bool preserveOrderMode, bool disableRevCompl, bool disableMismatches) {
        ExtendedReadsListWithConstantAccessOption *res = new ExtendedReadsListWithConstantAccessOption(pgh->getMaxReadLength());
        const uint_reads_cnt_max readsCount = pgh->getReadsCount();
        if (!preserveOrderMode)
            readCompressed(pgrcIn, res->off);
        if (!disableRevCompl)
            readCompressed(pgrcIn, res->revComp);
        if (!disableMismatches) {
            readCompressed(pgrcIn, res->misCnt);
            readCompressed(pgrcIn, res->misSymCode);
            uint8_t mismatchesCountSrcsLimit = 0;
            PgSAHelpers::readValue<uint8_t>(pgrcIn, mismatchesCountSrcsLimit, false);
            vector<uint8_t> misCnt2SrcIdx(UINT8_MAX, mismatchesCountSrcsLimit);
            for (uint8_t m = 1; m < mismatchesCountSrcsLimit; m++)
                PgSAHelpers::readValue<uint8_t>(pgrcIn, misCnt2SrcIdx[m], false);
            vector<uint8_t> srcs[UINT8_MAX];
            vector<uint_reads_cnt_std> srcCounter(UINT8_MAX, 0);
            for (uint8_t m = 1; m <= mismatchesCountSrcsLimit; m++) {
                cout << (int) m << ": ";
                readCompressed(pgrcIn, srcs[m]);
            }
            res->misOff.reserve(res->misSymCode.size());
            for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                uint8_t misCnt = res->misCnt[i];
                uint8_t srcIdx = misCnt2SrcIdx[misCnt];
                uint64_t misOffStartIdx = res->misOff.size();
                for (uint8_t m = 0; m < misCnt; m++)
                    res->misOff.push_back(srcs[srcIdx][srcCounter[srcIdx]++]);
                PgSAHelpers::convertMisRevOffsets2Offsets<uint8_t>(res->misOff.data() + misOffStartIdx,
                                                                   res->misCnt[i], res->readLength);
            }
        }
        if (!validationPgPrefix.empty()) {
            std::ifstream in((validationPgPrefix + SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX).c_str(), std::ifstream::binary);
            res->orgIdx.resize(readsCount);
            readArray(in, res->orgIdx.data(), readsCount * sizeof(uint_reads_cnt_std));
        }
        cout << "Loaded Pg reads list containing " << readsCount << " reads." << endl;
        return res;
    }

    bool ExtendedReadsListWithConstantAccessOption::moveNext() {
        if (++current < readsCount) {
            entry.advanceEntryByOffset(off[current], orgIdx[current], this->getRevComp(current));
            if (!misCnt.empty()) {
                uint8_t mismatchesCount = misCnt[current];
                for (uint8_t i = 0; i < mismatchesCount; i++)
                    entry.addMismatch(misSymCode[curMisCumCount], misOff[curMisCumCount++]);
            }
            return true;
        }
        return false;
    }

    void ExtendedReadsListWithConstantAccessOption::rewind() {
        current = -1;
        curMisCumCount = 0;
    }

    bool ExtendedReadsListWithConstantAccessOption::isRevCompEnabled() {
        return !this->revComp.empty();
    }

    bool ExtendedReadsListWithConstantAccessOption::areMismatchesEnabled() {
        return !this->misCumCount.empty() || !this->misOff.empty();
    }

    bool ExtendedReadsListWithConstantAccessOption::getRevComp(uint_reads_cnt_std idx) {
        return isRevCompEnabled()?revComp[idx]:false;
    }

    void ExtendedReadsListWithConstantAccessOption::enableConstantAccess(bool disableIterationMode) {
        pos.reserve(readsCount + 1);
        uint_pg_len_max currPos = 0;
        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
            currPos += off[i];
            this->pos.push_back(currPos);
        }
        this->pos.push_back(this->pos.back() + this->readLength);
        if (disableIterationMode)
            off.clear();
        if (!misCnt.empty()) {
            uint_reads_cnt_max cumCount = 0;
            misCumCount.reserve(readsCount + 1);
            misCumCount.push_back(0);
            for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                cumCount += misCnt[i];
                misCumCount.push_back(cumCount);
            }
            if (disableIterationMode)
                misCnt.clear();
        }
    }

    bool ExtendedReadsListWithConstantAccessOption::isConstantAccessEnalbed() {
        return !this->pos.empty();
    }

    template class SeparatedExtendedReadsListIterator<UINT8_MAX>;
    template class SeparatedExtendedReadsListIterator<0>;
}