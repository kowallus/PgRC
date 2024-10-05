#include "SeparatedExtendedReadsList.h"

#include "../../coders/CodersLib.h"
#include "../SeparatedPseudoGenomeBase.h"

namespace PgTools {

    template <int maxMismatches>
    SeparatedExtendedReadsListIterator<maxMismatches>::SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix) {
        ownProps = true;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pseudoGenomePrefix, pgh, rsProp, plainTextReadMode);
        if (PgReadsSet::isReadLengthMin(pgh->getMaxReadLength()))
            PgHelpers::bytePerReadLengthMode = true;
        initSrcs();
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
                src = nullptr;
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
    void SeparatedExtendedReadsListIterator<maxMismatches>::freeSrc(istream *&src) {
        if (src) {
/*            if (fromFileMode())
                ((ifstream*) src)->close();*/
            delete (src);
            src = nullptr;
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
            PgHelpers::readValue<uint_reads_cnt_std>(*rlOrgIdxSrc, idx, plainTextReadMode);
            if (rlRevCompSrc)
                PgHelpers::readValue<uint8_t>(*rlRevCompSrc, revComp, plainTextReadMode);
            if (rlOffSrc) {
                uint_read_len_max offset = 0;
                PgHelpers::readReadLengthValue(*rlOffSrc, offset, plainTextReadMode);
                entry.advanceEntryByOffset(offset, idx, revComp == 1);
            } else {
                uint_pg_len_max pos = 0;
                PgHelpers::readValue<uint_pg_len_max>(*rlPosSrc, pos, plainTextReadMode);
                entry.advanceEntryByPosition(pos, idx, revComp == 1);
            }
            if (maxMismatches != 0 && rlMisCntSrc) {
                uint8_t mismatchesCount = 0;
                PgHelpers::readValue<uint8_t>(*rlMisCntSrc, mismatchesCount, plainTextReadMode);
                for(uint8_t i = 0; i < mismatchesCount; i++) {
                    uint8_t mismatchCode = 0;
                    uint_read_len_max mismatchOffset = 0;
                    PgHelpers::readValue<uint8_t>(*rlMisSymSrc, mismatchCode, plainTextReadMode);
                    if (rlMisOffSrc)
                        PgHelpers::readReadLengthValue(*rlMisOffSrc, mismatchOffset, plainTextReadMode);
                    else
                        PgHelpers::readReadLengthValue(*rlMisRevOffSrc, mismatchOffset, plainTextReadMode);
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
            PgHelpers::readArray(*(rl.rlOrgIdxSrc), res->orgIdx.data(), sizeof(uint_reads_cnt_std) * readsCount);
        } else
            res->orgIdx.resize(readsCount, 0);
        if (rl.rlRevCompSrc) {
            res->revComp.resize(readsCount);
            PgHelpers::readArray(*(rl.rlRevCompSrc), res->revComp.data(), sizeof(uint8_t) * readsCount);
        }
        res->pos.reserve(readsCount + 1);
        if (rl.rlOffSrc) {
            uint_pg_len_max pos = 0;
            uint_read_len_max offset = 0;
            for(uint_reads_cnt_max i = 0; i < readsCount; i++) {
                PgHelpers::readReadLengthValue(*(rl.rlOffSrc), offset, rl.plainTextReadMode);
                pos = pos + offset;
                res->pos.push_back(pos);
            }
        } else {
            res->pos.resize(readsCount);
            PgHelpers::readArray(*(rl.rlPosSrc), res->pos.data(), sizeof(uint_pg_len_max) * readsCount);
        }
        if (pgLengthPosGuard)
            res->pos.push_back(pgLengthPosGuard);
        if (!skipMismatches && rl.rlMisCntSrc) {
            uint_reads_cnt_max cumCount = 0;
            uint8_t mismatchesCount = 0;
            res->misCumCount.reserve(readsCount + 1);
            res->misCumCount.push_back(0);
            for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                PgHelpers::readValue<uint8_t>(*(rl.rlMisCntSrc), mismatchesCount, rl.plainTextReadMode);
                cumCount += mismatchesCount;
                res->misCumCount.push_back(cumCount);
            }
            res->misSymCode.resize(cumCount);
            PgHelpers::readArray(*(rl.rlMisSymSrc), res->misSymCode.data(), sizeof(uint8_t) * cumCount);
            res->misOff.resize(cumCount);
            bool misRevOffMode = rl.rlMisOffSrc == nullptr;
            PgHelpers::readArray(misRevOffMode?*(rl.rlMisRevOffSrc):*(rl.rlMisOffSrc), res->misOff.data(),
                    sizeof(uint8_t) * cumCount);
            if (misRevOffMode) {
                for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                    PgHelpers::convertMisRevOffsets2Offsets<uint8_t>(res->misOff.data() + res->misCumCount[i],
                            res->getMisCount(i), res->readLength);
                }
            }
        }
        cout << "Loaded Pg reads list containing " << readsCount << " reads." << endl;
        return res;
    }

    ExtendedReadsListWithConstantAccessOption* ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(
            istream& pgrcIn, PseudoGenomeHeader* pgh, ReadsSetProperties* rsProp, const string validationPgPrefix,
            PgRCParams* params, bool disableRevCompl, bool disableMismatches, bool separateFirstOffsetMode) {
        ExtendedReadsListWithConstantAccessOption *res = new ExtendedReadsListWithConstantAccessOption(pgh->getMaxReadLength());
        const uint_reads_cnt_max readsCount = pgh->getReadsCount();
        vector<string *> destStrings;
        string offStr, revCompStr, misCntStr, nonZeros, misSymCodeStr, misPropsStr;
        string srcsFirstStr[UINT8_MAX];
        string srcsStr[UINT8_MAX];
        if (!params->preserveOrderMode)
            destStrings.push_back(&offStr);
        if (!disableRevCompl)
            destStrings.push_back(&revCompStr);
        uint8_t mismatchesCountSrcsLimit = 0;
        vector<uint8_t> misCnt2SrcIdx;
        if (!disableMismatches) {
            destStrings.push_back(&misCntStr);
            if (params->isVersionAtLeast(1, 3))
                destStrings.push_back(&nonZeros);
            destStrings.push_back(&misSymCodeStr);
            if (params->isVersionAtLeast(1, 3))
                destStrings.push_back(&misPropsStr);
        }
        readCompressedCollectiveParallel(pgrcIn, destStrings);
        destStrings.clear();
        if (!disableMismatches) {
            istream* misPropsIn = params->isVersionAtLeast(1, 3) ? new istringstream(misPropsStr) : &pgrcIn;
            PgHelpers::readValue<uint8_t>(*misPropsIn, mismatchesCountSrcsLimit, false);
            misCnt2SrcIdx.resize(UINT8_MAX, mismatchesCountSrcsLimit);
            for (uint8_t m = 1; m < mismatchesCountSrcsLimit; m++)
                PgHelpers::readValue<uint8_t>(*misPropsIn, misCnt2SrcIdx[m], false);
            if (params->isVersionAtLeast(1, 3))
                delete misPropsIn;
            for (uint8_t m = 1; m <= mismatchesCountSrcsLimit; m++) {
                if (separateFirstOffsetMode)
                    destStrings.push_back(&srcsFirstStr[m]);
                if (!separateFirstOffsetMode || m > 1)
                    destStrings.push_back(&srcsStr[m]);
            }
        }
        readCompressedCollectiveParallel(pgrcIn, destStrings);
        moveStringToVector(offStr, res->off);
        moveStringToVector(revCompStr, res->revComp);
        if (!disableMismatches) {
            moveStringToVector(misCntStr, res->misCnt);
            if (params->isVersionAtLeast(1, 3)) {
                size_t j = 0;
                for (size_t i = 0; i < res->misCnt.size(); i++)
                    res->misCnt[i] = ((uint8_t) res->misCnt[i]) ? 0 : nonZeros[j++];
            }
            moveStringToVector(misSymCodeStr, res->misSymCode);
            vector<uint8_t> srcsFirst[UINT8_MAX];
            vector<uint8_t> srcs[UINT8_MAX];
            vector<uint_reads_cnt_std> srcFirstCounter(UINT8_MAX, 0);
            vector<uint_reads_cnt_std> srcRestCounter(UINT8_MAX, 0);
            for (uint8_t m = 1; m <= mismatchesCountSrcsLimit; m++) {
                if (separateFirstOffsetMode)
                    moveStringToVector(srcsFirstStr[m], srcsFirst[m]);
                if (!separateFirstOffsetMode || m > 1)
                    moveStringToVector(srcsStr[m], srcs[m]);
            }
            res->misOff.reserve(res->misSymCode.size());
            for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                uint8_t misCnt = res->misCnt[i];
                uint8_t srcIdx = misCnt2SrcIdx[misCnt];
                uint64_t misOffStartIdx = res->misOff.size();
                uint8_t m = 0;
                if (separateFirstOffsetMode && misCnt > 0) {
                    res->misOff.push_back(srcsFirst[srcIdx][srcFirstCounter[srcIdx]++]);
                    m++;
                }
                for (; m < misCnt; m++)
                    res->misOff.push_back(srcs[srcIdx][srcRestCounter[srcIdx]++]);
                PgHelpers::convertMisRevOffsets2Offsets<uint8_t>(res->misOff.data() + misOffStartIdx,
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
                for (uint8_t i = 0; i < mismatchesCount; i++) {
                    entry.addMismatch(misSymCode[curMisCumCount], misOff[curMisCumCount]);
                    ++curMisCumCount;
                }
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

    void ExtendedReadsListWithConstantAccessOption::enableConstantAccess(bool disableIterationMode, bool skipPositions) {
#pragma parallel omp
        {
#pragma omp single
            {
#pragma parallel task
                {
                    if (pos.empty() && !skipPositions) {
                        pos.reserve(readsCount + 1);
                        uint_pg_len_max currPos = 0;
                        for (uint_reads_cnt_max i = 0; i < readsCount; i++) {
                            currPos += off[i];
                            this->pos.push_back(currPos);
                        }
                        this->pos.push_back(this->pos.back() + this->readLength);
                    }
                    if (disableIterationMode)
                        off.clear();
                }
#pragma parallel task
                {
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
            }
        }
    }

    bool ExtendedReadsListWithConstantAccessOption::isConstantAccessEnalbed() {
        return !this->pos.empty();
    }

    template class SeparatedExtendedReadsListIterator<UINT8_MAX>;
    template class SeparatedExtendedReadsListIterator<0>;
}