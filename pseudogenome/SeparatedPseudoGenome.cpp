#include "SeparatedPseudoGenome.h"

namespace PgTools {

    string &SeparatedPseudoGenome::getPgSequence() {
        return pgSequence;
    }

    ExtendedReadsListWithConstantAccessOption *SeparatedPseudoGenome::getReadsList() {
        return readsList;
    }

    SeparatedPseudoGenome::SeparatedPseudoGenome(uint_pg_len_max sequenceLength, ReadsSetProperties* properties):
              SeparatedPseudoGenomeBase(sequenceLength, properties),
              readsList(new ExtendedReadsListWithConstantAccessOption(properties->maxReadLength)){
        readsList->readsCount = properties->readsCount;
    }

    SeparatedPseudoGenome::SeparatedPseudoGenome(string &&pgSequence, ExtendedReadsListWithConstantAccessOption *readsList,
            ReadsSetProperties* properties)
            : SeparatedPseudoGenomeBase(pgSequence.length(), properties),
            pgSequence(std::move(pgSequence)), readsList(readsList) {
        if (readsList)
            readsList->readsCount = properties->readsCount;
    }


    SeparatedPseudoGenome::~SeparatedPseudoGenome() {
        disposeReadsList();
    }


    void SeparatedPseudoGenome::disposeReadsList() {
        if (readsList) {
            delete(readsList);
            readsList = 0;
        }
    }

    void SeparatedPseudoGenome::applyIndexesMapping(IndexesMapping *indexesMapping) {
        for(uint_reads_cnt_std& idx: readsList->orgIdx)
            idx = indexesMapping->getReadOriginalIndex(idx);
    }

    void SeparatedPseudoGenome::applyRevComplPairFile() {
        uint_reads_cnt_std readsCount = readsList->orgIdx.size();
        if (!readsList->isRevCompEnabled())
            readsList->revComp.resize(readsCount, false);
        for(uint_reads_cnt_std i = 0; i < readsCount; i++)
            if (readsList->orgIdx[i] % 2)
                readsList->revComp[i] = !readsList->revComp[i];
    }

    string SeparatedPseudoGenome::getTypeID() {
        return PGTYPE_SEPARATED;
    }

    void SeparatedPseudoGenome::write(std::ostream &dest) {
        fprintf(stderr, "ERROR: Separated Pseudogenome cannot be written to a ostream.");
        exit(EXIT_FAILURE);
    }

    const string SeparatedPseudoGenome::getPseudoGenomeVirtual() {
        return pgSequence;
    }

    const string SeparatedPseudoGenome::getRead(uint_reads_cnt_max idx) {
        string res;
        res.resize(this->readsList->readLength);
        getRead(idx, (char*) res.data());
        return res;
    }

    void SeparatedPseudoGenome::getRead_Unsafe(uint_reads_cnt_max idx, char *ptr) {
        getRead_RawSequence(idx, ptr);
        if (this->readsList->revComp[idx])
            PgSAHelpers::reverseComplementInPlace(ptr, this->readsList->readLength);
        for(uint8_t i = 0; i < this->readsList->getMisCount(idx); i++) {
            const uint8_t misPos = this->readsList->getMisOff(idx, i);
            ptr[misPos] = PgSAHelpers::code2mismatch(ptr[misPos],
                                               this->readsList->getMisSymCode(idx, i));
        }
    }

    void SeparatedPseudoGenome::getRead(uint_reads_cnt_max idx, char *ptr) {
        getRead_RawSequence(idx, ptr);
        if (this->readsList->isRevCompEnabled() && this->readsList->revComp[idx])
            PgSAHelpers::reverseComplementInPlace(ptr, this->readsList->readLength);
        if (this->readsList->areMismatchesEnabled()) {
            for (uint8_t i = 0; i < this->readsList->getMisCount(idx); i++) {
                const uint8_t misPos = this->readsList->getMisOff(idx, i);
                ptr[misPos] = PgSAHelpers::code2mismatch(ptr[misPos],
                                                         this->readsList->getMisSymCode(idx, i));
            }
        }
    }

    void SeparatedPseudoGenome::getNextRead_Unsafe(char *ptr, uint_pg_len_max pos) {
        getRawSequenceOfReadLength(ptr, pos);
        if (this->readsList->revComp[nextRlIdx])
            PgSAHelpers::reverseComplementInPlace(ptr, this->readsList->readLength);
        uint8_t mismatchesCount = this->readsList->misCnt[nextRlIdx];
        for (uint8_t i = 0; i < mismatchesCount; i++) {
            const uint8_t misPos = this->readsList->misOff[curMisCumCount];
            ptr[misPos] = PgSAHelpers::code2mismatch(ptr[misPos],
                                                     this->readsList->misSymCode[curMisCumCount++]);
        }
        nextRlIdx++;
    }

    void SeparatedPseudoGenome::getNextRead_Unsafe(char *ptr) {
        curPos += this->readsList->off[nextRlIdx];
        getNextRead_Unsafe(ptr, curPos);
    }

    void SeparatedPseudoGenome::getNextRead_RawSequence(char *ptr) {
        curPos += this->readsList->off[nextRlIdx++];
        getRawSequenceOfReadLength(ptr, curPos);
    }

    GeneratedSeparatedPseudoGenome::GeneratedSeparatedPseudoGenome(uint_pg_len_max sequenceLength,
                                                                   ReadsSetProperties *properties)
            : SeparatedPseudoGenome(sequenceLength, properties) {
        pgSequence.resize(sequenceLength);
        sequence = (char_pg*) pgSequence.data();
        readsList->off.reserve(properties->readsCount);
        readsList->orgIdx.reserve(properties->readsCount);
    }

    void GeneratedSeparatedPseudoGenome::append(const string &read, uint_read_len_max length, uint_read_len_max overlap,
                                                uint_reads_cnt_max orgIdx) {
        readsList->off.push_back(delta);
        readsList->orgIdx.push_back(orgIdx);

        delta = length - overlap;
        if (delta > 0) {
            strncpy(this->sequence + pos, read.data(), delta);
            pos += delta;
        }
    }

    void GeneratedSeparatedPseudoGenome::validate() {
    }
}