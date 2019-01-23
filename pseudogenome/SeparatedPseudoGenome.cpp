#include "SeparatedPseudoGenome.h"

namespace PgTools {

    string &SeparatedPseudoGenome::getPgSequence() {
        return pgSequence;
    }

    ConstantAccessExtendedReadsList *SeparatedPseudoGenome::getReadsList() {
        return readsList;
    }

    SeparatedPseudoGenome::SeparatedPseudoGenome(uint_pg_len_max sequenceLength, ReadsSetProperties* properties):
              PseudoGenomeBase(sequenceLength, properties),
              readsList(new ConstantAccessExtendedReadsList(properties->maxReadLength)){
        readsList->readsCount = properties->readsCount;
    }

    SeparatedPseudoGenome::SeparatedPseudoGenome(string &&pgSequence, ConstantAccessExtendedReadsList *readsList,
            ReadsSetProperties* properties)
            : PseudoGenomeBase(pgSequence.length(), properties),
            pgSequence(std::move(pgSequence)), readsList(readsList) {
        readsList->readsCount = properties->readsCount;
    }


    SeparatedPseudoGenome::~SeparatedPseudoGenome() {
        delete(readsList);
    }

    void SeparatedPseudoGenome::applyIndexesMapping(IndexesMapping *indexesMapping) {
        for(uint_reads_cnt_std& idx: readsList->orgIdx)
            idx = indexesMapping->getReadOriginalIndex(idx);
    }

    void SeparatedPseudoGenome::applyRevComplPairFile() {
        uint_reads_cnt_std readsCount = readsList->orgIdx.size();
        for(uint_reads_cnt_std i = 0; i < readsCount; i++)
            readsList->revComp[i] = readsList->orgIdx[i]%2?!readsList->revComp[i]:readsList->revComp[i];
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

    GeneratedSeparatedPseudoGenome::GeneratedSeparatedPseudoGenome(uint_pg_len_max sequenceLength,
                                                                   ReadsSetProperties *properties)
            : SeparatedPseudoGenome(sequenceLength, properties) {
        pgSequence.resize(sequenceLength);
        sequence = (char_pg*) pgSequence.data();
        readsList->pos.reserve(properties->readsCount);
        readsList->orgIdx.reserve(properties->readsCount);
        readsList->revComp.reserve(properties->readsCount);
        readsList->misCumCount.reserve(properties->readsCount);
    }

    void GeneratedSeparatedPseudoGenome::append(const string &read, uint_read_len_max length, uint_read_len_max overlap,
                                                uint_reads_cnt_max orgIdx) {
        readsList->pos.push_back(pos);
        readsList->orgIdx.push_back(orgIdx);

        uint_read_len_max len = length - overlap;
        if (len > 0) {
            strncpy(this->sequence + pos, read.data(), len);
            pos += len;
        }
    }

    void GeneratedSeparatedPseudoGenome::validate() {
    }
}