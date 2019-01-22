#include "SeparatedPseudoGenome.h"

namespace PgTools {

    string &SeparatedPseudoGenome::getPgSequence() {
        return pgSequence;
    }

    ConstantAccessExtendedReadsList *SeparatedPseudoGenome::getReadsList() {
        return readsList;
    }

    SeparatedPseudoGenome::SeparatedPseudoGenome(uint_read_len_max readLength):
              readsList(new ConstantAccessExtendedReadsList(readLength)){
    }

    SeparatedPseudoGenome::SeparatedPseudoGenome(string &&pgSequence, ConstantAccessExtendedReadsList *readsList)
            : pgSequence(std::move(pgSequence)), readsList(readsList) {}


    SeparatedPseudoGenome::~SeparatedPseudoGenome() {
        delete(readsList);
    }

    SeparatedPseudoGenome* SeparatedPseudoGenome::loadFromFile(string pgPrefix){
        PseudoGenomeHeader *pgh = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgPrefix, pgh, plainTextReadMode);
        string pgSequence = SeparatedPseudoGenomePersistence::getPseudoGenome(pgPrefix);
        ConstantAccessExtendedReadsList* caeRl =
                DefaultSeparatedExtendedReadsList::loadConstantAccessExtendedReadsList(pgPrefix, pgh->getPseudoGenomeLength());
        SeparatedPseudoGenome* spg = new SeparatedPseudoGenome(std::move(pgSequence), caeRl);
        return spg;
    }

    GeneratedSeparatedPseudoGenome::GeneratedSeparatedPseudoGenome(uint_pg_len_max sequenceLength,
                                                                   ReadsSetProperties *properties)
            : SeparatedPseudoGenome(properties->maxReadLength) {
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