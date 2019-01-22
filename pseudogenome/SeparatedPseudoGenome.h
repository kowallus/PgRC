#ifndef PGTOOLS_SEPARATEDPSEUDOGENOME_H
#define PGTOOLS_SEPARATEDPSEUDOGENOME_H

#include "readslist/SeparatedExtendedReadsList.h"

namespace PgTools {

    class SeparatedPseudoGenome {
    protected:
        string pgSequence;
        ConstantAccessExtendedReadsList* readsList;
    public:
        SeparatedPseudoGenome(uint_read_len_max readLength);

        SeparatedPseudoGenome(string &&pgSequence,
                              ConstantAccessExtendedReadsList *readsList);

        string &getPgSequence();
        ConstantAccessExtendedReadsList *getReadsList();

        static SeparatedPseudoGenome* loadFromFile(string pgPrefix);

        virtual ~SeparatedPseudoGenome();
    };

    class GeneratedSeparatedPseudoGenome: public SeparatedPseudoGenome {
    private:
        uint_pg_len_max pos = 0;
        char_pg* sequence;

    public:
        GeneratedSeparatedPseudoGenome(uint_pg_len_max sequenceLength, ReadsSetProperties* properties);

        void append(const string& read, uint_read_len_max length, uint_read_len_max overlap, uint_reads_cnt_max orgIdx);

        void validate();
    };


}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOME_H
