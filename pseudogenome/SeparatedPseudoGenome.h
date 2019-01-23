#ifndef PGTOOLS_SEPARATEDPSEUDOGENOME_H
#define PGTOOLS_SEPARATEDPSEUDOGENOME_H

#include "readslist/SeparatedExtendedReadsList.h"

namespace PgTools {

    const string PGTYPE_SEPARATED = "SEPARATED_PGEN";

    class SeparatedPseudoGenome: public PseudoGenomeBase {
    protected:
        string pgSequence;
        ConstantAccessExtendedReadsList* readsList;
    public:
        SeparatedPseudoGenome(uint_pg_len_max sequenceLength, ReadsSetProperties* properties);

        SeparatedPseudoGenome(string &&pgSequence,
                              ConstantAccessExtendedReadsList *readsList, ReadsSetProperties* properties);

        string &getPgSequence();
        ConstantAccessExtendedReadsList *getReadsList();

        virtual ~SeparatedPseudoGenome();

        void applyIndexesMapping(IndexesMapping *indexesMapping);

        void applyRevComplPairFile();

        string getTypeID() override;

        void write(std::ostream &dest) override;

        const string getPseudoGenomeVirtual() override;

        void disposeReadsList();
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
