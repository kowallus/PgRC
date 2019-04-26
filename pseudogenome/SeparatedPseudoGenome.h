#ifndef PGTOOLS_SEPARATEDPSEUDOGENOME_H
#define PGTOOLS_SEPARATEDPSEUDOGENOME_H

#include "readslist/SeparatedExtendedReadsList.h"

#include "SeparatedPseudoGenomeBase.h"

namespace PgTools {

    const string PGTYPE_SEPARATED = "SEPARATED_PGEN";

    class SeparatedPseudoGenome: public SeparatedPseudoGenomeBase {
    protected:
        string pgSequence;
        ExtendedReadsListWithConstantAccessOption* readsList;

        // iteration variables
        uint64_t nextRlIdx = 0;
        uint64_t curPos = 0;
        uint64_t curMisCumCount = 0;

    public:
        SeparatedPseudoGenome(uint_pg_len_max sequenceLength, ReadsSetProperties* properties);

        SeparatedPseudoGenome(string &&pgSequence,
                              ExtendedReadsListWithConstantAccessOption *readsList, ReadsSetProperties* properties);

        string &getPgSequence();
        ExtendedReadsListWithConstantAccessOption *getReadsList();

        virtual ~SeparatedPseudoGenome();

        void applyIndexesMapping(IndexesMapping *indexesMapping);

        void applyRevComplPairFile();

        string getTypeID() override;

        void write(std::ostream &dest) override;

        const string getPseudoGenomeVirtual() override;

        void disposeReadsList();

        // read access
        inline void getRawSequenceOfReadLength(char *ptr, uint_pg_len_max pos) {
            memcpy((void*) ptr, (void*) (pgSequence.data() + pos), this->readsList->readLength);
        }

        const string getRead(uint_reads_cnt_max idx);
        inline void getRead_RawSequence(uint_reads_cnt_max idx, char *ptr) {
            getRawSequenceOfReadLength(ptr, this->readsList->pos[idx]);
        }
        void getRead_Unsafe(uint_reads_cnt_max idx, char *ptr);
        void getRead(uint_reads_cnt_max idx, char *ptr);

        // iteration routines
        void getNextRead_RawSequence(char *ptr);
        void getNextRead_Unsafe(char *ptr, uint_pg_len_max pos);
        void getNextRead_Unsafe(char *ptr);

        void rewind();

    };

    class GeneratedSeparatedPseudoGenome: public SeparatedPseudoGenome {
    private:
        uint_pg_len_max pos = 0;
        uint_read_len_max delta = 0;
        char_pg* sequence;

    public:
        GeneratedSeparatedPseudoGenome(uint_pg_len_max sequenceLength, ReadsSetProperties* properties);

        void append(const string& read, uint_read_len_max length, uint_read_len_max overlap, uint_reads_cnt_max orgIdx);

        void validate();
    };


}

#endif //PGTOOLS_SEPARATEDPSEUDOGENOME_H
