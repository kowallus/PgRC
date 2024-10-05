#ifndef DEFAULTPSEUDOGENOME_H_INCLUDED
#define DEFAULTPSEUDOGENOME_H_INCLUDED

#include "PseudoGenomeInterface.h"
#include "readslist/ReadsListTypes.h"
#include "PseudoGenomeBase.h"

using namespace PgReadsSet;
using namespace PgHelpers;

namespace PgIndex {

    const string PGTYPE_DEFAULT = "DEFAULT_PGEN";

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    class DefaultPseudoGenome: public PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len, DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>>,
                               public PseudoGenomeBase
    {
        protected:
            char_pg* sequence;

            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* readsList = 0;

        public:

            DefaultPseudoGenome(uint_pg_len pgLength, ReadsSetProperties* properties);
            DefaultPseudoGenome(uint_pg_len pgLength, std::istream& src);
            ~DefaultPseudoGenome() override;
            
            static DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* castBase(PseudoGenomeBase* base);

            string getTypeID() override;

            void write(std::ostream& dest) override;
            
            bool validateUsing(DefaultReadsSet* readsSet) override;
            
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList();
            void unmanageReadsList() {
                readsList = 0;
            };

            uint_pg_len getLengthWithGuard();
            
            char getSymbol(uint_pg_len pos);

            const char_pg* getSuffix(const uint_pg_len pos);
            
            const char_pg* getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset);

            const char_pg* getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos);
            
            uint_read_len maxReadLength() override;

            uint_reads_cnt readsCount() override;

            bool isReadLengthConstant() override;

            const string getRead(const uint_reads_cnt originalIdx) override;
            
            uint_read_len readLength(const uint_reads_cnt originalIdx) override;

            const char_pg getSymbolImpl(const uint_pg_len posIdx);
            const string getPartImpl(const uint_pg_len posIdx, const uint_pg_len length);
            const uint_pg_len getLengthImpl();

            uint_read_len maxReadLengthVirtual();
            uint_reads_cnt readsCountVirtual();
            bool isReadLengthConstantVirtual();
            const string getReadVirtual(uint_reads_cnt i);
            uint_read_len readLengthVirtual(uint_reads_cnt i);

            const string getPseudoGenomeVirtual() override;
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    using DefaultPseudoGenomeOfConstantLengthReadsType = DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass >
    class GeneratedPseudoGenome: public DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass >
    {
        private:
            uint_pg_len pos = 0;
            GeneratedReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>* genReadsList = 0;

            const uint_read_len duplicateFilterKmerLength = 11;
            
            inline uint_max pgSuffixToTableIdx(const char_pg* suffix);
            
            uint_max getTableLength();;
            
        public:
            char_pg *getSequencePtr() const;

            GeneratedPseudoGenome(uint_pg_len sequenceLength, ReadsSetProperties* properties);

            ~GeneratedPseudoGenome() override;

            void append(const string& read, uint_read_len length, uint_read_len overlap, uint_reads_cnt orgIdx);
            void append(uint_read_len_max length, uint_read_len_max overlap, uint_reads_cnt_max orgIdx);

            void validate();
            
            void buildRepetitiveReadsFilter() override;
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len>
    using GeneratedPseudoGenomeOfConstantLengthReadsType = GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;

};

#endif // DEFAULTPSEUDOGENOME_H_INCLUDED
