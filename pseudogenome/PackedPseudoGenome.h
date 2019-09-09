#ifndef PACKEDPSEUDOGENOME_H_INCLUDED
#define PACKEDPSEUDOGENOME_H_INCLUDED

#include "DefaultPseudoGenome.h"
#include "PackedPseudoGenomeBase.h"
#include "readslist/ReadsListTypes.h"
#include "../readsset/tools/SymbolsPackingFacility.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    const string PGTYPE_PACKED = "PACKED_PGEN";

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass >
    class PackedPseudoGenome: public PseudoGenomeInterface<uint_read_len, uint_reads_cnt, uint_pg_len, PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>>,
            public PackedPseudoGenomeBase
    {
        private:
            uint_pg_element* sequence;
            
            SymbolsPackingFacility* sPacker;
            
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* readsList = 0;

        public:

            PackedPseudoGenome(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* srcPseudoGenome, uchar symbolsPerElement);

            PackedPseudoGenome(uint_pg_len pgLength, std::istream& src);

            ~PackedPseudoGenome();

            static PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* castBase(PseudoGenomeBase* base);

            string getTypeID();

            void write(std::ostream& dest);
            
            uint_pg_len getElementsCountWithGuard();
            
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList();

//            uint_pg_len getLength() { return this->length; };

            char getSymbol(uint_pg_len pos);

            const uint_pg_element* getRawSuffix(const uint_pg_len rawPos);
            
            const uint_pg_element* getRawSuffix(const uint_reads_cnt readsListIdx, const uint_read_len pos);
            
            uint_read_len maxReadLength();

            uint_reads_cnt readsCount();

            bool isReadLengthConstant();

            const string getSuffix(const uint_pg_len pos, const uint_pg_len length);

            const string getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset, const uint_pg_len length);
          
            void getKmerByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength, char_pg* kmerPtr);
            
            const string getRead(uint_reads_cnt originalIdx);
            
            uint_read_len readLength(uint_reads_cnt originalIdx);

            const char_pg getSymbolImpl(const uint_pg_len posIdx);
            const string getPartImpl(const uint_pg_len posIdx, const uint_pg_len length);
            const uint_pg_len getLengthImpl();
            
            SymbolsPackingFacility* getSymbolsPacker();

            uint_read_len maxReadLengthVirtual();
            uint_reads_cnt readsCountVirtual();
            bool isReadLengthConstantVirtual();
            const string getReadVirtual(uint_reads_cnt i);
            uint_read_len readLengthVirtual(uint_reads_cnt i);

            const string getPseudoGenomeVirtual();
    };

    template <typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element>
    using PackedPseudoGenomeOfConstantLengthReadsType = PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>;
    
};

#endif // PACKEDPSEUDOGENOME_H_INCLUDED
