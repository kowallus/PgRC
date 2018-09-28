#ifndef PSEUDOGENOMEINTERFACE_H_INCLUDED
#define PSEUDOGENOMEINTERFACE_H_INCLUDED

#include "../readsset/ReadsSetInterface.h"
#include "../readsset/ReadsSetBase.h"
#include "readslist/ReadsListInterface.h"

using namespace PgSAReadsSet;
using namespace PgSAIndex;

namespace PgSAIndex {
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass>
    class PseudoGenomeInterface: public ReadsSetInterface<uint_read_len, uint_reads_cnt>
    {
        public:
            virtual ~PseudoGenomeInterface() {};

            virtual const string getPseudoGenomeVirtual() = 0;

            template<class ReadsListClass>
            ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* getReadsList() { return static_cast<PseudoGenomeClass*>(this)->getReadsList(); };

            inline const char_pg getSymbol(const uint_pg_len posIdx) { return static_cast<PseudoGenomeClass*>(this)->getSymbolImpl(posIdx); };
            inline const string getPart(const uint_pg_len posIdx, const uint_pg_len length) { return static_cast<PseudoGenomeClass*>(this)->getPartImpl(posIdx, length); };
            inline const uint_pg_len getLength() { return static_cast<PseudoGenomeClass*>(this)->getLengthImpl(); };
    };

}

#endif // PSEUDOGENOMEINTERFACE_H_INCLUDED
