#ifndef READSLISTTYPES_H_INCLUDED
#define READSLISTTYPES_H_INCLUDED

#include "ListOfConstantLengthReads.h"

namespace PgIndex {

    class WrongListOfConstantLengthReads { };
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
    struct ListOfConstantLengthReadsTypeTemplate {
        typedef WrongListOfConstantLengthReads Type; 
    };

    template<typename uint_read_len>
    struct ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt_std, uint_pg_len_std> {
        typedef ListOfConstantLengthReads<uint_read_len, uint_reads_cnt_std, uint_pg_len_std, 9, 4> Type;
    };

    template<typename uint_read_len>
    struct ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt_std, uint_pg_len_max> {
        typedef ListOfConstantLengthReads<uint_read_len, uint_reads_cnt_std, uint_pg_len_max, 13, 8> Type;
    };

}

#endif // READSLISTTYPES_H_INCLUDED
