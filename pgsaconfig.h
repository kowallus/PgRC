#ifndef PGSACONFIG_H_INCLUDED
#define PGSACONFIG_H_INCLUDED

#include <climits>

using namespace std;

typedef long long int int_max;
typedef unsigned long long int uint_max;

namespace PgSAReadsSet {

    typedef unsigned char uint_read_len_min;
    typedef unsigned short uint_read_len_std;
    typedef uint_read_len_std uint_read_len_max;

    inline bool isReadLengthMin(uint_max value) { return value <= UCHAR_MAX; }
    inline bool isReadLengthStd(uint_max value) { return !isReadLengthMin(value) && value <= USHRT_MAX; }

    // NOTE: (overlapping requires indexing from 1)
    typedef unsigned int uint_reads_cnt_std; // support up to 32 bits - 1
    typedef unsigned int uint_reads_cnt_max; // support up to 32 bits - 1
    //typedef unsigned long long int uint_reads_cnt_max; // support up to 64 bits - 1

    inline bool isReadsCountStd(uint_max value) { return value <= UINT_MAX - 1; }
    inline bool isReadsCountMax(uint_max value) { return value <= UINT_MAX - 1; }
//    inline bool isReadsCountMax(uint_max value) { return !isReadsCountStd(value) && value <= ULLONG_MAX - 1; }

    typedef unsigned char uint_symbols_cnt;
    
    typedef unsigned long long int uint_flatten_occurrence_max;
    
}

using namespace PgSAReadsSet;

namespace PgSAIndex {

    typedef unsigned int uint_pg_len_std;
    typedef unsigned long long int uint_pg_len_max;

    // NOTE: up to USHRT_MAX pseudogenome symbols are added to the end of pseudogenome
    inline bool isPGLengthStd(uint_max value) { return value <= UINT_MAX - USHRT_MAX; }
    inline bool isPGLengthMax(uint_max value) { return !isPGLengthStd(value) && value <= ULLONG_MAX - USHRT_MAX; }

    typedef char char_pg;
    typedef unsigned char uchar;
    typedef unsigned short ushort;
    
    inline uint_max oldestSetBit(uint_max value) {
        for(uint_max i = 64; i > 0; i--)
            if (((uint_max) 1L << (i-1)) & value)
                return i;
        return 0;
    }

    inline uchar bytesPerValue(uint_max value) {
        return (oldestSetBit(value) + 7) / 8;
    }
    
    template <typename uint_reads_cnt>
    struct BytesPerReadIndex {
        static const uchar minimum = 0;
        static const uchar standard = 0;
    };
    
    template <>
    struct BytesPerReadIndex<uint_reads_cnt_std> {
        static const uchar minimum = 3;
        static const uchar standard = 4;
    };

    template<typename uint_read_len, typename uint_reads_cnt>
    struct RPGOffset {
        uint_reads_cnt readListIndex;
        uint_read_len offset;
    };
    
    // types for packed symbols elements
    typedef uchar uint_ps_element_min;
    typedef ushort uint_ps_element_std;
    
}

#endif // PGSACONFIG_H_INCLUDED
