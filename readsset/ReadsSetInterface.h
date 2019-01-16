#ifndef READSSETINTERFACE_H_INCLUDED
#define READSSETINTERFACE_H_INCLUDED

#include "../utils/helper.h"
#include "../pgsaconfig.h"

using namespace PgSAReadsSet;

namespace PgSAReadsSet {

    template < typename uint_read_len, typename uint_reads_cnt >
    class ReadsSetInterface
    {
        public:

            virtual ~ReadsSetInterface() {};

            virtual bool isReadLengthConstant() = 0;
            virtual uint_read_len maxReadLength() = 0;

            virtual uint_reads_cnt readsCount() = 0;

            virtual const string getRead(uint_reads_cnt) = 0;
            virtual uint_read_len readLength(uint_reads_cnt) = 0;
    };

    typedef ReadsSetInterface< uint_read_len_std, uint_reads_cnt_std > StandardReadsSetInterface;
}

#endif // READSSETINTERFACE_H_INCLUDED
