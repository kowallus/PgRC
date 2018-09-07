#ifndef READSSETINTERFACE_H_INCLUDED
#define READSSETINTERFACE_H_INCLUDED

#include "../helper.h"
#include "../pgsaconfig.h"

using namespace PgSAReadsSet;

namespace PgSAReadsSet {

    template < typename uint_read_len, typename uint_reads_cnt >
    class ReadsSetInterface
    {
        public:

            virtual ~ReadsSetInterface() {};

            virtual bool isReadLengthConstantVirtual() = 0;
            virtual uint_read_len maxReadLengthVirtual() = 0;

            virtual uint_reads_cnt readsCountVirtual() = 0;

            virtual const string getReadVirtual(uint_reads_cnt) = 0;
            virtual uint_read_len readLengthVirtual(uint_reads_cnt) = 0;
    };

    typedef ReadsSetInterface< uint_read_len_std, uint_reads_cnt_std > StandardReadsSetInterface;
}

#endif // READSSETINTERFACE_H_INCLUDED
