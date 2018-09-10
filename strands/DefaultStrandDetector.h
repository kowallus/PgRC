#ifndef PGTOOLS_DEFAULTSTRANDDETECTOR_H
#define PGTOOLS_DEFAULTSTRANDDETECTOR_H

#include "AbstractStrandDetector.h"

namespace PgTools {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class DefaultStrandDetector: public AbstractStrandDetector<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> {
    protected:
        void matchReadStrands() override;
    public:
        DefaultStrandDetector(ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList);
        ~DefaultStrandDetector();


    };

}

#endif //PGTOOLS_DEFAULTSTRANDDETECTOR_H
