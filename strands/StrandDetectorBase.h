#ifndef PGTOOLS_STRANDDETECTORBASE_H
#define PGTOOLS_STRANDDETECTORBASE_H

#include <stdint-gcc.h>
#include <vector>

using namespace std;

namespace PgTools {

    class StrandDetectorBase {

    public:
        virtual ~StrandDetectorBase() {};

        virtual vector<int8_t>
        detectStrands(int8_t groups_limit = INT8_MAX, bool paired_reads = true, bool concatanated_readssrc = false) = 0;
    };

}

#endif //PGTOOLS_STRANDDETECTORBASE_H
