#ifndef PGTOOLS_STRANDDETECTORBASE_H
#define PGTOOLS_STRANDDETECTORBASE_H

#include <stdint-gcc.h>
#include <vector>

#include "../pgsaconfig.h"

using namespace std;
using namespace PgSAReadsSet;

namespace PgTools {

    class StrandDetectorBase {

    public:
        virtual ~StrandDetectorBase() {};

        virtual vector<int8_t>
        detectStrands(uint_read_len_max overlap_threshold, bool paired_reads = true, bool concatanated_readssrc = false, int8_t groups_limit = INT8_MAX) = 0;
    };

}

#endif //PGTOOLS_STRANDDETECTORBASE_H
