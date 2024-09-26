#ifndef PGTOOLS_PGRC_ENCODER_H
#define PGTOOLS_PGRC_ENCODER_H

#include "../utils/helper.h"
#include "pgrc-params.h"
#include "pgrc-data.h"

namespace PgTools {

    class PgRCEncoder {
    private:

        PgRCParams* params;
        PgRCData data;
        fstream pgrcOut;

        // TESTING&REPORTING
        chrono::steady_clock::time_point start_t;
        chrono::steady_clock::time_point div_t;
        chrono::steady_clock::time_point pgDiv_t;
        chrono::steady_clock::time_point good_t;
        chrono::steady_clock::time_point match_t;
        chrono::steady_clock::time_point bad_t;
        chrono::steady_clock::time_point order_t;
        size_t pgRCSize = 0;

        void generateReport();

        // CHAIN METHODS
        void prepareChainData();

        void initCompressionParameters();

        void runQualityBasedDivision();

        void runPgGeneratorBasedReadsDivision();

        void runHQPgGeneration();

        void runMappingLQReadsOnHQPg();

        void runLQPgGeneration();

        void runNPgGeneration();

        void persistReadsQualityDivision();

        void persistHQPg();

        void compressMappedHQPgReadsList();

        void persistLQPg();

        void compressLQPgReadsList();

        void persistHQPgSequence();

        void persistNPg();

        void compressNPgReadsList();

        void compressMEMMappedPgSequences();

        void prepareForPgGeneratorBaseReadsDivision();

        void prepareForHqPgGeneration();

        void prepareForMappingLQReadsOnHQPg();

        void prepareForLQPgAndNPgGeneration();

        void persistMappedReadsQualityDivision();

        void prepareForPgMatching();

        void finalizeCompression();

    public:

        PgRCEncoder(PgRCParams* params): params(params){ }

        void executePgRCChain();

    };
}

#endif //PGTOOLS_PGRC_ENCODER_H
