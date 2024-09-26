#ifndef PGTOOLS_PGRC_DECODER_H
#define PGTOOLS_PGRC_DECODER_H

#include "pgrc-params.h"
#include "pgrc-data.h"

#include <condition_variable>
#include <queue>
#include <mutex>
#include <thread>

namespace PgTools {

    class PgRCDecoder {
    private:

        PgRCParams *params;
        PgRCData data;

        // TESTING&REPORTING
        chrono::steady_clock::time_point start_t;

        template<typename uint_pg_len>
        void applyRevComplPairFileToPgs(vector<uint_pg_len> &orgIdx2PgPos);

        void loadAllPgs(istream &pgrcIn);
        void loadAllPgs();

        void writeAllReadsInSEModeSequential(const string &outPrefix) const;
        void writeAllReadsInPEModeEachFileSequential(const string &outPrefix) const;
        template<typename uint_pg_len>
        void writeAllReadsInORDMode(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos) const;

        static const int PARALLEL_DECODING_THREADS_LIMIT = 4;
        const size_t CHUNK_SIZE_IN_BYTES = 1 << 17;

        int decodingThreadsCount;
        const int QUEUE_LIMIT_PER_DECODING_THREAD = 2;
        std::mutex mut[2];
        std::queue<string> out_queue[2];
        std::condition_variable data_cond[2];

        void writeAllReadsInSEModeParallelWritingThreads(const string &outPrefix);
        void writeAllReadsInPEModeParallelChunks(const string &outPrefix);
        void writeAllReadsInPEModeEachFileSequentialWithWritingThreads(const string &outPrefix);
        template<typename uint_pg_len>
        void writeAllReadsInORDModeWithWritingThreads(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos);
        template<typename uint_pg_len>
        void writeAllReadsInORDModeParallelChunks(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos);

        void pushOutToQueue(string &out, bool pairFile = false);
        void finishWritingParallel(bool pairFile = false);
        void writeFromQueue(const string outFileName, bool pairFile = false);

        void validateAllPgs();
        void validatePgsOrder();

        void preparePgsForValidation() const;

        uint_reads_cnt_max dnaStreamSize() const;

        const vector<uint_reads_cnt_max> getAllPgsOrgIdxs2RlIdx() const;
        const uint_reads_cnt_max getAllPgsOrgIdx(uint_reads_cnt_max idx) const;

    public:
        PgRCDecoder(PgRCParams* params): params(params){ }

        void decompressPgRC();
    };

}

#endif //PGTOOLS_PGRC_DECODER_H
