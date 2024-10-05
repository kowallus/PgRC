#include "pgrc-decoder.h"

#include "../matching/SimplePgMatcher.h"

namespace PgTools {

    void PgRCDecoder::decompressPgRC() {
#ifdef DEVELOPER_BUILD
        dump_after_decompression = params->extraFilesForValidation;
        dump_after_decompression_prefix = params->pgRCFileName + "_dump_";
#endif
        start_t = chrono::steady_clock::now();
        string tmpDirectoryPath = params->pgRCFileName + "/";
        ifstream pgrcIn(params->pgRCFileName);
        char pgrc_mode;
        if (pgrcIn) {
            for (int i = 0; i < strlen(PGRC_HEADER); i++) {
                if (PGRC_HEADER[i] != pgrcIn.get()) {
                    fprintf(stderr, "Error processing header.\n");
                    exit(EXIT_FAILURE);
                }
            }
            pgrc_mode = pgrcIn.get();
            if (pgrc_mode == PGRC_VERSION_MODE) {
                params->pgrcVersionMajor = pgrcIn.get();
                params->pgrcVersionMinor = pgrcIn.get();
                params->pgrcVersionRevision = pgrcIn.get();
                if (params->pgrcVersionMajor > PGRC_VERSION_MAJOR || (params->pgrcVersionMajor == PGRC_VERSION_MAJOR &&
                                                                   params->pgrcVersionMinor > PGRC_VERSION_MINOR)) {
                    fprintf(stderr, "ERROR: Archive is packed with a newer PgRC version %d.%d.\n",
                            (int) params->pgrcVersionMajor, (int) params->pgrcVersionMinor);
                    exit(EXIT_FAILURE);
                }
                params->compressionLevel = pgrcIn.get();
                pgrc_mode = pgrcIn.get();
            } else
                SimplePgMatcher::MATCH_MARK = (char) 128;
            if (pgrc_mode < PGRC_SE_MODE || pgrc_mode > PGRC_MIN_PE_MODE) {
                fprintf(stderr, "Unsupported decompression mode: %d.\n", pgrc_mode);
                exit(EXIT_FAILURE);
            }
            params->separateNReads = (bool) pgrcIn.get();
            if (pgrc_mode == PGRC_PE_MODE || pgrc_mode == PGRC_ORD_PE_MODE)
                params->revComplPairFile = (bool) pgrcIn.get();

            pgrcIn >> tmpDirectoryPath;
            tmpDirectoryPath = tmpDirectoryPath + "/";
            pgrcIn.get();
        }

        params->pgMappedHqPrefix = tmpDirectoryPath + GOOD_INFIX;
        params->pgMappedLqPrefix = tmpDirectoryPath + BAD_INFIX;
        params->pgSeqFinalHqPrefix = tmpDirectoryPath + GOOD_INFIX;
        params->pgSeqFinalLqPrefix = tmpDirectoryPath + BAD_INFIX;
        params->pgNPrefix = tmpDirectoryPath + N_INFIX;

        params->preserveOrderMode = pgrc_mode == PGRC_ORD_SE_MODE || pgrc_mode == PGRC_ORD_PE_MODE;
        params->ignorePairOrderInformation = pgrc_mode == PGRC_MIN_PE_MODE;
        params->singleReadsMode = pgrc_mode == PGRC_SE_MODE || pgrc_mode == PGRC_ORD_SE_MODE;
        if (pgrcIn)
            loadAllPgs(pgrcIn);
        else
            loadAllPgs();
        cout << "... loaded Pgs (checkpoint: " << time_millis(start_t) << " msec.)" << endl;

        uint_reads_cnt_max nonNPgReadsCount = data.hqPg->getReadsSetProperties()->readsCount
                                              + data.lqPg->getReadsSetProperties()->readsCount;
        uint_reads_cnt_max nPgReadsCount = data.nPg ? data.nPg->getReadsSetProperties()->readsCount : 0;
        uint_reads_cnt_max readsTotalCount = nonNPgReadsCount + nPgReadsCount;
        if (params->revComplPairFile) {
            if (params->isJoinedPgLengthStd) {
                applyRevComplPairFileToPgs<uint_pg_len_std>(data.orgIdx2StdPgPos);
            } else {
                applyRevComplPairFileToPgs<uint_pg_len_max>(data.orgIdx2PgPos);
            }
        }
        if (params->srcFastqFile.empty()) {
            if (params->singleReadsMode && !params->preserveOrderMode) {
                writeAllReadsInSEModeParallelWritingThreads(params->pgRCFileName);
            } else if (!params->preserveOrderMode) {
                writeAllReadsInPEModeParallelChunks(params->pgRCFileName);
            } else if (params->isJoinedPgLengthStd) {
                writeAllReadsInORDModeParallelChunks<uint_pg_len_std>(params->pgRCFileName, data.orgIdx2StdPgPos);
            } else {
                writeAllReadsInORDModeParallelChunks<uint_pg_len_max>(params->pgRCFileName, data.orgIdx2PgPos);
            }
            cout << "Decompressed ";
        } else {
            preparePgsForValidation();
            validateAllPgs();
            validatePgsOrder();
            cout << "Validated ";
        }

        cout << readsTotalCount << " reads in " << time_millis(start_t) << " msec." << endl;

        data.disposeChainData();
    }

    static constexpr timespec WRITE_THREAD_SLEEP_TIME[] {{0, 1000L}};

    void PgRCDecoder::pushOutToQueue(string &out, bool pairFile) {
        int id = pairFile ? 1 : 0;
        while (out_queue[id].size() >= decodingThreadsCount * QUEUE_LIMIT_PER_DECODING_THREAD)
            nanosleep(WRITE_THREAD_SLEEP_TIME, nullptr);
        std::lock_guard<std::mutex> _(this->mut[id]);
        out_queue[id].push(std::move(out));
        data_cond[id].notify_one();
        out.resize(0);
        out.reserve(CHUNK_SIZE_IN_BYTES);
    }

    void PgRCDecoder::finishWritingParallel(bool pairFile) {
        int id = pairFile ? 1 : 0;
        std::lock_guard<std::mutex> _(this->mut[id]);
        out_queue[id].push(std::move(""));
        data_cond[id].notify_one();
    }

    void PgRCDecoder::writeFromQueue(const string outFileName, bool pairFile) {
        int id = pairFile ? 1 : 0;
        fstream fout(outFileName, ios_base::out | ios_base::binary | std::ios::trunc);
        string out;
        do {
            std::unique_lock<std::mutex> lk(mut[id]);
            data_cond[id].wait(
                    lk, [&, this] { return !out_queue[id].empty(); });
            out = out_queue[id].front();
            out_queue[id].pop();
            lk.unlock();
            PgHelpers::writeArray(fout, (void*) out.data(), out.size());
        } while (!out.empty());
        fout.close();
    }


    void PgRCDecoder::writeAllReadsInSEModeParallelWritingThreads(const string &outPrefix) {
        int threadsLimit = PARALLEL_DECODING_THREADS_LIMIT + 1;
        threadsLimit = PgHelpers::numberOfThreads < threadsLimit? PgHelpers::numberOfThreads : threadsLimit;
        decodingThreadsCount = threadsLimit - 1;
        if (decodingThreadsCount <= 0 || dnaStreamSize() <= CHUNK_SIZE_IN_BYTES) {
            writeAllReadsInSEModeSequential(outPrefix);
            return;
        }
        data.hqPg->getReadsList()->enableConstantAccess(true);
        data.lqPg->getReadsList()->enableConstantAccess(true);
        if (data.nPg) data.nPg->getReadsList()->enableConstantAccess(true);
        *PgHelpers::logout << "... enabled constant access (checkpoint: " << time_millis(start_t) << " msec.)" << endl;

        std::thread writing(&PgRCDecoder::writeFromQueue, this, outPrefix + "_out", false);
        string res, read;
        read.resize(params->readLength);
        uint64_t reads_per_chunk = CHUNK_SIZE_IN_BYTES / (params->readLength + 1);
        uint64_t hq_full_chunks_count = params->hqReadsCount / reads_per_chunk;
#pragma omp parallel for ordered schedule(static, 1) firstprivate(read, res) num_threads(decodingThreadsCount)
        for (int64_t c = 0; c <= hq_full_chunks_count; c++) {
            res.resize(0);
            res.reserve(CHUNK_SIZE_IN_BYTES);
            uint64_t i_guard = c == hq_full_chunks_count ? params->hqReadsCount : (c + 1) * reads_per_chunk;
            for (uint_reads_cnt_max i = c * reads_per_chunk; i < i_guard; i++) {
                data.hqPg->getRead_Unsafe(i, (char *) read.data());
                res.append(read);
                res.push_back('\n');
            }
            pushOutToQueue(res);
        }

        uint64_t lq_full_chunks_count = params->lqReadsCount / reads_per_chunk;
#pragma omp parallel for ordered schedule(static, 1) firstprivate(read, res) num_threads(decodingThreadsCount)
        for (int64_t c = 0; c <= lq_full_chunks_count; c++) {
            res.resize(0);
            res.reserve(CHUNK_SIZE_IN_BYTES);
            uint64_t i_guard = c == lq_full_chunks_count ? params->lqReadsCount : (c + 1) * reads_per_chunk;
            for (uint_reads_cnt_max i = c * reads_per_chunk; i < i_guard; i++) {
                data.lqPg->getRead_RawSequence(i, (char *) read.data());
                res.append(read);
                res.push_back('\n');
            }
            pushOutToQueue(res);
        }

        uint64_t n_full_chunks_count = params->nPgReadsCount / reads_per_chunk;
#pragma omp parallel for ordered schedule(static, 1) firstprivate(read, res) num_threads(decodingThreadsCount)
        for (int64_t c = 0; c <= n_full_chunks_count; c++) {
            res.resize(0);
            res.reserve(CHUNK_SIZE_IN_BYTES);
            uint64_t i_guard = c == n_full_chunks_count ? params->nPgReadsCount : (c + 1) * reads_per_chunk;
            for (uint_reads_cnt_max i = c * reads_per_chunk; i < i_guard; i++) {
                data.nPg->getRead_RawSequence(i, (char *) read.data());
                res.append(read);
                res.push_back('\n');
            }
            pushOutToQueue(res);
        }
        pushOutToQueue(res);

        *PgHelpers::logout << "... finished loading queue (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        finishWritingParallel();
        writing.join();
    }


    void PgRCDecoder::writeAllReadsInSEModeSequential(const string &outPrefix) const {
        fstream fout(outPrefix + "_out", ios_base::out | ios_base::binary | std::ios::trunc);
        string res, read;
        read.resize(params->readLength);
        uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
        uint64_t totalSize = dnaStreamSize();
        res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (params->readLength + 1));
        for (uint_reads_cnt_max i = 0; i < params->hqReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            data.hqPg->getNextRead_Unsafe((char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        for (uint_reads_cnt_max i = 0; i < params->lqReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            data.lqPg->getNextRead_RawSequence((char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        for (uint_reads_cnt_max i = 0; i < params->nPgReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            data.nPg->getNextRead_RawSequence((char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        fout << res;
        fout.close();
    }

    void PgRCDecoder::writeAllReadsInPEModeParallelChunks(const string &outPrefix) {
        int threadsLimit = PARALLEL_DECODING_THREADS_LIMIT + 1;
        threadsLimit = PgHelpers::numberOfThreads < threadsLimit? PgHelpers::numberOfThreads : threadsLimit;
        decodingThreadsCount = threadsLimit - 1;
        if (decodingThreadsCount < 4) {
            writeAllReadsInPEModeEachFileSequentialWithWritingThreads(outPrefix);
            return;
        }
        data.hqPg->getReadsList()->enableConstantAccess(true);
        data.lqPg->getReadsList()->enableConstantAccess(true);
        if (data.nPg) data.nPg->getReadsList()->enableConstantAccess(true);
        *PgHelpers::logout << "... enabled constant access (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        const uint8_t PE_PARTS_COUNT = 2;
        const uint64_t pe_reads_per_chunk = CHUNK_SIZE_IN_BYTES / (params->readLength + 1) * PE_PARTS_COUNT;
        const uint64_t full_chunks_count = params->readsTotalCount / pe_reads_per_chunk;

        for (uint8_t p = 0; p < PE_PARTS_COUNT; p++) {
            bool pairFile = p == 1;
            std::thread writing(&PgRCDecoder::writeFromQueue, this, outPrefix + "_out_" + (pairFile ? "2" : "1"), pairFile);
            string res, read;
            read.resize(params->readLength);
#pragma omp parallel for ordered schedule(static, 1) firstprivate(read, res) num_threads(decodingThreadsCount)
            for (int64_t c = 0; c <= full_chunks_count; c++) {
                res.resize(0);
                res.reserve(CHUNK_SIZE_IN_BYTES);
                char *readPtr = (char *) read.data();
                uint64_t i_guard = c == full_chunks_count ? params->readsTotalCount : p + (c + 1) * pe_reads_per_chunk;
                for (uint_reads_cnt_max i = p + c * pe_reads_per_chunk; i < i_guard; i += PE_PARTS_COUNT) {
                    uint_reads_cnt_std idx = data.rlIdxOrder[i];
                    if (idx < params->hqReadsCount)
                        data.hqPg->getRead(idx, readPtr);
                    else {
                        if (idx < params->nonNPgReadsCount)
                            data.lqPg->getRead_RawSequence(idx - params->hqReadsCount, readPtr);
                        else
                            data.nPg->getRead_RawSequence(idx - params->nonNPgReadsCount, readPtr);
                        if (p)
                            PgHelpers::reverseComplementInPlace(readPtr, params->readLength);
                    }
                    res.append(read);
                    res.push_back('\n');
                }
#pragma omp ordered
                pushOutToQueue(res, pairFile);
            }
            *PgHelpers::logout << "... finished loading queue " << (pairFile ? "2" : "1") << " (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
            finishWritingParallel(pairFile);
            writing.join();
        }
    }

    void PgRCDecoder::writeAllReadsInPEModeEachFileSequentialWithWritingThreads(const string &outPrefix) {
        if (PgHelpers::numberOfThreads < 4) {
            writeAllReadsInPEModeEachFileSequential(outPrefix);
            return;
        }
        decodingThreadsCount = 1;

        data.hqPg->getReadsList()->enableConstantAccess(true);
        data.lqPg->getReadsList()->enableConstantAccess(true);
        if (data.nPg) data.nPg->getReadsList()->enableConstantAccess(true);
        *PgHelpers::logout << "... enabled constant access (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        const uint8_t PE_PARTS_COUNT = 2;

#pragma omp parallel for
        for (uint8_t p = 0; p < PE_PARTS_COUNT; p++) {
            bool pairFile = p == 1;
            std::thread writing(&PgRCDecoder::writeFromQueue, this, outPrefix + "_out_" + (pairFile ? "2" : "1"), pairFile);
            string res, read;
            read.resize(params->readLength);
            char *readPtr = (char *) read.data();
            uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
            uint64_t totalSize = dnaStreamSize();
            res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (params->readLength + 1));
            uint_reads_cnt_max i = 0;
            for (i = p; i < params->readsTotalCount; i += PE_PARTS_COUNT) {
                if (res.size() > res_size_guard) {
                    pushOutToQueue(res, pairFile);
                    res.resize(0);
                }
                uint_reads_cnt_std idx = data.rlIdxOrder[i];
                if (idx < params->hqReadsCount)
                    data.hqPg->getRead(idx, readPtr);
                else {
                    if (idx < params->nonNPgReadsCount)
                        data.lqPg->getRead_RawSequence(idx - params->hqReadsCount, readPtr);
                    else
                        data.nPg->getRead_RawSequence(idx - params->nonNPgReadsCount, readPtr);
                    if (p)
                        PgHelpers::reverseComplementInPlace(readPtr, params->readLength);
                }
                res.append(read);
                res.push_back('\n');
            }
            pushOutToQueue(res, pairFile);
            *PgHelpers::logout << "... finished loading queue " << (pairFile ? "2" : "1") << " (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
            finishWritingParallel(pairFile);
            writing.join();
        }
    }

    void PgRCDecoder::writeAllReadsInPEModeEachFileSequential(const string &outPrefix) const {
        data.hqPg->getReadsList()->enableConstantAccess(true);
        data.lqPg->getReadsList()->enableConstantAccess(true);
        if (data.nPg) data.nPg->getReadsList()->enableConstantAccess(true);
        *PgHelpers::logout << "... enabled constant access (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        const uint8_t PE_PARTS_COUNT = 2;

#pragma omp parallel for
        for (uint8_t p = 0; p < PE_PARTS_COUNT; p++) {
            fstream fout(outPrefix + "_out" + "_" + toString(p + 1),
                         ios_base::out | ios_base::binary | std::ios::trunc);
            string res, read;
            read.resize(params->readLength);
            char *readPtr = (char *) read.data();
            uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
            uint64_t totalSize = dnaStreamSize();
            res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (params->readLength + 1));
            uint_reads_cnt_max i = 0;
            for (i = p; i < params->readsTotalCount; i += PE_PARTS_COUNT) {
                if (res.size() > res_size_guard) {
                    PgHelpers::writeArray(fout, (void*) res.data(), res.size());
                    res.resize(0);
                }
                uint_reads_cnt_std idx = data.rlIdxOrder[i];
                if (idx < params->hqReadsCount)
                    data.hqPg->getRead(idx, readPtr);
                else {
                    if (idx < params->nonNPgReadsCount)
                        data.lqPg->getRead_RawSequence(idx - params->hqReadsCount, readPtr);
                    else
                        data.nPg->getRead_RawSequence(idx - params->nonNPgReadsCount, readPtr);
                    if (p)
                        PgHelpers::reverseComplementInPlace(readPtr, params->readLength);
                }
                res.append(read);
                res.push_back('\n');
            }
            PgHelpers::writeArray(fout, (void*) res.data(), res.size());
            fout.close();
        }
    }

    template<typename uint_pg_len>
    void PgRCDecoder::writeAllReadsInORDModeParallelChunks(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos)  {
        uint8_t parts = params->singleReadsMode ? 1 : 2;
        int threadsLimit = PARALLEL_DECODING_THREADS_LIMIT + 1;
        threadsLimit = PgHelpers::numberOfThreads < threadsLimit? PgHelpers::numberOfThreads : threadsLimit;
        decodingThreadsCount = threadsLimit - 1;
        if (decodingThreadsCount < 2) {
            writeAllReadsInORDModeWithWritingThreads<uint_pg_len>(outPrefix, orgIdx2PgPos);
            return;
        }
        data.hqPg->getReadsList()->enableConstantAccess(true, true);
        data.rlIdxOrder.resize(params->readsTotalCount);
        uint_reads_cnt_max idx = 0;
        for (uint_reads_cnt_max i = 0; i < params->readsTotalCount; i++) {
            uint_pg_len pos = orgIdx2PgPos[i];
            if (pos < params->hqPgLen)
                data.rlIdxOrder[i] = idx++;
        }
        *PgHelpers::logout << "... enabled constant access (checkpoint: " << time_millis(start_t) << " msec.)" << endl;

        const uint64_t reads_per_chunk = CHUNK_SIZE_IN_BYTES / (params->readLength + 1);
        const uint64_t full_chunks_count = params->readsTotalCount / reads_per_chunk / parts;
        for (uint8_t p = 0; p < parts; p++) {
            bool pairFile = p == 1;
            std::thread writing(&PgRCDecoder::writeFromQueue, this, outPrefix + "_out" + (parts == 1 ? "" : (pairFile ? "_2" : "_1")), pairFile);
            string res, read;
            read.resize(params->readLength);
#pragma omp parallel for ordered schedule(static, 1) firstprivate(read, res) num_threads(decodingThreadsCount)
            for (int64_t c = 0; c <= full_chunks_count; c++) {
                res.resize(0);
                res.reserve(CHUNK_SIZE_IN_BYTES);
                char *readPtr = (char *) read.data();
                uint_reads_cnt_std begI = (params->readsTotalCount / parts) * p;
                uint_reads_cnt_std endI = (params->readsTotalCount / parts) * (p + 1);
                uint64_t i_guard = c == full_chunks_count ? endI : begI + (c + 1) * reads_per_chunk;
                for (uint_reads_cnt_max i = begI + c * reads_per_chunk; i < i_guard; i++) {
                    uint_pg_len pos = orgIdx2PgPos[i];
                    if (pos < params->hqPgLen) {
                        idx = data.rlIdxOrder[i];
                        data.hqPg->getRead_Unsafe(idx, pos, readPtr);
                    } else {
                        if (pos < params->nonNPgLen)
                            data.lqPg->getRawSequenceOfReadLength(readPtr, pos - params->hqPgLen);
                        else
                            data.nPg->getRawSequenceOfReadLength(readPtr, pos - params->nonNPgLen);
                        if (p)
                            PgHelpers::reverseComplementInPlace(readPtr, params->readLength);
                    }
                    res.append(read);
                    res.push_back('\n');
                }
#pragma omp ordered
                pushOutToQueue(res, pairFile);
            }
            pushOutToQueue(res, pairFile);
            *PgHelpers::logout << "... finished loading queue" <<  (parts == 1 ? "" : (pairFile ? "2 " : "1 ")) << " (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
            finishWritingParallel(pairFile);
            writing.join();
        }
    }

    template<typename uint_pg_len>
    void PgRCDecoder::writeAllReadsInORDModeWithWritingThreads(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos)  {
        uint8_t parts = params->singleReadsMode ? 1 : 2;

        if (PgHelpers::numberOfThreads < 2) {
            writeAllReadsInORDMode<uint_pg_len>(outPrefix, orgIdx2PgPos);
            return;
        }
        decodingThreadsCount = 1;

        for (uint8_t p = 0; p < parts; p++) {
            bool pairFile = p == 1;
            std::thread writing(&PgRCDecoder::writeFromQueue, this, outPrefix + "_out" + (parts == 1 ? "" : (pairFile ? "_2" : "_1")), pairFile);
            string res, read;
            read.resize(params->readLength);
            char *readPtr = (char *) read.data();
            uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
            uint64_t totalSize = dnaStreamSize();
            res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (params->readLength + 1));
            uint_reads_cnt_std endI = (params->readsTotalCount / parts) * (p + 1);
            for (uint_reads_cnt_max i = (params->readsTotalCount / parts) * p; i < endI; i++) {
                if (res.size() > res_size_guard) {
                    pushOutToQueue(res, pairFile);
                    res.resize(0);
                }
                uint_pg_len pos = orgIdx2PgPos[i];
                if (pos < params->hqPgLen)
                    data.hqPg->getNextRead_Unsafe(readPtr, pos);
                else {
                    if (pos < params->nonNPgLen)
                        data.lqPg->getRawSequenceOfReadLength(readPtr, pos - params->hqPgLen);
                    else
                        data.nPg->getRawSequenceOfReadLength(readPtr, pos - params->nonNPgLen);
                    if (p)
                        PgHelpers::reverseComplementInPlace(readPtr, params->readLength);
                }
                res.append(read);
                res.push_back('\n');
            }
            pushOutToQueue(res, pairFile);
            *PgHelpers::logout << "... finished loading queue" <<  (parts == 1 ? "" : (pairFile ? "2 " : "1 ")) << " (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
            finishWritingParallel(pairFile);
            writing.join();
        }
    }

    template<typename uint_pg_len>
    void PgRCDecoder::writeAllReadsInORDMode(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos) const {
        uint8_t parts = params->singleReadsMode ? 1 : 2;

        for (uint8_t p = 0; p < parts; p++) {
            fstream fout(outPrefix + "_out" + (params->singleReadsMode ? "" : ("_" + toString(p + 1))),
                         ios_base::out | ios_base::binary | std::ios::trunc);
            string res, read;
            read.resize(params->readLength);
            char *readPtr = (char *) read.data();
            uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
            uint64_t totalSize = dnaStreamSize();
            res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (params->readLength + 1));
            uint_reads_cnt_std endI = (params->readsTotalCount / parts) * (p + 1);
            for (uint_reads_cnt_max i = (params->readsTotalCount / parts) * p; i < endI; i++) {
                if (res.size() > res_size_guard) {
                    fout << res;
                    res.resize(0);
                }
                uint_pg_len pos = orgIdx2PgPos[i];
                if (pos < params->hqPgLen)
                    data.hqPg->getNextRead_Unsafe(readPtr, pos);
                else {
                    if (pos < params->nonNPgLen)
                        data.lqPg->getRawSequenceOfReadLength(readPtr, pos - params->hqPgLen);
                    else
                        data.nPg->getRawSequenceOfReadLength(readPtr, pos - params->nonNPgLen);
                    if (p)
                        PgHelpers::reverseComplementInPlace(readPtr, params->readLength);
                }
                res.append(read);
                res.push_back('\n');
            }
            fout << res;
            fout.close();
        }
    }

    template void PgRCDecoder::writeAllReadsInORDModeParallelChunks<uint_pg_len_std>(const string &outPrefix,
                                                                                         vector<uint_pg_len_std> &orgIdx2PgPos);

    template void PgRCDecoder::writeAllReadsInORDModeParallelChunks<uint_pg_len_max>(const string &outPrefix,
                                                                                         vector<uint_pg_len_max> &orgIdx2PgPos);

    template void PgRCDecoder::writeAllReadsInORDModeWithWritingThreads<uint_pg_len_std>(const string &outPrefix,
                                                                       vector<uint_pg_len_std> &orgIdx2PgPos);

    template void PgRCDecoder::writeAllReadsInORDModeWithWritingThreads<uint_pg_len_max>(const string &outPrefix,
                                                                       vector<uint_pg_len_max> &orgIdx2PgPos);

    template void PgRCDecoder::writeAllReadsInORDMode<uint_pg_len_std>(const string &outPrefix,
                                                                       vector<uint_pg_len_std> &orgIdx2PgPos) const;

    template void PgRCDecoder::writeAllReadsInORDMode<uint_pg_len_max>(const string &outPrefix,
                                                                       vector<uint_pg_len_max> &orgIdx2PgPos) const;

    uint_reads_cnt_max PgRCDecoder::dnaStreamSize() const {
        return (data.hqPg->getReadsSetProperties()->readsCount + data.lqPg->getReadsSetProperties()->readsCount
                + data.nPg->getReadsSetProperties()->readsCount) * (data.hqPg->getReadsSetProperties()->maxReadLength + 1);
    }

    void PgRCDecoder::preparePgsForValidation() const {
        if (params->preserveOrderMode) {
            data.lqPg->getReadsList()->orgIdx.clear();
            data.nPg->getReadsList()->orgIdx.clear();
            uint8_t parts = params->singleReadsMode ? 1 : 2;
            for (uint_reads_cnt_std i = 0; i < params->readsTotalCount; i++) {
                uint_pg_len_max pos = params->isJoinedPgLengthStd ? data.orgIdx2StdPgPos[i] : data.orgIdx2PgPos[i];
                uint_reads_cnt_std orgIdx =
                        i < params->readsTotalCount / parts ? i * parts : (i - params->readsTotalCount / parts) * parts + 1;
                if (pos < params->hqPgLen)
                    data.hqPg->getReadsList()->pos.push_back(pos);
                else if (pos < params->nonNPgLen) {
                    data.lqPg->getReadsList()->pos.push_back(pos - params->hqPgLen);
                    data.lqPg->getReadsList()->orgIdx.push_back(orgIdx);
                } else {
                    data.nPg->getReadsList()->pos.push_back(pos - params->nonNPgLen);
                    data.nPg->getReadsList()->orgIdx.push_back(orgIdx);
                }
            }
            data.hqPg->getReadsList()->enableConstantAccess(true);
        } else {
            data.hqPg->getReadsList()->enableConstantAccess(true);
            data.lqPg->getReadsList()->enableConstantAccess(true);
            if (data.nPg) data.nPg->getReadsList()->enableConstantAccess(true);
        }
    }

    void PgRCDecoder::validateAllPgs() {
        vector<uint_reads_cnt_max> orgIdx2rlIdx = getAllPgsOrgIdxs2RlIdx();
        ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                params->srcFastqFile, params->pairFastqFile);

        vector<bool> validated(params->readsTotalCount, false);
        uint_reads_cnt_max notValidatedCount = 0;
        uint_reads_cnt_max errorsCount = 0;
        for (uint_reads_cnt_max i = 0; i < params->readsTotalCount; i++) {
            if (!allReadsIterator->moveNext()) {
                cout << "The number of compressed reads is too big (" << (params->readsTotalCount - i) << " reads more)."
                     << endl;
                break;
            }
            string read;
            uint_reads_cnt_max idx = orgIdx2rlIdx[i];
            //uint_reads_cnt_max idx = (i % 2) * (readsTotalCount / 2) + (i / 2);
            if (validated[idx]) {
                notValidatedCount++;
                continue;
            }
            validated[idx] = true;
            if (idx < params->hqReadsCount)
                read = data.hqPg->getRead(idx);
            else {
                if (idx < params->nonNPgReadsCount)
                    read = data.lqPg->getRead(idx - params->hqReadsCount);
                else
                    read = data.nPg->getRead(idx - params->nonNPgReadsCount);
                if (!params->singleReadsMode && !params->ignorePairOrderInformation && (i % 2 == 1))
                    PgHelpers::reverseComplementInPlace(read);
            }
            if (read != allReadsIterator->getRead())
                errorsCount++;
        }
        uint_reads_cnt_max missingCount = 0;
        while (allReadsIterator->moveNext())
            missingCount++;
        if (missingCount)
            cout << "The number of compressed reads is too small (" << (missingCount) << " reads missing)." << endl;
        if (notValidatedCount)
            cout << notValidatedCount << " of compressed reads could not been properly validated." << endl;
        if (errorsCount)
            cout << "Found " << errorsCount << " errors in compressed reads." << endl;
        if (!missingCount && !notValidatedCount && !errorsCount)
            cout << "Validation successful!" << endl;

        delete (allReadsIterator);

    }

    const vector<uint_reads_cnt_max> PgRCDecoder::getAllPgsOrgIdxs2RlIdx() const {
        vector<uint_reads_cnt_max> orgIdx2rlIdx;
        orgIdx2rlIdx.resize(params->readsTotalCount);
        for (uint_reads_cnt_max i = 0; i < data.hqPg->getReadsSetProperties()->readsCount; i++)
            orgIdx2rlIdx[data.hqPg->getReadsList()->orgIdx[i]] = i;
        for (uint_reads_cnt_max i = 0; i < data.lqPg->getReadsSetProperties()->readsCount; i++)
            orgIdx2rlIdx[data.lqPg->getReadsList()->orgIdx[i]] = data.hqPg->getReadsSetProperties()->readsCount + i;
        for (uint_reads_cnt_max i = 0; i < params->nPgReadsCount; i++)
            orgIdx2rlIdx[data.nPg->getReadsList()->orgIdx[i]] = params->nonNPgReadsCount + i;
        return orgIdx2rlIdx;
    }

    const uint_reads_cnt_max PgRCDecoder::getAllPgsOrgIdx(uint_reads_cnt_max idx) const {
        if (idx < params->hqReadsCount)
            return data.hqPg->getReadsList()->orgIdx[idx];
        else if (idx < params->nonNPgReadsCount)
            return data.lqPg->getReadsList()->orgIdx[idx - params->hqReadsCount];
        else
            return data.nPg->getReadsList()->orgIdx[idx - params->nonNPgReadsCount];
    }

    void PgRCDecoder::validatePgsOrder() {
        if (!params->preserveOrderMode && params->singleReadsMode)
            return;

        vector<bool> validated(params->readsTotalCount, false);
        uint_reads_cnt_max notValidatedCount = 0;
        uint_reads_cnt_max errorsCount = 0;
        if (params->preserveOrderMode) {
            data.rlIdxOrder = getAllPgsOrgIdxs2RlIdx();
            for (uint_reads_cnt_std i = 0; i < params->readsTotalCount; i++) {
                uint_reads_cnt_std rlIdx = data.rlIdxOrder[i];
                if (validated[rlIdx]) {
                    notValidatedCount++;
                    continue;
                }
                if (i != getAllPgsOrgIdx(rlIdx))
                    errorsCount++;
            }
        } else {
            for (uint_reads_cnt_std p = 0; p < params->readsTotalCount / 2; p++) {
                uint_reads_cnt_std rlIdx = data.rlIdxOrder[p * 2];
                uint_reads_cnt_std rlPairIdx = data.rlIdxOrder[p * 2 + 1];
                if (validated[rlIdx]) notValidatedCount++;
                if (validated[rlPairIdx]) notValidatedCount++;
                if (!validated[rlIdx] && !validated[rlPairIdx]) {
                    validated[rlIdx] = true;
                    validated[rlPairIdx] = true;
                    uint_reads_cnt_std orgIdx = getAllPgsOrgIdx(rlIdx);
                    uint_reads_cnt_std orgPairIdx = getAllPgsOrgIdx(rlPairIdx);
                    uint_reads_cnt_std smallerIdx = orgIdx < orgPairIdx ? orgIdx : orgPairIdx;
                    uint_reads_cnt_std largerIdx = orgIdx >= orgPairIdx ? orgIdx : orgPairIdx;
                    if (largerIdx - smallerIdx != 1 || smallerIdx % 2)
                        errorsCount++;
                    else if (!params->ignorePairOrderInformation && smallerIdx != orgIdx)
                        errorsCount++;
                }
            }
        }
        if (notValidatedCount)
            cout << "Order of " << notValidatedCount << " compressed reads could not been properly validated." << endl;
        if (errorsCount)
            cout << "Found " << errorsCount << " errors in compressed reads order." << endl;
        if (!notValidatedCount && !errorsCount)
            cout << "Order validation successful!" << endl;
    }

    template<typename uint_pg_len>
    void PgRCDecoder::applyRevComplPairFileToPgs(vector<uint_pg_len> &orgIdx2PgPos) {
        if (params->preserveOrderMode) {
            uint_reads_cnt_std hqRlIdx = 0;
            const uint_reads_cnt_max pairsCount = params->readsTotalCount / 2;
            for (uint_reads_cnt_max i = 0; i < pairsCount; i++) {
                uint_pg_len pgPos = orgIdx2PgPos[i];
                if (pgPos < params->hqPgLen)
                    hqRlIdx++;
            }
            for (uint_reads_cnt_max i = pairsCount; i < params->readsTotalCount; i++) {
                uint_pg_len pgPos = orgIdx2PgPos[i];
                if (pgPos < params->hqPgLen) {
                    data.hqPg->getReadsList()->revComp[hqRlIdx] = !data.hqPg->getReadsList()->revComp[hqRlIdx];
                    hqRlIdx++;
                }
            }
        } else {
            for (uint_reads_cnt_max i = 1; i < params->readsTotalCount; i += 2) {
                uint_reads_cnt_std idx = data.rlIdxOrder[i];
                if (idx < params->hqReadsCount)
                    data.hqPg->getReadsList()->revComp[idx] = !data.hqPg->getReadsList()->revComp[idx];
            }
        }
    }

    template void PgRCDecoder::applyRevComplPairFileToPgs<uint_pg_len_std>(vector<uint_pg_len_std> &orgIdx2PgPos);

    template void PgRCDecoder::applyRevComplPairFileToPgs<uint_pg_len_max>(vector<uint_pg_len_max> &orgIdx2PgPos);

    void PgRCDecoder::loadAllPgs(istream &pgrcIn) {
        chrono::steady_clock::time_point start_t = chrono::steady_clock::now();

        istream* propsIn = &pgrcIn;
        string propsString;
        if (params->isVersionAtLeast(1, 3)) {
            char basesOrder[5];
            PgHelpers::readArray(pgrcIn, basesOrder, sizeof(basesOrder));
            PgHelpers::reorderSymAndVal(basesOrder);
            prefetchCompressedCollectiveParallel(pgrcIn);
            readCompressed(pgrcIn, propsString);
            propsIn = new istringstream(propsString);
        }
        PseudoGenomeHeader hqPgh(*propsIn);
        ReadsSetProperties hqRsProp(*propsIn);
        if (confirmTextReadMode(*propsIn)) {
            cout << "Reads list text mode unsupported during decompression." << endl;
            exit(EXIT_FAILURE);
        }
        ExtendedReadsListWithConstantAccessOption *hqCaeRl =
                ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgrcIn,
                                                                                               &hqPgh, &hqRsProp,
                                                                                               params->srcFastqFile.empty()
                                                                                               ? ""
                                                                                               : params->pgSeqFinalHqPrefix,
                                                                                               params);
        if (params->isVersionAtLeast(1, 3)) {
            delete(propsIn);
            readCompressed(pgrcIn, propsString);
            propsIn = new istringstream(propsString);
        }
        PseudoGenomeHeader lqPgh(*propsIn);
        ReadsSetProperties lqRsProp(*propsIn);
        if (confirmTextReadMode(*propsIn)) {
            cout << "Reads list text mode unsupported during decompression." << endl;
            exit(EXIT_FAILURE);
        }
        ExtendedReadsListWithConstantAccessOption *lqCaeRl =
                ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgrcIn,
                                                                                               &lqPgh, &lqRsProp,
                                                                                               params->srcFastqFile.empty()
                                                                                               ? ""
                                                                                               : params->pgSeqFinalLqPrefix,
                                                                                               params, true, true);
        ExtendedReadsListWithConstantAccessOption *nCaeRl = nullptr;
        ReadsSetProperties nRsProp;
        PseudoGenomeHeader nPgh;
        if (params->separateNReads) {
            if (params->isVersionAtLeast(1, 3)) {
                delete(propsIn);
                readCompressed(pgrcIn, propsString);
                propsIn = new istringstream(propsString);
            }
            nPgh = PseudoGenomeHeader(*propsIn);
            nRsProp = ReadsSetProperties(*propsIn);
            if (confirmTextReadMode(*propsIn)) {
                cout << "Reads list text mode unsupported during decompression." << endl;
                exit(EXIT_FAILURE);
            }
            nCaeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgrcIn,
                                                                                                    &nPgh, &nRsProp,
                                                                                                    params->srcFastqFile.empty()
                                                                                                    ? ""
                                                                                                    : params->pgNPrefix,
                                                                                                    params, true, true);
        }
        if (params->isVersionAtLeast(1, 3))
            delete(propsIn);
        params->readLength = hqRsProp.maxReadLength;
        params->hqReadsCount = hqRsProp.readsCount;
        params->lqReadsCount = lqRsProp.readsCount;
        params->nonNPgReadsCount = params->hqReadsCount + params->lqReadsCount;
        params->nPgReadsCount = params->separateNReads ? nRsProp.readsCount : 0;
        params->readsTotalCount = params->nonNPgReadsCount + params->nPgReadsCount;

        params->hqPgLen = hqPgh.getPseudoGenomeLength();
        params->nonNPgLen = params->hqPgLen + lqPgh.getPseudoGenomeLength();
        if (params->preserveOrderMode) {
            params->isJoinedPgLengthStd = params->nonNPgLen + nPgh.getPseudoGenomeLength() <= UINT32_MAX;
            if (params->isJoinedPgLengthStd)
                SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_std>(pgrcIn,
                                                                                              data.orgIdx2StdPgPos,
                                                                                              params);
            else
                SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_max>(pgrcIn,
                                                                                              data.orgIdx2PgPos,
                                                                                              params);
        } else {
            SeparatedPseudoGenomePersistence::decompressReadsOrder(pgrcIn, data.rlIdxOrder,
                                                                   params->preserveOrderMode,
                                                                   params->ignorePairOrderInformation,
                                                                   params->singleReadsMode);
        }
        cout << "... loaded Pgs Reads Lists (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        string hqPgSeq, lqPgSeq, nPgSeq;
        SimplePgMatcher::restoreMatchedPgs(pgrcIn, params->hqPgLen, hqPgSeq, lqPgSeq, nPgSeq,
                                           params);
        data.hqPg = new SeparatedPseudoGenome(move(hqPgSeq), hqCaeRl, &hqRsProp);
        data.lqPg = new SeparatedPseudoGenome(move(lqPgSeq), lqCaeRl, &lqRsProp);
        data.nPg = new SeparatedPseudoGenome(move(nPgSeq), nCaeRl, &nRsProp);
    }

    void PgRCDecoder::loadAllPgs() {
        string hqPgSeq = SimplePgMatcher::restoreAutoMatchedPg(params->pgSeqFinalHqPrefix, true);
        PseudoGenomeHeader *pgh = nullptr;
        ReadsSetProperties *prop = nullptr;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(params->pgSeqFinalHqPrefix, pgh, prop, plainTextReadMode);
        ExtendedReadsListWithConstantAccessOption *caeRl =
                ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(params->pgSeqFinalHqPrefix,
                                                                                               pgh->getPseudoGenomeLength());
        data.hqPg = new SeparatedPseudoGenome(move(hqPgSeq), caeRl, prop);
        delete (pgh);
        delete (prop);

        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(params->pgSeqFinalLqPrefix, pgh, prop, plainTextReadMode);
        string lqPgSeq = SimplePgMatcher::restoreMatchedPg(data.hqPg->getPgSequence(), params->pgSeqFinalLqPrefix,
                                                           true, plainTextReadMode);
        caeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(params->pgSeqFinalLqPrefix,
                                                                                               pgh->getPseudoGenomeLength());
        data.lqPg = new SeparatedPseudoGenome(move(lqPgSeq), caeRl, prop);
        delete (pgh);
        delete (prop);

        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(params->pgNPrefix, pgh, prop, plainTextReadMode);
        string nPgSeq = SimplePgMatcher::restoreMatchedPg(data.hqPg->getPgSequence(), params->pgNPrefix, true,
                                                          plainTextReadMode);
        caeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(params->pgNPrefix,
                                                                                               pgh->getPseudoGenomeLength());
        data.nPg = new SeparatedPseudoGenome(move(nPgSeq), caeRl, prop);
        delete (pgh);
        delete (prop);

        params->readLength = data.hqPg->getReadsSetProperties()->maxReadLength;
    }

}