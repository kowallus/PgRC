#ifndef PGTOOLS_PGRCMANAGER_H
#define PGTOOLS_PGRCMANAGER_H

#include "utils/helper.h"
#include "pgsaconfig.h"
#include "utils/LzmaLib.h"

#include "readsset/DividedPCLReadsSets.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

#include <condition_variable>
#include <queue>
#include <mutex>
#include <thread>

#define ENABLE_PARALLEL_DECOMPRESSION false

namespace PgTools {

    static const char PGRC_SE_MODE = 0;
    static const char PGRC_PE_MODE = 1;
    static const char PGRC_ORD_SE_MODE = 2;
    static const char PGRC_ORD_PE_MODE = 3;
    static const char PGRC_MIN_PE_MODE = 4;

    class PgRCManager {
    private:
        static const int MIN_CHARS_PER_PGMATCH = 20;
        static const int MIN_CHARS_PER_MISMATCH = 4;
        static const int MIN_READS_EXACT_MATCHING_CHARS = 20;

        static const char DEFAULT_CHAR_PARAM = CHAR_MAX;
        static const uint16_t DEFAULT_UINT16_PARAM = UINT16_MAX;
        static constexpr double DEFAULT_DOUBLE_PARAM = -1;

        // INPUT PARAMETERS
        bool singleReadsMode = false;
        bool preserveOrderMode = false;
        bool ignorePairOrderInformation = false;
        bool nReadsLQ = false;
        bool separateNReads = true;
        bool extraFilesForValidation = false;
        string pgRCFileName = "";
        bool disableRevComplPairFileMode = false;

        string srcFastqFile = "";
        string pairFastqFile = "";

        // COMPRESSION PARAMETERS
        uint8_t compressionLevel = PGRC_CODER_LEVEL_NORMAL;
        bool forceConstantParamsMode = true;
        uint16_t error_limit_in_promils = DEFAULT_UINT16_PARAM;
        string gen_quality_str;
        double gen_quality_coef = DEFAULT_DOUBLE_PARAM;
        uint16_t preReadsExactMatchingChars = DEFAULT_UINT16_PARAM;
        uint16_t readsExactMatchingChars = DEFAULT_UINT16_PARAM;
        uint16_t minCharsPerMismatch = DEFAULT_UINT16_PARAM;
        char preMatchingMode = CHAR_MAX;
        char matchingMode = CHAR_MAX;
        uint16_t targetPgMatchLength = DEFAULT_UINT16_PARAM;

        // CHAIN MANAGEMENT
        uint8_t skipStages = 0;
        uint8_t endAtStage = UINT8_MAX;

        // TESTING&REPORTING
        bool disableInMemoryMode = false;
        clock_t start_t;
        clock_t div_t;
        clock_t pgDiv_t;
        clock_t good_t;
        clock_t match_t;
        clock_t bad_t;
        clock_t pgSeqs_t;
        size_t pgRCSize = 0;

        void generateReport();

        // CHAIN VARIABLES
        uint_read_len_max readLength;
        uint8_t stageCount;
        fstream pgrcOut;

        DividedPCLReadsSets *divReadsSets = 0;
        SeparatedPseudoGenome *hqPg = 0;
        SeparatedPseudoGenome *lqPg = 0;
        SeparatedPseudoGenome *nPg = 0;

        vector<uint_reads_cnt_std> orgIdxs;

        bool revComplPairFile;
        bool qualityDivision;
        string lqDivisionFile;
        string nDivisionFile;
        string pgHqPrefix;
        string pgFilesPrefixesWithM;
        string pgMappedHqPrefix;
        string pgMappedLqPrefix;
        string pgSeqFinalHqPrefix;
        string pgSeqFinalLqPrefix;
        string pgNPrefix;
        string mappedLqDivisionFile;

        // Decompression chain variables
        uint_reads_cnt_max hqReadsCount;
        uint_reads_cnt_max lqReadsCount;
        uint_reads_cnt_max nonNPgReadsCount;
        uint_reads_cnt_max nPgReadsCount;
        uint_reads_cnt_max readsTotalCount;

        // CHAIN METHODS
        void prepareChainData();

        void initCompressionParameters();

    public:

        PgRCManager() { }

        void executePgRCChain();

        void setCompressionLevel(int compressionLevel) {
            if (compressionLevel > PGRC_CODER_LEVEL_MAX || compressionLevel < PGRC_CODER_LEVEL_FAST) {
                fprintf(stderr, "Generate quality coefficient should be between %d and %d.\n",
                        PGRC_CODER_LEVEL_FAST, PGRC_CODER_LEVEL_MAX);
                exit(EXIT_FAILURE);
            }
            PgRCManager::compressionLevel = compressionLevel;
        }

        void setPreserveOrderMode() {
            if (singleReadsMode) {
                fprintf(stderr, "Single reads and preserve order modes cannot be used together.");
                exit(EXIT_FAILURE);
            }
            PgRCManager::preserveOrderMode = true;
        }

        void setIgnorePairOrderInformation() {
            if (singleReadsMode || preserveOrderMode) {
                fprintf(stderr, "Ignore pair order works only with default PE mode.");
                exit(EXIT_FAILURE);
            }
            PgRCManager::ignorePairOrderInformation = true;
        };

        void setSingleReadsMode() {
            if (preserveOrderMode || ignorePairOrderInformation) {
                fprintf(stderr, "Single reads and ordering parameters cannot be used together.");
                exit(EXIT_FAILURE);
            }
            PgRCManager::singleReadsMode = true;
        }

        void allowVariableParams() {
            PgRCManager::forceConstantParamsMode = false;
        }

        void setError_limit_in_promils(uint16_t error_limit_in_promils) {
            if (PgRCManager::error_limit_in_promils != DEFAULT_UINT16_PARAM)
                return;
            if (error_limit_in_promils > 1000) {
                fprintf(stderr, "Error limit should not be greater than 1000.\n");
                exit(EXIT_FAILURE);
            }
            PgRCManager::error_limit_in_promils = error_limit_in_promils;
        }

        void setGen_quality_str(const string &gen_quality_str) {
            if (PgRCManager::gen_quality_coef != DEFAULT_DOUBLE_PARAM)
                return;
            gen_quality_coef = atoi(gen_quality_str.c_str()) / 100.0;
            if (gen_quality_coef > 1 || gen_quality_coef <= 0) {
                fprintf(stderr, "Generate quality coefficient should be between 1 and 100.\n");
                exit(EXIT_FAILURE);
            }
            PgRCManager::gen_quality_str = gen_quality_str;
        }

        void setNReadsLQ() {
            PgRCManager::separateNReads = false;
            PgRCManager::nReadsLQ = true;
        }

        void doNotSeparateNReads() {
            PgRCManager::separateNReads = false;
        }

        void setValidationOutputMode() {
            PgRCManager::extraFilesForValidation = true;
        }

        void setReadsExactMatchingChars(uint16_t readsExactMatchingChars) {
            if (PgRCManager::readsExactMatchingChars != DEFAULT_UINT16_PARAM)
                return;
            if (readsExactMatchingChars < MIN_READS_EXACT_MATCHING_CHARS) {
                fprintf(stderr, "Chars per reads exact matching cannot be lower than %d.\n",
                        MIN_READS_EXACT_MATCHING_CHARS);
                exit(EXIT_FAILURE);
            }
            PgRCManager::readsExactMatchingChars = readsExactMatchingChars;
        }

        void setPreReadsExactMatchingChars(uint16_t preReadsExactMatchingChars) {
            if (PgRCManager::preReadsExactMatchingChars != DEFAULT_UINT16_PARAM)
                return;
            if (preReadsExactMatchingChars < MIN_READS_EXACT_MATCHING_CHARS &&
                preReadsExactMatchingChars > 0) {
                fprintf(stderr, "Chars per reads exact matching cannot be lower than %d.\n",
                        MIN_READS_EXACT_MATCHING_CHARS);
                exit(EXIT_FAILURE);
            }
            PgRCManager::preReadsExactMatchingChars = preReadsExactMatchingChars;
        }

        void setMinCharsPerMismatch(uint16_t minCharsPerMismatch) {
            if (PgRCManager::minCharsPerMismatch != DEFAULT_UINT16_PARAM)
                return;
            if (minCharsPerMismatch < MIN_CHARS_PER_MISMATCH) {
                fprintf(stderr, "Chars per mismatch cannot be lower than %d.\n", MIN_CHARS_PER_MISMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCManager::minCharsPerMismatch = minCharsPerMismatch;
        }

        void setPreMatchingMode(char matchingMode) {
            if (PgRCManager::preMatchingMode != DEFAULT_CHAR_PARAM)
                return;
            PgRCManager::preMatchingMode = matchingMode;
        }

        void setMatchingMode(char matchingMode) {
            if (PgRCManager::matchingMode != DEFAULT_CHAR_PARAM)
                return;
            PgRCManager::matchingMode = matchingMode;
        }

        void setMinimalPgMatchLength(uint16_t targetPgMatchLength) {
            if (PgRCManager::targetPgMatchLength != DEFAULT_UINT16_PARAM)
                return;
            if (targetPgMatchLength < MIN_CHARS_PER_PGMATCH) {
                fprintf(stderr, "Target Pg match length cannot be lower than %d.\n", MIN_CHARS_PER_PGMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCManager::targetPgMatchLength = targetPgMatchLength;
        }

        void setPgRCFileName(const string &pgRCFileName) {
            PgRCManager::pgRCFileName = pgRCFileName;
        }

        void disableRevComplPairFile() {
            PgRCManager::disableRevComplPairFileMode = true;
        }

        void setSrcFastqFile(const string &srcFastqFile) {
            PgRCManager::srcFastqFile = srcFastqFile;
        }

        void setPairFastqFile(const string &pairFastqFile) {
            PgRCManager::pairFastqFile = pairFastqFile;
        }

        void setSkipStages(uint8_t skipStages) {
            if (skipStages >= endAtStage) {
                fprintf(stdout,
                        "Number of stages to skip (%d) should be smaller than a number of a stage to finish (%d).\n",
                        skipStages, endAtStage);
                exit(EXIT_FAILURE);
            }
            PgRCManager::skipStages = skipStages;
        }

        void setEndAtStage(uint8_t endAtStage) {
            if (skipStages >= endAtStage) {
                fprintf(stdout,
                        "Number of stages to skip (%d) should be smaller than a number of a stage to finish (%d).\n",
                        skipStages, endAtStage);
                exit(EXIT_FAILURE);
            }
            PgRCManager::endAtStage = endAtStage;
        }

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

        void compressMEMMappedPgSequences();

        void prepareForPgGeneratorBaseReadsDivision();

        void disposeChainData();

        void prepareForHqPgGeneration();

        void prepareForMappingLQReadsOnHQPg();

        void prepareForLQPgAndNPgGeneration();

        void persistMappedReadsQualityDivision();

        void prepareForPgMatching();

        void finalizeCompression();

        void loadAllPgs(istream &pgrcIn, vector<uint_reads_cnt_std>& rlIdxOrder);
        void loadAllPgs();

        void decompressPgRC();

        void applyRevComplPairFileToPgs(vector<uint_reads_cnt_std>& rlIdxOrder);

        const size_t CHUNK_SIZE_IN_BYTES = 100000;

        void writeAllReadsInSEMode(const string &outPrefix) const;
        void writeAllReads(const string &outPrefix, vector<uint_reads_cnt_std> &rlIdxOrder) const;

        std::mutex mut;
        std::queue<string> out_queue;
        std::condition_variable data_cond;

        void writeAllReadsInSEModeParallel(const string &outPrefix);

        void pushOutToQueue(string &out);

        void finishWritingParallel();

        void writeFromQueue(const string &outPrefix);

        void validateAllPgs();
        void validatePgsOrder(vector<uint_reads_cnt_std>& rlIdxOrder);

        uint_reads_cnt_max dnaStreamSize() const;

        const vector<uint_reads_cnt_max> getAllPgsOrgIdxs2RlIdx() const;
        const uint_reads_cnt_max getAllPgsOrgIdx(uint_reads_cnt_max idx) const;
    };
}


#endif //PGTOOLS_PGRCMANAGER_H
