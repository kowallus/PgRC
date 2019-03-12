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

    class PgRCManager {
    private:
        static const int MIN_CHARS_PER_PGMATCH = 20;
        static const int MIN_CHARS_PER_MISMATCH = 4;
        static const int MIN_READS_EXACT_MATCHING_CHARS = 20;

        // INPUT PARAMETERS
        uint8_t compressionLevel = PGRC_CODER_LEVEL_NORMAL;
        bool singleReadsMode = false;
        bool preserveOrderMode = false;
        uint16_t error_limit_in_promils = 1000;
        string gen_quality_str = "50";
        double gen_quality_coef = 0.5;
        bool nReadsLQ = false;
        bool separateNReads = false;
        bool extraFilesForValidation = false;
        uint16_t preReadsExactMatchingChars = 0;
        uint16_t readsExactMatchingChars = UINT16_MAX;
        uint16_t minCharsPerMismatch = UINT16_MAX;
        char preMatchingMode = 'd';
        char matchingMode = 'd';
        uint32_t targetPgMatchLength = 50;
        string pgRCFileName = "";
        bool revComplPairFile = false;

        string srcFastqFile = "";
        string pairFastqFile = "";

        // CHAIN MANAGEMENT
        uint8_t skipStages = 0;
        uint8_t endAtStage = UINT8_MAX;

        // TESTING
        bool disableInMemoryMode = false;
        clock_t start_t;
        clock_t div_t;
        clock_t pgDiv_t;
        clock_t good_t;
        clock_t match_t;
        clock_t bad_t;
        clock_t gooder_t;

        void reportTimes();

        // CHAIN VARIABLES
        uint_read_len_max readLength;
        uint8_t stageCount;
        fstream pgrcOut;

        DividedPCLReadsSets *divReadsSets = 0;
        SeparatedPseudoGenome *hqPg = 0;
        SeparatedPseudoGenome *lqPg = 0;
        SeparatedPseudoGenome *nPg = 0;

        vector<uint_reads_cnt_std> orgIdxs;

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

        // CHAIN METHODS
        void prepareChainData();

    public:

        PgRCManager() {
            setError_limit_in_promils(1000);
            setGen_quality_str("50");
        }

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

        void setSingleReadsMode() {
            if (preserveOrderMode) {
                fprintf(stderr, "Single reads and preserve order modes cannot be used together.");
                exit(EXIT_FAILURE);
            }
            PgRCManager::singleReadsMode = true;
        }

        void setError_limit_in_promils(uint16_t error_limit_in_promils) {
            if (error_limit_in_promils > 1000) {
                fprintf(stderr, "Error limit should not be greater than 1000.\n");
                exit(EXIT_FAILURE);
            }
            PgRCManager::error_limit_in_promils = error_limit_in_promils;
        }

        void setGen_quality_str(const string &gen_quality_str) {
            gen_quality_coef = atoi(gen_quality_str.c_str()) / 100.0;
            if (gen_quality_coef > 1 || gen_quality_coef <= 0) {
                fprintf(stderr, "Generate quality coefficient should be between 1 and 100.\n");
                exit(EXIT_FAILURE);
            }
            PgRCManager::gen_quality_str = gen_quality_str;
        }

        void setNReadsLQ() {
            if (separateNReads) {
                fprintf(stderr, "Reads containing N should be processed either separately or consider low quality");
                exit(EXIT_FAILURE);
            }
            PgRCManager::nReadsLQ = true;
        }

        void setSeparateNReads() {
            if (nReadsLQ) {
                fprintf(stderr, "Reads containing N should be processed either separately or considered low quality");
                exit(EXIT_FAILURE);
            }
            PgRCManager::separateNReads = true;
        }

        void setValidationOutputMode() {
            PgRCManager::extraFilesForValidation = true;
        }

        void setReadsExactMatchingChars(uint16_t readsExactMatchingChars) {
            if (readsExactMatchingChars < MIN_READS_EXACT_MATCHING_CHARS) {
                fprintf(stderr, "Chars per reads exact matching cannot be lower than %d.\n",
                        MIN_READS_EXACT_MATCHING_CHARS);
                exit(EXIT_FAILURE);
            }
            PgRCManager::readsExactMatchingChars = readsExactMatchingChars;
        }

        void setPreReadsExactMatchingChars(uint16_t preReadsExactMatchingChars) {
            if (preReadsExactMatchingChars < MIN_READS_EXACT_MATCHING_CHARS) {
                fprintf(stderr, "Chars per reads exact matching cannot be lower than %d.\n",
                        MIN_READS_EXACT_MATCHING_CHARS);
                exit(EXIT_FAILURE);
            }
            PgRCManager::preReadsExactMatchingChars = preReadsExactMatchingChars;
        }

        void setMaxCharsPerMismatch(uint16_t minCharsPerMismatch) {
            if (minCharsPerMismatch < MIN_CHARS_PER_MISMATCH) {
                fprintf(stderr, "Chars per mismatch cannot be lower than %d.\n", MIN_CHARS_PER_MISMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCManager::minCharsPerMismatch = minCharsPerMismatch;
        }

        void setPreMatchingMode(char matchingModde) {
            PgRCManager::preMatchingMode = matchingModde;
        }

        void setMatchingMode(char matchingMode) {
            PgRCManager::matchingMode = matchingMode;
        }

        void setMinimalPgMatchLength(uint32_t targetPgMatchLength) {
            if (targetPgMatchLength < MIN_CHARS_PER_PGMATCH) {
                fprintf(stderr, "Target Pg match length cannot be lower than %d.\n", MIN_CHARS_PER_PGMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCManager::targetPgMatchLength = targetPgMatchLength;
        }

        void setPgRCFileName(const string &pgRCFileName) {
            PgRCManager::pgRCFileName = pgRCFileName;
        }

        void setRevComplPairFile() {
            PgRCManager::revComplPairFile = true;
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

        void loadAllPgs(istream &pgrcIn, vector<uint_reads_cnt_std>& rlIdxOrder,
                bool completeOrderInfo, bool singleFileMode);

        void loadAllPgs();

        void decompressPgRC();

        const size_t CHUNK_SIZE_IN_BYTES = 100000;

        void writeAllReadsInSEMode(const string &tmpDirectoryPath) const;
        void writeAllReads(const string &tmpDirectoryPath, vector<uint_reads_cnt_std> &rlIdxOrder, bool singleFileMode) const;

        std::mutex mut;
        std::queue<string> out_queue;
        std::condition_variable data_cond;

        void writeAllReadsInSEModeParallel(const string &tmpDirectoryPath);

        void pushOutToQueue(string &out);

        void finishWritingParallel();

        void writeFromQueue(const string &tmpDirectoryPath);

        void validateAllPgs();

        uint_reads_cnt_max dnaStreamSize() const;
    };
}


#endif //PGTOOLS_PGRCMANAGER_H
