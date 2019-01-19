#ifndef PGTOOLS_PGRCMANAGER_H
#define PGTOOLS_PGRCMANAGER_H

#include "utils/helper.h"
#include "pgsaconfig.h"

#include "readsset/DividedPCLReadsSets.h"

namespace PgTools {

    class PgRCManager {
    private:
        static const int MIN_CHARS_PER_PGMATCH = 20;
        static const int MIN_CHARS_PER_MISMATCH = 4;

        // INPUT PARAMETERS
        uint16_t error_limit_in_promils = 1000;
        string gen_quality_str = "50";
        double gen_quality_coef = 0.5;
        bool nReadsLQ = false;
        bool separateNReads = false;
        uint16_t targetCharsPerMismatch = UINT16_MAX;
        uint16_t maxCharsPerMismatch = UINT16_MAX;
        char mismatchesMode = 'd';
        uint32_t minimalPgMatchLength = 50;
        string pgFilesPrefixes = "";
        bool revComplPairFile = false;

        string srcFastqFile = "";
        string pairFastqFile = "";

        // CHAIN MANAGEMENT
        bool skipIntermediateOutput = true;
        uint8_t skipStages = 0;
        uint8_t endAtStage = UINT8_MAX;

        // TESTING
        bool disableInMemoryMode = true;
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
        uint8_t targetMismatches;
        uint8_t maxMismatches;
        uint8_t stageCount;

        DividedPCLReadsSets* divReadsSets = 0;

        bool qualityDivision;
        string lqDivisionFile;
        string nDivisionFile;
        string pgGoodPrefix;
        string pgFilesPrefixesWithM;
        string pgMappedGoodPrefix;
        string pgMappedBadPrefix;
        string pgNPrefix;
        string mappedBadDivisionFile;

        // CHAIN METHODS
        void prepareChainData();

    public:

        PgRCManager() {
            setError_limit_in_promils(1000);
            setGen_quality_str("50");
        }

        void executePgRCChain();

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

        void setNReadsLQ(bool nReadsLQ) {
            if (separateNReads) {
                fprintf(stderr, "Reads containing N should be processed either separately or consider low quality");
                exit(EXIT_FAILURE);
            }
            PgRCManager::nReadsLQ = nReadsLQ;
        }

        void setSeparateNReads(bool separateNReads) {
            if (nReadsLQ) {
                fprintf(stderr, "Reads containing N should be processed either separately or consider low quality");
                exit(EXIT_FAILURE);
            }
            PgRCManager::separateNReads = separateNReads;
        }

        void setTargetCharsPerMismatch(uint16_t targetCharsPerMismatch) {
            if (targetCharsPerMismatch < MIN_CHARS_PER_MISMATCH || maxCharsPerMismatch < MIN_CHARS_PER_MISMATCH) {
                fprintf(stderr, "Chars per mismatch cannot be lower than %d.\n", MIN_CHARS_PER_MISMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCManager::targetCharsPerMismatch = targetCharsPerMismatch;
        }

        void setMaxCharsPerMismatch(uint16_t maxCharsPerMismatch) {
            if (maxCharsPerMismatch > targetCharsPerMismatch) {
                fprintf(stdout, "allowedMaxMismatches cannot be smaller than targetMaxMismatches.\n");
                exit(EXIT_FAILURE);
            }
            PgRCManager::maxCharsPerMismatch = maxCharsPerMismatch;
        }

        void setMismatchesMode(char mismatchesMode) {
            PgRCManager::mismatchesMode = mismatchesMode;
        }

        void setMinimalPgMatchLength(uint32_t minimalPgMatchLength) {
            if (minimalPgMatchLength < MIN_CHARS_PER_PGMATCH) {
                fprintf(stderr, "Minimal Pg match length cannot be lower than %d.\n", MIN_CHARS_PER_PGMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCManager::minimalPgMatchLength = minimalPgMatchLength;
        }

        void setPgFilesPrefixes(const string &pgFilesPrefixes) {
            PgRCManager::pgFilesPrefixes = pgFilesPrefixes;
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

        void generateIntermediateOutput() {
            PgRCManager::skipIntermediateOutput = false;
        }

        void setSkipStages(uint8_t skipStages) {
            if (skipStages >= endAtStage) {
                fprintf(stdout, "Number of stages to skip (%d) should be smaller than a number of a stage to finish (%d).\n",
                        skipStages, endAtStage);
                exit(EXIT_FAILURE);
            }
            PgRCManager::skipStages = skipStages;
        }

        void setEndAtStage(uint8_t endAtStage) {
            if (skipStages >= endAtStage) {
                fprintf(stdout, "Number of stages to skip (%d) should be smaller than a number of a stage to finish (%d).\n",
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

        void persistQualityBasedDivision();

        void persistPgGeneratorBasedReadsDivision();

        void persistHQPg();

        void saveHQPgReadsList();

        void saveLQPg();

        void saveLQPgReadsList();

        void saveHQPgSequence();

        void saveNPg();

        void saveMEMMappedPgSequences();

        void extractHQPgSequence();

        void freeHQPg();

        void extractLQPgSequence();

        void freeLQPg();

        void prepareForPgGeneratorBaseReadsDivision();

        void disposeChainData();

        void prepareForHqPgGeneration();
    };
}


#endif //PGTOOLS_PGRCMANAGER_H
