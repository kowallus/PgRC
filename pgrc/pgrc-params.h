#ifndef PGTOOLS_PGRCPARAMS_H
#define PGTOOLS_PGRCPARAMS_H

#include "../coders/CodersLib.h"
#include "pg-config.h"

namespace PgTools {

    using namespace PgIndex;

    static const char PGRC_SE_MODE = 0;
    static const char PGRC_PE_MODE = 1;
    static const char PGRC_ORD_SE_MODE = 2;
    static const char PGRC_ORD_PE_MODE = 3;
    static const char PGRC_MIN_PE_MODE = 4;

    static const char PGRC_VERSION_MODE = '#';
    static const char PGRC_VERSION_MAJOR = 2;
    static const char PGRC_VERSION_MINOR = 0;
    static const char PGRC_VERSION_REVISION = 1;

    static const char *const BAD_INFIX = "bad";
    static const char *const GOOD_INFIX = "good";
    static const char *const N_INFIX = "N";
    static const char *const DIVISION_EXTENSION = ".div";
    static const char *const TEMPORARY_FILE_SUFFIX = ".temp";
    static const char *const PGRC_HEADER = "PgRC";

    static const int MIN_CHARS_PER_PGMATCH = 20;
    static const int MIN_CHARS_PER_MISMATCH = 2;
    static const int MIN_READS_EXACT_MATCHING_CHARS = 20;

    static const char DEFAULT_CHAR_PARAM = CHAR_MAX;
    static const uint16_t DEFAULT_UINT16_PARAM = UINT16_MAX;
    static constexpr double DEFAULT_DOUBLE_PARAM = -1;

    class PgRCParams {
    public:

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
        uint8_t compressionLevel = CODER_LEVEL_NORMAL;
        bool forceConstantParamsMode = true;
        uint16_t error_limit_in_promils = DEFAULT_UINT16_PARAM;
        bool simplified_suffix_mode = true;
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

        // CHAIN PARAMS
        uint_read_len_max readLength;

        bool revComplPairFile;
        bool qualityDivision;
        bool generatorDivision;
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
        char pgrcVersionMajor = 1;
        char pgrcVersionMinor = 0;
        char pgrcVersionRevision = 0;

        bool isVersion(char major, char minor) {
            return pgrcVersionMajor == major && pgrcVersionMinor == minor;
        }

        bool isVersionAtLeast(char major, char minor) {
            return (pgrcVersionMajor > major ||
                    (pgrcVersionMajor == major && pgrcVersionMinor >= minor));
        }

        // CHAIN VARIABLES
        uint_reads_cnt_max hqReadsCount;
        uint_reads_cnt_max lqReadsCount;
        uint_reads_cnt_max nonNPgReadsCount;
        uint_reads_cnt_max nPgReadsCount;
        uint_reads_cnt_max readsTotalCount;
        uint_pg_len_max hqPgLen;
        uint_pg_len_max nonNPgLen;
        bool isJoinedPgLengthStd;

        void initCompressionParameters() {
            setPreMatchingMode('c');
            switch (compressionLevel) {
                case CODER_LEVEL_FAST:
/*                    setQualityBasedDivisionErrorLimitInPromils(1);
                    disableSimplifiedSuffixMode4QualityBasedDivision();
                    setPgGeneratorBasedDivisionOverlapThreshold_str("0");
                    setMatchingMode('C');
                    setPreReadsExactMatchingChars(0);
                    setReadSeedLength(54);
                    setMinCharsPerMismatch(2);
                    setMinimalPgReverseComplementedRepeatLength(25);
                    break;*/
                case CODER_LEVEL_MAX:
/*                    setQualityBasedDivisionErrorLimitInPromils(200);
                    disableSimplifiedSuffixMode4QualityBasedDivision();
                    setPgGeneratorBasedDivisionOverlapThreshold_str("65");
                    setPreReadsExactMatchingChars(64);
                    setMatchingMode('c');
                    setReadSeedLength(30);
                    setMinCharsPerMismatch(4);
                    setMinimalPgReverseComplementedRepeatLength(50);
                    break;*/
                case CODER_LEVEL_NORMAL:
                    setQualityBasedDivisionErrorLimitInPromils(120);
                    setPgGeneratorBasedDivisionOverlapThreshold_str("65");
                    setPreReadsExactMatchingChars(0);
                    setMatchingMode('c');
                    setReadSeedLength(38);
                    setMinCharsPerMismatch(3);
                    setMinimalPgReverseComplementedRepeatLength(45);
                    break;
                default:
                    fprintf(stderr, "Error: unknown compression level: %d.", compressionLevel);
                    exit(EXIT_FAILURE);
            }
        }

        void setCompressionLevel(int compressionLevel) {
            if (compressionLevel > CODER_LEVEL_MAX || compressionLevel < CODER_LEVEL_FAST) {
                fprintf(stderr, "Generate quality coefficient should be between %d and %d.\n",
                        CODER_LEVEL_FAST, CODER_LEVEL_MAX);
                exit(EXIT_FAILURE);
            }
            PgRCParams::compressionLevel = compressionLevel;
        }

        void setPreserveOrderMode() {
            if (singleReadsMode) {
                fprintf(stderr, "Single reads and preserve order modes cannot be used together.");
                exit(EXIT_FAILURE);
            }
            PgRCParams::preserveOrderMode = true;
        }

        void setIgnorePairOrderInformation() {
            if (singleReadsMode || preserveOrderMode) {
                fprintf(stderr, "Ignore pair order works only with default PE mode.");
                exit(EXIT_FAILURE);
            }
            PgRCParams::ignorePairOrderInformation = true;
        };

        void setSingleReadsMode() {
            if (preserveOrderMode || ignorePairOrderInformation) {
                fprintf(stderr, "Single reads and ordering parameters cannot be used together.");
                exit(EXIT_FAILURE);
            }
            PgRCParams::singleReadsMode = true;
        }

        void allowVariableParams() {
            PgRCParams::forceConstantParamsMode = false;
        }

        void setQualityBasedDivisionErrorLimitInPromils(uint16_t error_limit_in_promils) {
            if (PgRCParams::error_limit_in_promils != DEFAULT_UINT16_PARAM)
                return;
            if (error_limit_in_promils > 1000) {
                fprintf(stderr, "Error limit should not be greater than 1000.\n");
                exit(EXIT_FAILURE);
            }
            PgRCParams::error_limit_in_promils = error_limit_in_promils;
        }

        void disableSimplifiedSuffixMode4QualityBasedDivision() {
            PgRCParams::simplified_suffix_mode = false;
        }

        void setPgGeneratorBasedDivisionOverlapThreshold_str(const string &gen_quality_str) {
            if (PgRCParams::gen_quality_coef != DEFAULT_DOUBLE_PARAM)
                return;
            gen_quality_coef = atoi(gen_quality_str.c_str()) / 100.0;
            if (gen_quality_coef > 1 || gen_quality_coef < 0) {
                fprintf(stderr, "Generate quality coefficient should be between 0 and 100.\n");
                exit(EXIT_FAILURE);
            }
            PgRCParams::gen_quality_str = gen_quality_str;
        }

        void setNReadsLQ() {
            PgRCParams::separateNReads = false;
            PgRCParams::nReadsLQ = true;
        }

        void doNotSeparateNReads() {
            PgRCParams::separateNReads = false;
        }

        void setValidationOutputMode() {
            PgRCParams::extraFilesForValidation = true;
        }

        void setReadSeedLength(uint16_t readsExactMatchingChars) {
            if (PgRCParams::readsExactMatchingChars != DEFAULT_UINT16_PARAM)
                return;
            if (readsExactMatchingChars < MIN_READS_EXACT_MATCHING_CHARS) {
                fprintf(stderr, "Chars per reads exact matching cannot be lower than %d.\n",
                        MIN_READS_EXACT_MATCHING_CHARS);
                exit(EXIT_FAILURE);
            }
            PgRCParams::readsExactMatchingChars = readsExactMatchingChars;
        }

        void setPreReadsExactMatchingChars(uint16_t preReadsExactMatchingChars) {
            if (PgRCParams::preReadsExactMatchingChars != DEFAULT_UINT16_PARAM)
                return;
            if (preReadsExactMatchingChars < MIN_READS_EXACT_MATCHING_CHARS &&
                preReadsExactMatchingChars > 0) {
                fprintf(stderr, "Chars per reads exact matching cannot be lower than %d.\n",
                        MIN_READS_EXACT_MATCHING_CHARS);
                exit(EXIT_FAILURE);
            }
            PgRCParams::preReadsExactMatchingChars = preReadsExactMatchingChars;
        }

        void setMinCharsPerMismatch(uint16_t minCharsPerMismatch) {
            if (PgRCParams::minCharsPerMismatch != DEFAULT_UINT16_PARAM)
                return;
            if (minCharsPerMismatch < MIN_CHARS_PER_MISMATCH) {
                fprintf(stderr, "Chars per mismatch cannot be lower than %d.\n", MIN_CHARS_PER_MISMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCParams::minCharsPerMismatch = minCharsPerMismatch;
        }

        void setPreMatchingMode(char matchingMode) {
            if (PgRCParams::preMatchingMode != DEFAULT_CHAR_PARAM)
                return;
            PgRCParams::preMatchingMode = matchingMode;
        }

        void setMatchingMode(char matchingMode) {
            if (PgRCParams::matchingMode != DEFAULT_CHAR_PARAM)
                return;
            PgRCParams::matchingMode = matchingMode;
        }

        void setMinimalPgReverseComplementedRepeatLength(uint16_t targetPgMatchLength) {
            if (PgRCParams::targetPgMatchLength != DEFAULT_UINT16_PARAM)
                return;
            if (targetPgMatchLength < MIN_CHARS_PER_PGMATCH) {
                fprintf(stderr, "Target Pg match length cannot be lower than %d.\n", MIN_CHARS_PER_PGMATCH);
                exit(EXIT_FAILURE);
            }
            PgRCParams::targetPgMatchLength = targetPgMatchLength;
        }

        void setPgRCFileName(const string &pgRCFileName) {
            PgRCParams::pgRCFileName = pgRCFileName;
        }

        void disableRevComplPairFile() {
            PgRCParams::disableRevComplPairFileMode = true;
        }

        void setSrcFastqFile(const string &srcFastqFile) {
            PgRCParams::srcFastqFile = srcFastqFile;
        }

        void setPairFastqFile(const string &pairFastqFile) {
            PgRCParams::pairFastqFile = pairFastqFile;
        }

        void setBeginAfterStage(uint8_t skipStages) {
            if (skipStages >= endAtStage) {
                fprintf(stdout,
                        "Number of stages to skip (%d) should be smaller than a number of a stage to finish (%d).\n",
                        skipStages, endAtStage);
                exit(EXIT_FAILURE);
            }
            PgRCParams::skipStages = skipStages;
        }

        void setEndAtStage(uint8_t endAtStage) {
            if (skipStages >= endAtStage) {
                fprintf(stdout,
                        "Number of stages to skip (%d) should be smaller than a number of a stage to finish (%d).\n",
                        skipStages, endAtStage);
                exit(EXIT_FAILURE);
            }
            PgRCParams::endAtStage = endAtStage;
        }
    };

}

#endif //PGTOOLS_PGRCPARAMS_H
