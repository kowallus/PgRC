#include <cstdlib>
#include <unistd.h>

#include "pgrc/pgrc-encoder.h"
#include "pgrc/pgrc-decoder.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include <omp.h>

#define RELEASE_DATE "2024-09-26"

using namespace std;
using namespace PgTools;

void printVersion(bool details) {
    string date = RELEASE_DATE;
    if (!details)
        date.resize(4);
    fprintf(stderr, "PgRC %d.%d: Copyright (c) Tomasz Kowalski, Szymon Grabowski: %s\n\n",
                        (int) PGRC_VERSION_MAJOR, (int) PGRC_VERSION_MINOR, date.c_str());
}

int main(int argc, char *argv[])
{
    int opt; // current option
    PgRCParams* params = new PgRCParams();
    bool expectedPairFile = false;
    bool srcFilePresent = false;
    bool pairFilePresent = false;
    bool compressionParamPresent = false;
    bool decompressMode = false;

    numberOfThreads = omp_get_num_procs();

#ifndef DEVELOPER_BUILD
    NullBuffer null_buffer;
    std::ostream null_stream(&null_buffer);
    logout = &null_stream;
#endif

#ifdef DEVELOPER_BUILD
    while ((opt = getopt(argc, argv, "c:t:i:q:g:s:M:p:l:B:E:C:doSIrNRVTaAQvh?")) != -1) {
        char* valPtr;
#else
    while ((opt = getopt(argc, argv, "t:i:q:g:s:M:p:doQvh?")) != -1) {
#endif
        switch (opt) {
            case 'i':
                params->setSrcFastqFile(optarg);
                srcFilePresent = true;
                if (optind < (argc - 1) && argv[optind][0] != '-') {
                    params->setPairFastqFile(argv[optind++]);
                    pairFilePresent = true;
                }
                break;
            case 'o':
                compressionParamPresent = true;
                params->setPreserveOrderMode();
                break;
            case 'd':
                decompressMode = true;
                break;
            case 't':
                numberOfThreads = atoi(optarg);
                break;
            case 'q':
                compressionParamPresent = true;
                params->setQualityBasedDivisionErrorLimitInPromils(atoi(optarg));
                break;
            case 'Q':
                compressionParamPresent = true;
                params->disableSimplifiedSuffixMode4QualityBasedDivision();
            case 'g':
                compressionParamPresent = true;
                params->setPgGeneratorBasedDivisionOverlapThreshold_str(optarg);
                break;
            case 's':
                compressionParamPresent = true;
#ifdef DEVELOPER_BUILD
                valPtr = optarg + 1;
                switch (*optarg) {
                    case 'd':case 'i':case 'c':
                        if(*valPtr == 's') {
                            valPtr++;
                            params->setMatchingMode(toupper(*optarg));
                        } else
                            params->setMatchingMode(*optarg);
                        break;
                    default: valPtr--;
                }
                params->setReadSeedLength(atoi(valPtr));
#else
                params->setReadSeedLength(atoi(optarg));
#endif
                break;
            case 'M':
                compressionParamPresent = true;
                params->setMinCharsPerMismatch(atoi(optarg));
                break;
            case 'p':
                compressionParamPresent = true;
                params->setMinimalPgReverseComplementedRepeatLength(atoi(optarg));
                break;
#ifdef DEVELOPER_BUILD
            case 'c':
                compressionParamPresent = true;
                params->setCompressionLevel(atoi(optarg));
                break;
            case 'C':
                compressionParamPresent = true;
                setAutoSelectorLevel(atoi(optarg));
                break;
            case 'l':
                compressionParamPresent = true;
                valPtr = optarg + 1;
                switch (*optarg) {
                    case 'd':case 'i':case 'c':
                        if(*valPtr == 's') {
                            valPtr++;
                            params->setPreMatchingMode(toupper(*optarg));
                        } else
                            params->setPreMatchingMode(*optarg);
                        break;
                    default: valPtr--;
                }
                params->setPreReadsExactMatchingChars(atoi(valPtr));
                break;
            case 'S':
                compressionParamPresent = true;
                params->setSingleReadsMode();
                break;
            case 'I':
                compressionParamPresent = true;
                expectedPairFile = true;
                params->setIgnorePairOrderInformation();
                break;
            case 'r':
                compressionParamPresent = true;
                expectedPairFile = true;
                params->disableRevComplPairFile();
                break;
            case 'N':
                compressionParamPresent = true;
                params->doNotSeparateNReads();
                break;
            case 'R':
                compressionParamPresent = true;
                params->allowVariableParams();
                break;
            case 'V':
                params->setValidationOutputMode();
                break;
            case 'T':
                compressionParamPresent = true;
                plainTextWriteMode = true;
                break;
            case 'a':
                compressionParamPresent = true;
                SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
                break;
            case 'A':
                compressionParamPresent = true;
                SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = false;
                break;
            case 'B':
                compressionParamPresent = true;
                params->setBeginAfterStage(atoi(optarg));
                break;
            case 'E':
                compressionParamPresent = true;
                params->setEndAtStage(atoi(optarg));
                break;
#endif
            case 'v':
                printVersion(true);
                exit(EXIT_SUCCESS);
            case '?':
            case 'h':
            default: /* '?' */
                printVersion(false);
                fprintf(stderr, "Usage: %s [-i seqSrcFile [pairSrcFile]] [-t noOfThreads]"
                                "\n[-o] [-d] archiveName\n\n", argv[0]);
                fprintf(stderr, "\t-d decompression mode\n");
                fprintf(stderr, "\t-o preserve original read order information\n");
                fprintf(stderr, "\t-t number of threads used (%d - default)\n", numberOfThreads);
                fprintf(stderr, "\t-h print full command help and exit\n");
                fprintf(stderr, "\t-v print version number and exit\n");
                fprintf(stderr, "\n------------------ EXPERT OPTIONS ----------------\n");
                fprintf(stderr, "[-q qualityStreamErrorProbability*1000] (1000=>disable)\n"
                                "[-Q] disable simplified quality estimation mode\n"
                                "[-g generatorBasedQualityCoefficientIn_%%] (0=>disable; 'ov' param in the paper)\n"
                                "[-s "
#ifdef DEVELOPER_BUILD
                                "[matchingMode]"
#endif
                                "lengthOfReadSeedPartForReadsAlignmentPhase]\n"
                                "[-M minimalNumberOfCharsPerMismatchForReadsAlignmentPhase]\n"
                                "[-p minimalReverseComplementedRepeatLength]\n\n");
#ifdef DEVELOPER_BUILD
                fprintf(stderr, "Matching modes: d[s]:default; i[s]:interleaved; c[s]:copMEM ('s' suffix: shortcut after first read match)\n");
                fprintf(stderr, "------------------ DEVELOPER OPTIONS ----------------\n");
                fprintf(stderr, "[-c backendCompressionLevel] 1 - fast; 2 - default; 3 - max\n");
                fprintf(stderr, "[-C backendCompressionAutoSelectorLevel] 0 - default\n");
                fprintf(stderr, "[-l [matchingMode]lengthOfReadSeedPartForReadsAlignmentPhase] (enables preliminary reads matching stage)\n"
                                "[-S] [-I] [-r] [-N] [-V] [-v] [-t] [-a] [-A]\n"
                                "[-B numberOfStagesToSkip] [-E numberOfAStageToEnd]\n\n");
                fprintf(stderr, "-S ignore pair information (explicit single reads mode)\n");
                fprintf(stderr, "-I ignore order of reads in a pair (works when pairSrcFile is specified)\n");
                fprintf(stderr, "-r disable reverse compliment reads in a pair file for all PE modes\n");
                fprintf(stderr, "-N reads containing N are not processed separately\n");
                fprintf(stderr, "-R allow variable (auto-adjusting) parameters during processing\n");
                fprintf(stderr, "-V dump extra files for validation mode and development purposes "
                                "(decompression supports -i parameter in validation mode)\n"
                                "-T write numbers in text mode\n");
                fprintf(stderr, "-a write absolute read position \n-A write mismatches as positions\n");
                fprintf(stderr, "Stages: 1:QualDivision; 2:PgGenDivision; 3:Pg(HQ); 4:ReadsMatching; 5:Pg(LQ&N); 6:OrderInfo; 7:PgSequences\n\n");
#endif
                fprintf(stderr, "The order of all selected options is arbitrary.\n\n");
                exit(EXIT_FAILURE);
        }
    }
    if (optind > (argc - 1) || optind < (argc - 1)) {
        fprintf(stderr, "%s: Expected 1 argument after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -h' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (decompressMode && compressionParamPresent) {
        fprintf(stderr, "Cannot use compression options in decompression mode.\n");
        fprintf(stderr, "try '%s -h' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!srcFilePresent && !decompressMode) {
        fprintf(stderr, "Input file(s) not specified.\n");
        fprintf(stderr, "try '%s -h' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
#ifdef DEVELOPER_BUILD
    if (expectedPairFile && !pairFilePresent) {
        fprintf(stderr, "Cannot use -r or -I option without specifying a pair file.\n");
        fprintf(stderr, "try '%s -h' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
#endif
    if (numberOfThreads <= 0) {
        fprintf(stderr, "The number of threads must be positive.\n");
        exit(EXIT_FAILURE);
    }
    omp_set_num_threads(numberOfThreads);

    params->setPgRCFileName(argv[optind++]);

    if (decompressMode) {
        PgRCDecoder decoder(params);
        decoder.decompressPgRC();
    }
    else {
        PgRCEncoder encoder(params);
        encoder.executePgRCChain();
    }

    delete(params);

    exit(EXIT_SUCCESS);
}