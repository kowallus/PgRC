#include <cstdlib>
#include <unistd.h>

#include "PgRCManager.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

#define RELEASE_DATE "2019-09-13"

using namespace std;
using namespace PgTools;

int main(int argc, char *argv[])
{
    int opt; // current option
    PgRCManager* pgRC = new PgRCManager();
    bool expectedPairFile = false;
    bool srcFilePresent = false;
    bool pairFilePresent = false;
    bool compressionParamPresent = false;
    bool decompressMode = false;

#ifndef DEVELOPER_BUILD
    logout = &null_stream;
#endif

#ifdef DEVELOPER_BUILD
    while ((opt = getopt(argc, argv, "c:i:q:g:s:M:p:l:B:E:doSIrNVvtaA?")) != -1) {
        char* valPtr;
#else
    while ((opt = getopt(argc, argv, "c:i:q:g:s:M:p:do?")) != -1) {
#endif
        switch (opt) {
            case 'c':
                compressionParamPresent = true;
                pgRC->setCompressionLevel(atoi(optarg));
                break;
            case 'i':
                pgRC->setSrcFastqFile(optarg);
                srcFilePresent = true;
                if (optind < (argc - 1) && argv[optind][0] != '-') {
                    pgRC->setPairFastqFile(argv[optind++]);
                    pairFilePresent = true;
                }
                break;
            case 'o':
                compressionParamPresent = true;
                pgRC->setPreserveOrderMode();
                break;
            case 'd':
                decompressMode = true;
                break;

            case 'q':
                compressionParamPresent = true;
                pgRC->setQualityBasedDivisionErrorLimitInPromils(atoi(optarg));
                break;
            case 'g':
                compressionParamPresent = true;
                pgRC->setPgGeneratorBasedDivisionOverlapThreshold_str(optarg);
                break;
            case 's':
                compressionParamPresent = true;
#ifdef DEVELOPER_BUILD
                valPtr = optarg + 1;
                switch (*optarg) {
                    case 'd':case 'i':case 'c':
                        if(*valPtr == 's') {
                            valPtr++;
                            pgRC->setMatchingMode(toupper(*optarg));
                        } else
                            pgRC->setMatchingMode(*optarg);
                        break;
                    default: valPtr--;
                }
                pgRC->setReadSeedLength(atoi(valPtr));
#else
                pgRC->setReadSeedLength(atoi(optarg));
#endif
                break;
            case 'M':
                compressionParamPresent = true;
                pgRC->setMinCharsPerMismatch(atoi(optarg));
                break;
            case 'p':
                compressionParamPresent = true;
                pgRC->setMinimalPgReverseComplementedRepeatLength(atoi(optarg));
                break;
#ifdef DEVELOPER_BUILD
            case 'l':
                compressionParamPresent = true;
                valPtr = optarg + 1;
                switch (*optarg) {
                    case 'd':case 'i':case 'c':
                        if(*valPtr == 's') {
                            valPtr++;
                            pgRC->setPreMatchingMode(toupper(*optarg));
                        } else
                            pgRC->setPreMatchingMode(*optarg);
                        break;
                    default: valPtr--;
                }
                pgRC->setPreReadsExactMatchingChars(atoi(valPtr));
                break;
            case 'S':
                compressionParamPresent = true;
                pgRC->setSingleReadsMode();
                break;
            case 'I':
                compressionParamPresent = true;
                expectedPairFile = true;
                pgRC->setIgnorePairOrderInformation();
                break;
            case 'r':
                compressionParamPresent = true;
                expectedPairFile = true;
                pgRC->disableRevComplPairFile();
                break;
            case 'N':
                compressionParamPresent = true;
                pgRC->doNotSeparateNReads();
                break;
            case 'V':
                compressionParamPresent = true;
                pgRC->allowVariableParams();
                break;
            case 'v':
                compressionParamPresent = true;
                pgRC->setValidationOutputMode();
                break;
            case 't':
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
                pgRC->setBeginAfterStage(atoi(optarg));
                break;
            case 'E':
                compressionParamPresent = true;
                pgRC->setEndAtStage(atoi(optarg));
                break;
#endif
            case '?':
            default: /* '?' */
                fprintf(stderr, "PgRC %d.%d: Copyright (c) 2019 Tomasz Kowalski, Szymon Grabowski: %s\n\n",
                        (int) PGRC_VERSION_MAJOR, (int) PGRC_VERSION_MINOR, RELEASE_DATE);
                fprintf(stderr, "Usage: %s [-c compressionLevel] [-i inputSrcFile [pairSrcFile]] [-o] [-d] "
                                "outputName\n\n", argv[0]);
                fprintf(stderr, "-c compression levels: 1 - fast; 2 - default; 3 - max\n");
                fprintf(stderr, "-d decompression mode\n");
                fprintf(stderr, "-o preserve original read order information\n\n");
                fprintf(stderr, "------------------ EXPERT OPTIONS ----------------\n");
                fprintf(stderr, "[-q qualityStreamErrorProbability*1000] (1000=>disable)\n"
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
                fprintf(stderr, "[-l [matchingMode]lengthOfReadSeedPartForReadsAlignmentPhase] (enables preliminary reads matching stage)\n"
                                "[-S] [-I] [-r] [-N] [-V] [-v] [-t] [-a] [-A]\n"
                                "[-B numberOfStagesToSkip] [-E numberOfAStageToEnd]\n\n");
                fprintf(stderr, "-S ignore pair information (explicit single reads mode)\n");
                fprintf(stderr, "-I ignore order of reads in a pair (works when pairSrcFile is specified)\n");
                fprintf(stderr, "-r disable reverse compliment reads in a pair file for all PE modes\n");
                fprintf(stderr, "-N reads containing N are not processed separately\n");
                fprintf(stderr, "-V allow variable (auto-adjusting) parameters during processing\n");
                fprintf(stderr, "-v dump extra files for validation mode purposes "
                                "(decompression supports -i parameter in validation mode)\n"
                                "-t write numbers in text mode\n");
                fprintf(stderr, "-a write absolute read position \n-A write mismatches as positions\n");
                fprintf(stderr, "Stages: 1:QualDivision; 2:PgGenDivision; 3:Pg(HQ); 4:ReadsMatching; 5:Pg(LQ&N); 6:OrderInfo; 7:PgSequences\n\n");
#endif
                fprintf(stderr, "The order of all selected options is arbitrary.\n\n");
                exit(EXIT_FAILURE);
        }
    }
    if (optind > (argc - 1) || optind < (argc - 1)) {
        fprintf(stderr, "%s: Expected 1 argument after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (decompressMode && compressionParamPresent) {
        fprintf(stderr, "Cannot use compression options in decompression mode.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!srcFilePresent && !decompressMode) {
        fprintf(stderr, "Input file(s) not specified.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
#ifdef DEVELOPER_BUILD
    if (expectedPairFile && !pairFilePresent) {
        fprintf(stderr, "Cannot use -r or -I option without specifying a pair file.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
#endif

    pgRC->setPgRCFileName(argv[optind++]);

    if (decompressMode)
        pgRC->decompressPgRC();
    else
        pgRC->executePgRCChain();

    delete(pgRC);

    exit(EXIT_SUCCESS);
}