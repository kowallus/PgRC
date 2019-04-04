#include <cstdlib>
#include <unistd.h>

#include "PgRCManager.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

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

    while ((opt = getopt(argc, argv, "c:i:l:m:M:p:q:g:S:E:dsoIvrNtaAV?")) != -1) {
        char* valPtr;
        switch (opt) {
            case 'c':
                compressionParamPresent = true;
                pgRC->setCompressionLevel(atoi(optarg));
                break;
            case 's':
                compressionParamPresent = true;
                pgRC->setSingleReadsMode();
                break;
            case 'o':
                compressionParamPresent = true;
                pgRC->setPreserveOrderMode();
                break;
            case 'I':
                compressionParamPresent = true;
                expectedPairFile = true;
                pgRC->setIgnorePairOrderInformation();
                break;
            case 'V':
                compressionParamPresent = true;
                pgRC->allowVariableParams();
                break;
            case 'd':
                decompressMode = true;
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
            case 'n':
                compressionParamPresent = true;
                pgRC->setNReadsLQ();
                break;
            case 'v':
                compressionParamPresent = true;
                pgRC->setValidationOutputMode();
                break;
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
            case 'm':
                compressionParamPresent = true;
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
                pgRC->setReadsExactMatchingChars(atoi(valPtr));
                break;
            case 'M':
                compressionParamPresent = true;
                pgRC->setMinCharsPerMismatch(atoi(optarg));
                break;
            case 'q':
                compressionParamPresent = true;
                pgRC->setError_limit_in_promils(atoi(optarg));
                break;
            case 'p':
                compressionParamPresent = true;
                pgRC->setMinimalPgMatchLength(atoi(optarg));
                break;
            case 'S':
                compressionParamPresent = true;
                pgRC->setSkipStages(atoi(optarg));
                break;
            case 'E':
                compressionParamPresent = true;
                pgRC->setEndAtStage(atoi(optarg));
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
            case 'i':
                pgRC->setSrcFastqFile(optarg);
                srcFilePresent = true;
                if (optind < (argc - 1) && argv[optind][0] != '-') {
                    pgRC->setPairFastqFile(argv[optind++]);
                    pairFilePresent = true;
                }
                break;
            case 'g':
                compressionParamPresent = true;
                pgRC->setGen_quality_str(optarg);
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-c 1<=compressionLevel<=3 (2 - default)] [-d] [-s] [-o] [-I]\n"
                                "[-i readssrcfile [pairsrcfile]] pgRCFileName\n\n"
                                ,
                        argv[0]);
                fprintf(stderr, "-d for decompression mode (supports only -i parameter for validation mode)\n");
                fprintf(stderr, "-s ignore pair information (explicit single reads mode)\n");
                fprintf(stderr, "-o preserve original order information\n\n");
                fprintf(stderr, "-I ignore order of reads in a pair (works when pairsrcfile is specified) \n\n");
                fprintf(stderr, "------------------ EXPERT(/DEVELOPER) OPTIONS ----------------\n");
                fprintf(stderr, "[-m [matchingMode]lengthOfExactMatchedReadPart] [-M minCharsPerMismatch]\n"
                                "[-l [matchingMode]lengthOfExactMatchedReadPart]\n[-p minimalExactMatchingLength]\n"
                                "[-q error_probability*1000]\n[-g gen_quality_coef_in_%%]\n"
                                "[-N] [-V] [-r] [-a] [-A] [-t]\n\n");
                fprintf(stderr, "-N reads containing N are not processed separately\n"); // -n reads containing N are low quality
                fprintf(stderr, "-l enables preliminary reads matching stage\n");
                fprintf(stderr, "Matching modes: d[s]:default; i[s]:interleaved; c[s]:copMEM ('s' suffix: shortcut after first read match)\n");
                fprintf(stderr, "(Stages: 1:division; 2:PgGenDivision; 3:Pg(good); 4:ReadsMatching; 5:Pg(bad&N); 6:PgSeqsCompression; 7:orderInfo\n");
                fprintf(stderr, "-r disable reverse compliment reads in a pair file for all PE modes\n");
                fprintf(stderr, "-v dump extra files for validation mode purposes\n-t write numbers in text mode\n");
                fprintf(stderr, "-a write absolute read position \n-A write mismatches as positions\n");
                fprintf(stderr, "-S number of stages to skip \n-E number of a stage to finish\n");
                fprintf(stderr, "\n\n");
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
    if (expectedPairFile && !pairFilePresent) {
        fprintf(stderr, "Cannot use -r or -I option without specifying a pair file.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    pgRC->setPgRCFileName(argv[optind++]);

    if (decompressMode)
        pgRC->decompressPgRC();
    else
        pgRC->executePgRCChain();

    delete(pgRC);

    exit(EXIT_SUCCESS);
}