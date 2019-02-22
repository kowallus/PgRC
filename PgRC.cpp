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

    while ((opt = getopt(argc, argv, "i:l:m:M:p:q:g:S:E:drnNtaA?")) != -1) {
        char* valPtr;
        switch (opt) {
            case 'd':
                decompressMode = true;
                break;
            case 'r':
                compressionParamPresent = true;
                expectedPairFile = true;
                pgRC->setRevComplPairFile();
                break;
            case 'N':
                compressionParamPresent = true;
                pgRC->setSeparateNReads(true);
                break;
            case 'n':
                compressionParamPresent = true;
                pgRC->setNReadsLQ(true);
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
                pgRC->setMaxCharsPerMismatch(atoi(optarg));
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
                fprintf(stderr, "Usage: %s [-m [matchingMode]exactMatchingCharsCount] [-M maxCharsPerMismatch]\n"
                                "[-d] [-r] [-n] [-N] [-a] [-A] [-t] [-s]\n"
                                "[-l [matchingMode]exactMatchingCharsCount] [-q error_probability*1000]\n"
                                "[-g gen_quality_coef_in_%%] -i readssrcfile [pairsrcfile] pgRCFileName\n\n",
                        argv[0]);
                fprintf(stderr, "-d for decompression mode (supports only -i parameter for validation)\n");
                fprintf(stderr, "-r reverse compliment reads in a pair file\n");
                fprintf(stderr, "-n reads containing N are low quality\n-N reads containing N are processed separately\n");
                fprintf(stderr, "-t write numbers in text mode\n");
                fprintf(stderr, "-a write absolute read position \n-A write mismatches as positions\n");
                fprintf(stderr, "-S number of stages to skip \n-E number of a stage to finish\n");
                fprintf(stderr, "-l enables preliminary reads matching stage\n");
                fprintf(stderr, "Matching modes: d[s]:default; i[s]:interleaved; c[s]:copMEM ('s' suffix: shortcut after first read match)\n");
                fprintf(stderr, "(Stages: 1:division; 2:PgGenDivision; 3:Pg(good); 4:ReadsMatching; 5:Pg(bad); 6:PgMatching; 7:pairDump\n");
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
    if (!srcFilePresent) {
        fprintf(stderr, "Input file(s) not specified.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (expectedPairFile && !pairFilePresent) {
        fprintf(stderr, "Cannot use -r option without specifying a pair file.\n");
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