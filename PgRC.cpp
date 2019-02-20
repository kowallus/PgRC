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

    while ((opt = getopt(argc, argv, "l:m:M:p:q:S:E:rnNtaA?")) != -1) {
        char* valPtr;
        switch (opt) {
            case 'r':
                expectedPairFile = true;
                pgRC->setRevComplPairFile();
                break;
            case 'N':
                pgRC->setSeparateNReads(true);
                break;
            case 'n':
                pgRC->setNReadsLQ(true);
                break;
            case 'l':
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
                pgRC->setMaxCharsPerMismatch(atoi(optarg));
                break;
            case 'q':
                pgRC->setError_limit_in_promils(atoi(optarg));
                break;
            case 'p':
                pgRC->setMinimalPgMatchLength(atoi(optarg));
                break;
            case 'S':
                pgRC->setSkipStages(atoi(optarg));
                break;
            case 'E':
                pgRC->setEndAtStage(atoi(optarg));
                break;
            case 't':
                plainTextWriteMode = true;
                break;
            case 'a':
                SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
                break;
            case 'A':
                SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = false;
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-m [matchingMode]exactMatchingCharsCount] [-M maxCharsPerMismatch]\n"
                                "[-r] [-n] [-N] [-a] [-A] [-t] [-s]\n"
                                "[-l [matchingMode]exactMatchingCharsCount] [-q error_probability*1000]\n"
                                "gen_quality_coef_in_%% readssrcfile [pairsrcfile] outputFileName\n\n",
                        argv[0]);
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
    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 or 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (expectedPairFile && optind != argc - 4) {
        fprintf(stderr, "Cannot use -r option without specifying a pair file.\n");
        exit(EXIT_FAILURE);
    }

    pgRC->setGen_quality_str(argv[optind++]);
    pgRC->setSrcFastqFile(argv[optind++]);
    if (optind == argc - 2)
        pgRC->setPairFastqFile(argv[optind++]);
    pgRC->setPgRCFileName(argv[optind++]);

    pgRC->executePgRCChain();

    delete(pgRC);

    exit(EXIT_SUCCESS);
}