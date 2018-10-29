#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultReadsMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "readsset/tools/division.h"

using namespace std;
using namespace PgTools;

int main(int argc, char *argv[])
{
    int opt; // current option
    uint8_t maxMismatches = 0;

    while ((opt = getopt(argc, argv, "m:itae?")) != -1) {
        switch (opt) {
        case 'm':
            maxMismatches = atoi(optarg);
            break;
        case 't':
            plainTextWriteMode = true;
            break;
        case 'a':
            SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
            break;
        case 'e':
            SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-m maxMismatches] [-a] [-e] [-t] \n"
                            "error_probability readssrcfile [pairsrcfile] pgFilesPrefixes\n\n",
                    argv[0]);
            fprintf(stderr, "-t write numbers in text mode\n");
            fprintf(stderr, "-a write absolute read position \n-e write mismatches as offsets from end\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 or 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string err_limit_str = argv[optind++];
    double error_limit = atof(err_limit_str.c_str());
    string srcFastqFile(argv[optind++]);
    string pairFastqFile = "";
    if (optind == argc - 4)
        pairFastqFile = argv[optind++];
    string pgFilesPrefixes(argv[optind++]);

    ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile);

    string divisionFile = pgFilesPrefixes + "_" + err_limit_str + "_bad.idxs";
    divideReads(readsIterator, divisionFile, error_limit);

    delete(readsIterator);

    // ...

    exit(EXIT_SUCCESS);
}