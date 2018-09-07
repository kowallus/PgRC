#include <cstdlib>
#include <unistd.h>

#include "helper.h"
#include "pghelper.h"

using namespace std;

int main(int argc, char *argv[])
{

    int opt; // current option
    int overlap_threshold = 40;

    while ((opt = getopt(argc, argv, "k?")) != -1) {
        switch (opt) {
            case 'k':
                overlap_threshold = atoi(optarg);
                break;
            case '?':
                default: /* '?' */
                fprintf(stderr, "Usage: %s [-k overlap_threshold] pgfile readssrcfile [pairsrcfile] outputprefix\n\n",
                    argv[0]);
                fprintf(stderr, "\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 or 4 arguments after options (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string pgFile(argv[optind++]);
    string srcFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 2)
        pairFile = argv[optind++];
    string outPrefix(argv[optind++]);

    PseudoGenomeBase* pgb = pgTools::openPg(pgFile);
    pgb->getReadsSetProperties()->printout();

    delete pgb;

    exit(EXIT_SUCCESS);
}