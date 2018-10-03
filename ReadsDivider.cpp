#include "helper.h"
#include <stdlib.h>    /* for exit */
#include <unistd.h>

using namespace PgSAHelpers;

int main(int argc, char *argv[]) {

    int opt; // current option
    bool tFlag = false;

    while ((opt = getopt(argc, argv, "t?")) != -1) {
        switch (opt) {
            case 't':
                plainTextWriteMode = true;
                plainTextReadMode = true;
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-t] readssrcfile [pairsrcfile] divisionfile\n\n",
                        argv[0]);
                fprintf(stderr, "-t write numbers in text mode\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2) || optind < (argc - 3)) {
        fprintf(stderr, "%s: Expected 2 or 3 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    string srcFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 2)
        pairFile = argv[optind++];
    string divisionFile(argv[optind++]);


    exit(EXIT_SUCCESS);
}