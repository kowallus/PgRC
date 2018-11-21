#include <cstdlib>
#include <unistd.h>

#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedExtendedReadsListIterator.h"

using namespace std;
using namespace PgTools;

int main(int argc, char *argv[])
{
    int opt; // current option

    while ((opt = getopt(argc, argv, "tae?")) != -1) {
        switch (opt) {
        case 't':
            plainTextWriteMode = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-t] pgfileprefix [pgfileprefix_2 [ pgfileprefix_3 [ ... ]]]\n\n",
                    argv[0]);
            fprintf(stderr, "-t write numbers in text mode\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 1)) {
        fprintf(stderr, "%s: Expected at least 1 argument after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    vector<string> pgFilePrefixes;
    do {
        pgFilePrefixes.push_back(argv[optind++]);
    } while (optind < argc);

    SeparatedPseudoGenomePersistence::dumpPgPairs(pgFilePrefixes);

    exit(EXIT_SUCCESS);
}