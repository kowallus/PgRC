#include <cstdlib>
#include <unistd.h>

#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedExtendedReadsListIterator.h"

using namespace std;
using namespace PgTools;

void rewritePg(string pgFileSrcPrefix, string pgFileDestPrefix, bool compactBytes) {
    clock_checkpoint();
    SeparatedExtendedReadsListIterator* rlIt = new SeparatedExtendedReadsListIterator(pgFileSrcPrefix);
    SeparatedPseudoGenomeOutputBuilder* builder = new SeparatedPseudoGenomeOutputBuilder(pgFileDestPrefix,
            !rlIt->isRevCompEnabled(), !rlIt->areMismatchesEnabled());
    builder->setReadsSourceIterator(rlIt);
    builder->copyPseudoGenomeHeader(pgFileSrcPrefix);
    builder->writeReadsFromIterator();
    delete(rlIt);
    builder->build();
    delete(builder);
    cout << "... rewriting completed in " << clock_millis() << " msec. " << endl;
}

int main(int argc, char *argv[])
{
    int opt; // current option
    bool compactBytes = false;

    while ((opt = getopt(argc, argv, "tcae?")) != -1) {
        switch (opt) {
        case 't':
            plainTextWriteMode = true;
            break;
        case 'c':
            compactBytes = true;
            break;
        case 'a':
            SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
            break;
        case 'e':
            SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-t] [-c] [-a] [-e] pgsrcfileprefix [pgdestfileprefix]\n\n",
                    argv[0]);
            fprintf(stderr, "-t write numbers in text mode\n-c compact output (minimal bytes used)\n");
            fprintf(stderr, "-a write absolute read position \n-e write mismatches as offsets from end\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 1) || optind < (argc - 2)) {
        fprintf(stderr, "%s: Expected 1 or 2 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string pgFileSrcPrefix(argv[optind++]);
    string pgFileDestPrefix = pgFileSrcPrefix;
    if (optind == argc - 1)
        pgFileDestPrefix = argv[optind++];

    rewritePg(pgFileSrcPrefix, pgFileDestPrefix, compactBytes);

    exit(EXIT_SUCCESS);
}