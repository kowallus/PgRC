#include <cstdlib>
#include <unistd.h>

#include "helper.h"
#include "pghelper.h"

using namespace std;

static const string RAWPSEUDOGENOME_EXTENSION = ".pg";

int main(int argc, char *argv[])
{

    int opt; // current option

    while ((opt = getopt(argc, argv, "?")) != -1) {
        switch (opt) {
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s pgfile outputprefix\n\n",
                    argv[0]);
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind != (argc - 2)) {
        fprintf(stderr, "%s: Expected 2 arguments (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }
        
    string pgFile(argv[optind++]);
    string outPrefix(argv[optind++]);

    PseudoGenomeBase* pgb = PgTools::openPg(pgFile);
    string pg = pgb->getPseudoGenomeVirtual();
    delete pgb;

    string rawPgFile = outPrefix + RAWPSEUDOGENOME_EXTENSION;
    PgSAHelpers::writeArrayToFile(rawPgFile, (void*) pg.data(), pg.length());

    exit(EXIT_SUCCESS);
}