#include <cstdlib>
#include <unistd.h>

#include "helper.h"

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

    if (optind > (argc - 2) || optind < (argc - 3)) {
        fprintf(stderr, "%s: Expected 2 arguments (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }
        
    string pgFile(argv[optind++]);
    string outPrefix(argv[optind++]);

    std::ifstream pgSrc(pgFile, std::ios::in | std::ios::binary);
    if (pgSrc.fail()) {
        fprintf(stderr, "cannot open pseudogenome file %s\n", pgFile.c_str());
        exit(EXIT_FAILURE);
    }
    string pg = pgTools::getPgFromPgenFile(pgSrc);
    pgSrc.close();

    string rawPgFile = outPrefix + RAWPSEUDOGENOME_EXTENSION;
    std::ofstream rawPgDest(rawPgFile, std::ios::out | std::ios::binary);
    if (rawPgDest.fail()) {
        fprintf(stderr, "cannot write to raw pseudogenome file %s\n", rawPgFile.c_str());
        exit(EXIT_FAILURE);
    }

    rawPgDest << pg.data();
    rawPgDest.close();

    exit(EXIT_SUCCESS);
}