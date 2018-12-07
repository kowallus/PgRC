#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultPgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;

int main(int argc, char *argv[])
{
    int opt; // current option
    bool revComplPg = false;
    bool dumpInfo = false;

    while ((opt = getopt(argc, argv, "rt?")) != -1) {
        switch (opt) {
        case 'r':
            revComplPg = true;
            break;
        case 't':
            plainTextWriteMode = true;
            break;
        case 'i':
            dumpInfo = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r] [-t] targetMatchLength srcPgFilePrefix targetPgFilePrefix destPgFilePrefix\n\n",
                    argv[0]);
                fprintf(stderr, "-r match reverse compliment of pseudogenome\n-t write numbers in text mode\n\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind != (argc - 4)) {
        fprintf(stderr, "%s: Expected 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    uint_pg_len_max targetMatchLength = atoi(argv[optind++]);
    string srcPgFilePrefix(argv[optind++]);
    string targetPgFilePrefix(argv[optind++]);
    string destPgFilePrefix(argv[optind++]);

    PgTools::matchPgInPgFile(srcPgFilePrefix, targetPgFilePrefix, targetMatchLength, destPgFilePrefix, revComplPg, dumpInfo);

    exit(EXIT_SUCCESS);
}