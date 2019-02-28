#include <cstdlib>
#include <unistd.h>

#include "matching/SimplePgMatcher.h"
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
    uint32_t minMatchLength = UINT32_MAX;

    while ((opt = getopt(argc, argv, "m:rt?")) != -1) {
        switch (opt) {
        case 'm':
            minMatchLength = atoi(optarg);
            break;
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
            fprintf(stderr, "Usage: %s [-r] [-t] [-m minMatchLength] targetMatchLength srcPgFilePrefix targetPgFilePrefix destPgFilePrefix\n\n",
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
    clock_t pgMatcher_start = clock();
    string pgSeq1 = PgTools::SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence("1");
    string pgSeq2 = PgTools::SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence("2");
    PgTools::SimplePgMatcher::matchPgsInPg(pgSeq1, pgSeq2, srcPgFilePrefix, targetPgFilePrefix, targetMatchLength,
                                           minMatchLength);
    cout << "... pg matching completed in " << clock_millis(pgMatcher_start) << " msec. " << endl;
//    PgTools::DefaultPgMatcher::matchPgInPgFile(srcPgFilePrefix, targetPgFilePrefix, targetMatchLength, destPgFilePrefix, revComplPg, dumpInfo);

    exit(EXIT_SUCCESS);
}