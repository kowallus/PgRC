#include <cstdlib>
#include <unistd.h>

#include "helper.h"
#include "pghelper.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "strands/DefaultStrandDetector.h"

using namespace std;
using namespace PgTools;

int main(int argc, char *argv[])
{

    int opt; // current option
    uint_read_len_max overlap_threshold = 40;
    bool concatenated_readssrc = false;
    bool paired_reads = false;

    while ((opt = getopt(argc, argv, "k:c?")) != -1) {
        switch (opt) {
            case 'k':
                overlap_threshold = atoi(optarg);
                break;
            case 'c':
                concatenated_readssrc = true;
                paired_reads = true;
                break;
            case '?':
                default: /* '?' */
                fprintf(stderr, "Usage: %s [-k overlap_threshold] [-c concatenated source files] pgfile readssrcfile [pairsrcfile] outputprefix\n\n",
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
    if (optind == argc - 2) {
        pairFile = argv[optind++];
        paired_reads = true;
        if (concatenated_readssrc) {
            fprintf(stderr, "Cannot use -c option with [pairsrcfile] parameter present.\n");
            exit(EXIT_FAILURE);
        }
    }
    string outPrefix(argv[optind++]);

    PseudoGenomeBase* pgb = PgTools::openPg(pgFile);
    pgb->getReadsSetProperties()->printout();

    StrandDetectorBase* sdb = TemplateUserGenerator::generateReadsListUser<DefaultStrandDetector, StrandDetectorBase>(pgb);

    sdb->detectStrands(overlap_threshold, paired_reads, concatenated_readssrc);

    delete pgb;

    exit(EXIT_SUCCESS);
}