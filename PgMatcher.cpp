#include <cstdlib>
#include <unistd.h>

#include "matching/matcher.h"
#include "pghelper.h"

using namespace std;

static const string OFFSETS_SUFFIX = "_matched_offsets.txt";
static const string MISSED_READS_SUFFIX = "_missed.txt";

void matchReadsInPgFile(const string &pgFile, const string &readsFile, const string &outPrefix,
        uint8_t max_mismatches, uint_read_len_max matchPrefixLength, bool revComplPg = false) {
    std::ifstream readsSrc(readsFile, std::ios::in | std::ios::binary);
    if (readsSrc.fail()) {
        fprintf(stderr, "cannot open readsfile %s\n", readsFile.c_str());
        exit(EXIT_FAILURE);
    }
    readsSrc.close();
    PseudoGenomeBase* pgb = PgTools::openPg(pgFile);
    string pg = pgb->getPseudoGenomeVirtual();
    delete pgb;
    string offsetsFile = outPrefix + OFFSETS_SUFFIX;
    std::ofstream offsetsDest(offsetsFile, std::ios::out | std::ios::binary);
    if (offsetsDest.fail()) {
        fprintf(stderr, "cannot write to offsets file %s\n", offsetsFile.c_str());
        exit(EXIT_FAILURE);
    }
    string missedReadsFile = outPrefix + MISSED_READS_SUFFIX;
    std::ofstream missedReadsDest(missedReadsFile, std::ios::out | std::ios::binary);
    if (missedReadsDest.fail()) {
        fprintf(stderr, "cannot write to missed reads file %s\n", missedReadsFile.c_str());
        exit(EXIT_FAILURE);
    }

    if (revComplPg)
        pg = pg + "XXXXXX" + PgSAHelpers::reverseComplement(pg);

    if (max_mismatches)
        approxMatchConstantLengthPatterns(pg, readsFile, offsetsDest, max_mismatches, matchPrefixLength, missedReadsDest);
    else
        exactMatchConstantLengthPatterns(pg, readsFile, offsetsDest, matchPrefixLength, missedReadsDest);

    offsetsDest.close();
    missedReadsDest.close();
}

int main(int argc, char *argv[])
{

    int opt; // current option
    bool revComplPg = false;
    uint8_t max_mismatches = 0;
    uint_read_len_max matchPrefixLength = (uint_read_len_max) -1;

    while ((opt = getopt(argc, argv, "m:p:r?")) != -1) {
        switch (opt) {
        case 'm':
            max_mismatches = atoi(optarg);
            break;
        case 'p':
            matchPrefixLength = atoi(optarg);
            break;
        case 'r':
            revComplPg = true;
        break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r] [-m max_mismatches] [-p match_prefix_length] pgfile readsfile outputprefix\n\n",
                    argv[0]);
                fprintf(stderr, "-r match reverse compliment of pseudogenome\n\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind != (argc - 3)) {
        fprintf(stderr, "%s: Expected 3 arguments (found %d) after options\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }
        
    string pgFile(argv[optind++]);
    string readsFile(argv[optind++]);
    string outPrefix(argv[optind++]);

    matchReadsInPgFile(pgFile, readsFile, outPrefix, max_mismatches, matchPrefixLength, revComplPg);
   
    exit(EXIT_SUCCESS);
}