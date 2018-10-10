#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultReadsMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;

void matchReadsInPgFile(const string &pgFilePrefix, ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, const string &outPrefix,
        uint8_t maxMismatches, uint_read_len_max matchPrefixLength, bool revComplPg = false, bool infoDump = false) {

    cout << "Reading reads set\n";
    PackedReadsSet *readsSet = new PackedReadsSet(readsIterator);
    readsSet->printout();

    PgTools::DefaultReadsMatcher matcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength, maxMismatches);

    matcher.matchConstantLengthReads();

    if (infoDump)
        matcher.writeMatchesInfo(outPrefix);

    delete (readsSet);
}

int main(int argc, char *argv[])
{

    int opt; // current option
    bool revComplPg = false;
    uint8_t maxMismatches = 0;
    uint_read_len_max matchPrefixLength = PgTools::DefaultReadsMatcher::DISABLED_PREFIX_MODE;
    string divisionFile = "";
    bool divisionComplement = false;
    bool infoDump = false;

    while ((opt = getopt(argc, argv, "m:p:d:ctri?")) != -1) {
        switch (opt) {
        case 'm':
            maxMismatches = atoi(optarg);
            break;
        case 'p':
            matchPrefixLength = atoi(optarg);
            break;
        case 'r':
            revComplPg = true;
            break;
        case 'c':
            divisionComplement = true;
            break;
        case 'd':
            divisionFile = optarg;
            break;
        case 't':
            plainTextWriteMode = true;
            plainTextReadMode = true;
            break;
        case 'i':
            infoDump = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r] [-m maxMismatches] [-p match_prefix_length] [-t] [-c] [-d divisionfile] readssrcfile [pairsrcfile] pgfileprefix outputdivisionfile\n\n",
                    argv[0]);
                fprintf(stderr, "-r match reverse compliment of pseudogenome\n-c use complement of reads division\n-t write numbers in text mode\n\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 or 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string readsFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 3)
        pairFile = argv[optind++];
    string pgFilePrefix(argv[optind++]);
    string outPrefix(argv[optind++]);

    ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
            readsFile, pairFile, divisionFile, divisionComplement);
    matchReadsInPgFile(pgFilePrefix, readsIterator, outPrefix, maxMismatches, matchPrefixLength, revComplPg, infoDump);
    delete(readsIterator);
   
    exit(EXIT_SUCCESS);
}