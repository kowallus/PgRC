#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultReadsMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;

int main(int argc, char *argv[])
{
    int opt; // current option
    bool revComplPg = false;
    uint8_t maxMismatches = 0;
    uint_read_len_max matchPrefixLength = PgTools::DefaultReadsMatcher::DISABLED_PREFIX_MODE;
    string divisionFile = "";
    bool divisionComplement = false;
    string infoPrefix = "";

    while ((opt = getopt(argc, argv, "m:p:i:d:ctr?")) != -1) {
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
            break;
        case 'i':
            infoPrefix = optarg;
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
    string outDivisionFile(argv[optind++]);

    ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
            readsFile, pairFile, divisionFile, divisionComplement);
    cout << "Reading reads set\n";
    PackedReadsSet *readsSet = new PackedReadsSet(readsIterator);
    readsSet->printout();
    PgTools::DefaultReadsMatcher* matcher;
    if (maxMismatches == 0)
        matcher = new PgTools::DefaultReadsExactMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength);
    else
        matcher = new PgTools::DefaultReadsApproxMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength, maxMismatches);

    matcher->matchConstantLengthReads();
    if (infoPrefix != "")
        matcher->writeMatchesInfo(infoPrefix);

    const std::vector<uint32_t> &readsMatchPos = matcher->getReadMatchPos();
    const vector<uint8_t> &readsMismatches = matcher->getReadMismatches();
    const vector<uint_reads_cnt_max> orgIndexesMapping = ReadsSetPersistence::getReadsOriginalIndexes(divisionFile,
            divisionComplement, readsSet->getReadsSetProperties()->readsCount);

    ReadsSetPersistence::writeOutputDivision(orgIndexesMapping, readsMatchPos,
            PgTools::DefaultReadsMatcher::NOT_MATCHED_VALUE, outDivisionFile, divisionComplement);

    if (matchPrefixLength == PgTools::DefaultReadsMatcher::DISABLED_PREFIX_MODE)
        matcher->writeIntoPseudoGenome(orgIndexesMapping);

    delete(matcher);
    delete(readsSet);
    delete(readsIterator);
   
    exit(EXIT_SUCCESS);
}