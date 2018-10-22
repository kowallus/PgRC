#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultReadsMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;
using namespace PgTools;

int main(int argc, char *argv[])
{
    int opt; // current option
    bool revComplPg = false;
    uint8_t maxMismatches = 0;
    uint8_t minMismatches = 0;
    uint_read_len_max matchPrefixLength = DefaultReadsMatcher::DISABLED_PREFIX_MODE;
    string divisionFile = "";
    bool divisionComplement = false;
    bool dumpInfo = false;

    while ((opt = getopt(argc, argv, "m:n:p:d:citaer?")) != -1) {
        switch (opt) {
        case 'm':
            maxMismatches = atoi(optarg);
            break;
        case 'n':
            minMismatches = atoi(optarg);
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
            dumpInfo = true;
            break;
        case 'a':
            SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
            break;
        case 'e':
            SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r] [-m maxMismatches] [-n expectedMinMismatches]\n"
                            "[-p match_prefix_length] [-a] [-e] [-t] [-i] [-c] [-d divisionfile]\n"
                            "readssrcfile [pairsrcfile] pgsrcfileprefix outputdivisionfile pgdestfileprefix\n\n",
                    argv[0]);
            fprintf(stderr, "-r match reverse compliment of pseudogenome\n-c use complement of reads division\n"
                            "-t write numbers in text mode\n-i dump verbose matching info\n");
            fprintf(stderr, "-a write absolute read position \n-e write mismatches as offsets from end\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (minMismatches > maxMismatches) {
        fprintf(stderr, "Min mismatches (%d) should not be higher than max mismatches (%d)\n", minMismatches, maxMismatches);

        exit(EXIT_FAILURE);
    }

    if (optind > (argc - 4) || optind < (argc - 5)) {
        fprintf(stderr, "%s: Expected 4 or 5 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string readsFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 4)
        pairFile = argv[optind++];
    string pgFilePrefix(argv[optind++]);
    string outDivisionFile(argv[optind++]);
    string pgDestFilePrefix(argv[optind++]);

    ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
            readsFile, pairFile, divisionFile, divisionComplement);
    cout << "Reading reads set\n";
    PackedReadsSet *readsSet = new PackedReadsSet(readsIterator);
    readsSet->printout();
    DefaultReadsMatcher* matcher;
    if (maxMismatches == 0)
        matcher = new DefaultReadsExactMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength);
    else
        matcher = new DefaultReadsApproxMatcher(pgFilePrefix, revComplPg, readsSet, matchPrefixLength, maxMismatches, minMismatches);

    matcher->matchConstantLengthReads();
    if (dumpInfo)
        matcher->writeMatchesInfo(pgDestFilePrefix);

    const std::vector<uint32_t> &readsMatchPos = matcher->getReadMatchPos();
    const vector<uint8_t> &readsMismatches = matcher->getReadMismatches();
    const vector<uint_reads_cnt_max> orgIndexesMapping = ReadsSetPersistence::getReadsOriginalIndexes(divisionFile,
            divisionComplement, readsSet->getReadsSetProperties()->readsCount);

    ReadsSetPersistence::writeOutputDivision(orgIndexesMapping, readsMatchPos,
            DefaultReadsMatcher::NOT_MATCHED_VALUE, outDivisionFile, divisionComplement);

    if (matchPrefixLength == DefaultReadsMatcher::DISABLED_PREFIX_MODE)
        matcher->writeIntoPseudoGenome(pgDestFilePrefix, orgIndexesMapping);

    delete(matcher);
    delete(readsSet);
    delete(readsIterator);
   
    exit(EXIT_SUCCESS);
}