#include <cstdlib>
#include <unistd.h>

#include "matching/matcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;

static const string OFFSETS_SUFFIX = "_matched_offsets.txt";
static const string SUFFIXES_SUFFIX = "_matched_suffixes.txt";
static const string MISSED_READS_SUFFIX = "_missed.txt";

const uint_read_len_max DISABLED_PREFIX_MODE = (uint_read_len_max) -1;

void matchReadsInPgFile(const string &pgFilePrefix, ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, const string &outPrefix,
        uint8_t maxMismatches, uint_read_len_max matchPrefixLength, bool revComplPg = false) {

    string pg = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(pgFilePrefix);
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

    string suffixesFile = outPrefix + SUFFIXES_SUFFIX;
    std::ofstream suffixesDest;
    if (matchPrefixLength != DISABLED_PREFIX_MODE) {
        suffixesDest.open(suffixesFile, std::ios::out | std::ios::binary);
        if (suffixesDest.fail()) {
            fprintf(stderr, "cannot write to suffixes file %s\n", suffixesFile.c_str());
            exit(EXIT_FAILURE);
        }
    }
    if (revComplPg)
        pg = pg + "XXXXXX" + PgSAHelpers::reverseComplement(pg);

    if (maxMismatches)
        PgTools::approxMatchConstantLengthReads(pg, readsIterator, offsetsDest, maxMismatches, matchPrefixLength,
                                                missedReadsDest, suffixesDest);
    else
        PgTools::exactMatchConstantLengthReads(pg, readsIterator, offsetsDest, matchPrefixLength, missedReadsDest,
                                               suffixesDest);

    offsetsDest.close();
    missedReadsDest.close();
}

int main(int argc, char *argv[])
{

    int opt; // current option
    bool revComplPg = false;
    uint8_t maxMismatches = 0;
    uint_read_len_max matchPrefixLength = DISABLED_PREFIX_MODE;
    string divisionFile = "";
    bool divisionComplement = false;

    while ((opt = getopt(argc, argv, "m:p:d:ctr?")) != -1) {
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
    matchReadsInPgFile(pgFilePrefix, readsIterator, outPrefix, maxMismatches, matchPrefixLength, revComplPg);
    delete(readsIterator);
   
    exit(EXIT_SUCCESS);
}