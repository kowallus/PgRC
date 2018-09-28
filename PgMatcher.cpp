#include <cstdlib>
#include <unistd.h>

#include "matching/matcher.h"
#include "pghelper.h"
#include "matching/DefaultPgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"

using namespace std;

static const string OFFSETS_SUFFIX = "_matched_offsets.txt";
static const string SUFFIXES_SUFFIX = "_matched_suffixes.txt";
static const string MISSED_READS_SUFFIX = "_missed.txt";

const uint_read_len_max DISABLED_PREFIX_MODE = (uint_read_len_max) -1;

void matchReadsInPgFile(const string &pgFile, const string &readsFile, const string &outPrefix,
        uint8_t maxMismatches, uint_read_len_max matchPrefixLength, bool revComplPg = false) {
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
        PgTools::approxMatchConstantLengthReads(pg, readsFile, offsetsDest, maxMismatches, matchPrefixLength,
                                                missedReadsDest, suffixesDest);
    else
        PgTools::exactMatchConstantLengthReads(pg, readsFile, offsetsDest, matchPrefixLength, missedReadsDest,
                                               suffixesDest);

    offsetsDest.close();
    missedReadsDest.close();
}

void matchPgInPgFile(const string &pgFile, const string &pgReadsFile, const string &outPrefix,
        bool revComplPg = false) {
    bool samePg = pgFile == pgReadsFile;
    PseudoGenomeBase* pgb = PgTools::openPg(pgFile);
    if (samePg)
        cout << "Reading pseudogenome..." << endl;
    else
        cout << "Reading base pseudogenome..." << endl;
    cout << "Pseudogenome length: " << pgb->getPseudoGenomeLength() << endl;
    pgb->getReadsSetProperties()->printout();
    string pg = pgb->getPseudoGenomeVirtual();
    delete pgb;
    string offsetsFile = outPrefix + OFFSETS_SUFFIX;
    std::ofstream offsetsDest(offsetsFile, std::ios::out | std::ios::binary);
    if (offsetsDest.fail()) {
        fprintf(stderr, "cannot write to offsets file %s\n", offsetsFile.c_str());
        exit(EXIT_FAILURE);
    }

    if (revComplPg) {
        //pg = pg + "XXXXXX" + PgSAHelpers::reverseComplement(pg);
        pg = PgSAHelpers::reverseComplement(pg);
        samePg = false;
    }

    pgb = PgTools::openPg(pgReadsFile);
    if (pgFile != pgReadsFile) {
        cout << "Reading pattern pseudogenome..." << endl;
        cout << "Pseudogenome length: " << pgb->getPseudoGenomeLength() << endl;
        pgb->getReadsSetProperties()->printout();
    }
    uint_pg_len_max minMatchLength = pgb->getReadsSetProperties()->maxReadLength * 1.5;

    using namespace PgTools;

    PgMatcherBase* pgmb = TemplateUserGenerator::generatePseudoGenomeUser<DefaultPgMatcher, PgMatcherBase>(pgb);

    pgmb->exactMatchPg(pg, offsetsDest, minMatchLength, samePg);
    delete pgb;
    offsetsDest.close();
}


int main(int argc, char *argv[])
{

    int opt; // current option
    bool revComplPg = false;
    uint8_t maxMismatches = 0;
    uint_read_len_max matchPrefixLength = DISABLED_PREFIX_MODE;

    while ((opt = getopt(argc, argv, "m:p:r?")) != -1) {
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
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r] [-m maxMismatches] [-p match_prefix_length] pgfile readsfile outputprefix\n\n",
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

    if (PseudoGenomePersistence::isValidPseudoGenome(readsFile))
        matchPgInPgFile(pgFile, readsFile, outPrefix, revComplPg);
    else
        matchReadsInPgFile(pgFile, readsFile, outPrefix, maxMismatches, matchPrefixLength, revComplPg);
   
    exit(EXIT_SUCCESS);
}