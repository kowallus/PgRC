#include <cstdlib>
#include <unistd.h>

#include "matching/ReadsMatchers.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "readsset/tools/division.h"
#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

static const char *const BAD_INFIX = "_bad";
static const char *const GOOD_INFIX = "_good";
static const char *const DIVISION_EXTENSION = ".div";

clock_t getTimeInSec(clock_t end_t, clock_t begin_t) { return ((end_t - begin_t) / CLOCKS_PER_SEC); }

using namespace std;
using namespace PgTools;

void divideGenerateAndMatch(string err_limit_str, string srcFastqFile, string pairFastqFile,
                            uint8_t targetMismatches, uint8_t maxMismatches,
                            string pgFilesPrefixes, bool revComplPairFile, bool skipIntermediateOutput) {

    string outputfile = "pgrc_res.txt";
    bool hasHeader = (bool) std::ifstream(outputfile);
    fstream fout(outputfile, ios::out | ios::binary | ios::app);
    if (!hasHeader)
        fout << "srcFastq\tpairFastq\trcPairFile\tpgPrefix\tq[%]\tm\tM\ttotal[s]\tdiv[s]\tgood[s]\tmatch[s]\tbad[s]\tpost[s]" << endl;
    clock_t start_t = clock();
    fout << srcFastqFile << "\t" << pairFastqFile << "\t" << (revComplPairFile?"yes":"no") << "\t"
        << pgFilesPrefixes << "\t" << err_limit_str << "\t"
        << (int) targetMismatches << "\t" << (int) maxMismatches << "\t";

    double error_limit = atoi(err_limit_str.c_str()) / 100.0;
    pgFilesPrefixes = pgFilesPrefixes + "_q" + (err_limit_str.length()==1?"0":"") + err_limit_str;
    string badDivisionFile = pgFilesPrefixes + BAD_INFIX + DIVISION_EXTENSION;
    string pgGoodPrefix = pgFilesPrefixes + GOOD_INFIX;
    string pgFilesPrefixesWithM = pgFilesPrefixes + "_m" + toString(targetMismatches)
            + (maxMismatches>targetMismatches?("_M" + toString(maxMismatches)):"");
    string pgMappedGoodPrefix = pgFilesPrefixesWithM + GOOD_INFIX;
    string pgMappedBadPrefix = pgFilesPrefixesWithM + BAD_INFIX;
    string mappedBadDivisionFile = pgFilesPrefixesWithM + DIVISION_EXTENSION;
    if (skipIntermediateOutput) {
        pgGoodPrefix = pgMappedGoodPrefix;
        badDivisionFile = mappedBadDivisionFile;
    }

    ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile);
    divideReads(allReadsIterator, badDivisionFile, error_limit);
    delete(allReadsIterator);
    clock_t div_t = clock();

    ReadsSourceIteratorTemplate<uint_read_len_max> *goodReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile, badDivisionFile, true, revComplPairFile);
    PseudoGenomeBase* goodPgb = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(goodReadsIterator);
    SeparatedPseudoGenomePersistence::writePseudoGenome(goodPgb, pgGoodPrefix, badDivisionFile, true, revComplPairFile);
    delete(goodPgb);
    delete(goodReadsIterator);
    clock_t good_t = clock();

    ReadsSourceIteratorTemplate<uint_read_len_max> *badReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile, badDivisionFile, false);
    cout << "Reading (div: " << badDivisionFile << ") reads set\n";
    PackedReadsSet *badReadsSet = new PackedReadsSet(badReadsIterator);
    badReadsSet->printout();
    mapReadsIntoPg(
            pgGoodPrefix, true, badReadsSet, DefaultReadsMatcher::DISABLED_PREFIX_MODE,
            targetMismatches, maxMismatches, 0,
            false, pgMappedGoodPrefix, badDivisionFile, false, mappedBadDivisionFile);
    delete(badReadsSet);
    delete(badReadsIterator);
    clock_t match_t = clock();

    ReadsSourceIteratorTemplate<uint_read_len_max> *mappedBadReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile, mappedBadDivisionFile, false, revComplPairFile);
    PseudoGenomeBase* badPgb = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(mappedBadReadsIterator);
    SeparatedPseudoGenomePersistence::writePseudoGenome(badPgb, pgMappedBadPrefix, mappedBadDivisionFile, false, revComplPairFile);
    delete(badPgb);
    delete(mappedBadReadsIterator);
    clock_t bad_t = clock();

    SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedGoodPrefix, pgMappedBadPrefix});

    fout << getTimeInSec(clock(), start_t) << "\t";
    fout << getTimeInSec(div_t, start_t) << "\t";
    fout << getTimeInSec(good_t, div_t) << "\t";
    fout << getTimeInSec(match_t, good_t) << "\t";
    fout << getTimeInSec(bad_t, match_t) << "\t";
    fout << getTimeInSec(clock(), bad_t) << endl;
}

int main(int argc, char *argv[])
{
    int opt; // current option
    uint8_t maxMismatches = 0;
    uint8_t targetMismatches = 0;
    bool skipIntermediateOutput = false;
    bool revComplPairFile = false;

    while ((opt = getopt(argc, argv, "m:M:rsitae?")) != -1) {
        switch (opt) {
        case 'r':
            revComplPairFile = true;
            break;
        case 's':
            skipIntermediateOutput = true;
            break;
        case 'm':
            targetMismatches = atoi(optarg);
            break;
        case 'M':
            maxMismatches = atoi(optarg);
            break;
        case 't':
            plainTextWriteMode = true;
            break;
        case 'a':
            SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
            break;
        case 'e':
            SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-m targetMaxMismatches] [-M allowedMaxMismatches] [-r] [-a] [-e] [-t] [-s] \n"
                            "error_probability_in_%% readssrcfile [pairsrcfile] pgFilesPrefixes\n\n",
                    argv[0]);
            fprintf(stderr, "-r reverse compliment reads in a pair file\n");
            fprintf(stderr, "-t write numbers in text mode\n-c skip intermediate output files");
            fprintf(stderr, "-a write absolute read position \n-e write mismatches as offsets from end\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 or 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    if (revComplPairFile && optind != argc - 4) {
        fprintf(stderr, "Cannot use -r option without specifying a pair file.\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    if (maxMismatches < targetMismatches) {
        fprintf(stdout, "INFO: allowedMaxMismatches set to targetMaxMismatches.\n");
        maxMismatches = targetMismatches;
    }

    string error_limit = argv[optind++];
    string srcFastqFile(argv[optind++]);
    string pairFastqFile = "";
    if (optind == argc - 2)
        pairFastqFile = argv[optind++];
    string pgFilesPrefixes(argv[optind++]);

    divideGenerateAndMatch(error_limit, srcFastqFile, pairFastqFile, targetMismatches, maxMismatches, pgFilesPrefixes,
            revComplPairFile, skipIntermediateOutput);

    exit(EXIT_SUCCESS);
}