#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultPgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;

static const string OFFSETS_SUFFIX = "_matched_offsets.txt";

const uint_read_len_max DISABLED_PREFIX_MODE = (uint_read_len_max) -1;

void matchPgInPgFile(const string &destPgFile, const string &srcPgFile, bool revComplPg = false) {
    bool samePg = destPgFile == srcPgFile;
    PseudoGenomeBase* pgb = PgSAIndex::PseudoGenomePersistence::checkAndReadPseudoGenome(destPgFile);
    if (samePg)
        cout << "Reading pseudogenome..." << endl;
    else
        cout << "Reading base pseudogenome..." << endl;
    cout << "Pseudogenome length: " << pgb->getPseudoGenomeLength() << endl;
    pgb->getReadsSetProperties()->printout();
    string pg = pgb->getPseudoGenomeVirtual();
    delete pgb;
    string offsetsFile = srcPgFile + OFFSETS_SUFFIX;
    std::ofstream offsetsDest(offsetsFile, std::ios::out | std::ios::binary);
    if (offsetsDest.fail()) {
        fprintf(stderr, "cannot write to offsets file %s\n", srcPgFile.c_str());
        exit(EXIT_FAILURE);
    }

    if (revComplPg) {
        //pg = pg + "XXXXXX" + PgSAHelpers::reverseComplement(pg);
        pg = PgSAHelpers::reverseComplement(pg);
        samePg = false;
    }

    pgb = PgSAIndex::PseudoGenomePersistence::checkAndReadPseudoGenome(srcPgFile);
    if (destPgFile != srcPgFile) {
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

    while ((opt = getopt(argc, argv, "m:tr?")) != -1) {
        switch (opt) {
        case 'm':
            maxMismatches = atoi(optarg);
            break;
        case 'r':
            revComplPg = true;
            break;
        case 't':
            plainTextWriteMode = true;
            plainTextReadMode = true;
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-r] [-t] srcPgFile destPgFile\n\n",
                    argv[0]);
                fprintf(stderr, "-r match reverse compliment of pseudogenome\n-t write numbers in text mode\n\n");
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind != (argc - 2)) {
        fprintf(stderr, "%s: Expected 2 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string srcPgFile(argv[optind++]);
    string destPgFile(argv[optind++]);

    matchPgInPgFile(destPgFile, srcPgFile, revComplPg);

    exit(EXIT_SUCCESS);
}