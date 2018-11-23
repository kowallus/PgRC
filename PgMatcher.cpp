#include <cstdlib>
#include <unistd.h>

#include "matching/DefaultPgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"

using namespace std;

static const string OFFSETS_SUFFIX = "_matched_offsets.txt";

void matchPgInPgFile(const string &destPgPrefix, const string &srcPgPrefix, bool revComplPg = false) {
    bool samePg = destPgPrefix == srcPgPrefix;
    if (samePg)
        cout << "Reading pseudogenome..." << endl;
    else
        cout << "Reading base pseudogenome..." << endl;
    PseudoGenomeHeader* pgh = 0;
    bool plainTextReadMode = false;
    PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(destPgPrefix, pgh, plainTextReadMode);
    cout << "Pseudogenome length: " << pgh->getPseudoGenomeLength() << endl;
    string pg = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(destPgPrefix);

    string offsetsFile = srcPgPrefix + OFFSETS_SUFFIX;
    std::ofstream offsetsDest(offsetsFile, std::ios::out | std::ios::binary);
    if (offsetsDest.fail()) {
        fprintf(stderr, "cannot write to offsets file %s\n", srcPgPrefix.c_str());
        exit(EXIT_FAILURE);
    }

    if (revComplPg)
        pg = PgSAHelpers::reverseComplement(pg);

    PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(srcPgPrefix, pgh, plainTextReadMode);
    if (destPgPrefix != srcPgPrefix) {
        cout << "Reading pattern pseudogenome..." << endl;
        cout << "Pseudogenome length: " << pgh->getPseudoGenomeLength() << endl;
    }


    uint_pg_len_max minMatchLength = pgh->getMaxReadLength() * 1.5;

    using namespace PgTools;

    DefaultPgMatcher matcher(srcPgPrefix);

    matcher.exactMatchPg(pg, offsetsDest, minMatchLength, samePg, revComplPg);

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