#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "readsset/iterator/DivisionReadsSetDecorators.h"

#include "utils/helper.h"
#include <stdlib.h>    /* for exit */
#include <unistd.h>

using namespace PgTools;

void prepareAndWritePg(string srcFile, string pairFile, string divisionFile, bool divisionComplement, string idxPrefix) {
    PseudoGenomeBase* pgb = 0;

    vector<uint_reads_cnt_max> indexesMapping = {};
    if (PseudoGenomePersistence::isValidPseudoGenome(srcFile)) {
        pgb = PseudoGenomePersistence::readPseudoGenome(srcFile);
        cout << "Pseudogenome length: " << pgb->getPseudoGenomeLength() << endl;
        pgb->getReadsSetProperties()->printout();

    } else {
        ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
                srcFile, pairFile, divisionFile, divisionComplement);
        pgb = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(readsIterator);
        indexesMapping = readsIterator->retainVisitedIndexesMapping();
        delete (readsIterator);
    }

    PgTools::SeparatedPseudoGenomePersistence::writePseudoGenome(pgb, idxPrefix, indexesMapping);
    delete(pgb);
}

int main(int argc, char *argv[]) {

    int opt; // current option
    bool tFlag = false;
    string divisionFile = "";
    bool divisionComplement = false;

    while ((opt = getopt(argc, argv, "d:cta?")) != -1) {
        switch (opt) {
            case 'c':
                divisionComplement = true;
                break;
            case 'd':
                divisionFile = optarg;
                break;
            case 't':
                plainTextWriteMode = true;
                break;
            case 'a':
                SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = true;
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-t] [-c] [-a] [-d divisionfile] readssrcfile [pairsrcfile] pgprefix\n\n",
                        argv[0]);
                fprintf(stderr, "-c use complement of reads division\n-t write numbers in text mode\n");
                fprintf(stderr, "-a write absolute read position \n");
                fprintf(stderr, "\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (divisionFile == "" && divisionComplement) {
        fprintf(stderr, "%s: Division complement option (-c) cannot be used without division file (-d).\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    if (optind > (argc - 2) || optind < (argc - 3)) {
        fprintf(stderr, "%s: Expected 2 or 3 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    string srcFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 2)
        pairFile = argv[optind++];
    string idxPrefix(argv[optind++]);

    prepareAndWritePg(srcFile, pairFile, divisionFile, divisionComplement, idxPrefix);

    exit(EXIT_SUCCESS);
}