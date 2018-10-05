#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "readsset/iterator/DivisionReadsSetDecorators.h"

#include "helper.h"
#include <stdlib.h>    /* for exit */
#include <unistd.h>

using namespace PgTools;

PseudoGenomeBase *preparePg(string srcFile, string pairFile, string divisionFile, bool divisionComplement) {
    PseudoGenomeBase* pgb = 0;

    if (PseudoGenomePersistence::isValidPseudoGenome(srcFile)) {
        pgb = PseudoGenomePersistence::readPseudoGenome(srcFile);
        cout << "Pseudogenome length: " << pgb->getPseudoGenomeLength() << endl;
        pgb->getReadsSetProperties()->printout();
        return pgb;
    }

    PseudoGenomeGeneratorFactory* pggf = new GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
    ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFile, pairFile, divisionFile, divisionComplement);

    PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsIterator);
    pgb = pggb->generatePseudoGenomeBase();
    delete(pggb);
    delete(readsIterator);
    delete(pggf);
    
    if (pgb == 0) {
        fprintf(stderr, "Failed generating Pg\n");
        exit(EXIT_FAILURE);    
    }
    
    return pgb;
}

int main(int argc, char *argv[]) {

    int opt; // current option
    bool tFlag = false;
    string divisionFile = "";
    bool divisionComplement = false;

    while ((opt = getopt(argc, argv, "d:ct?")) != -1) {
        switch (opt) {
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
                fprintf(stderr, "Usage: %s [-t] [-c] [-d divisionfile] readssrcfile [pairsrcfile] pgprefix\n\n",
                        argv[0]);
                fprintf(stderr, "-c use complement of reads division\n-t write numbers in text mode\n\n");
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

    PseudoGenomeBase *pgb = preparePg(srcFile, pairFile, divisionFile, divisionComplement);

    PgTools::SeparatedPseudoGenomePersistence::writePseudoGenome(pgb, idxPrefix, divisionFile, divisionComplement);

    delete(pgb);
    
    exit(EXIT_SUCCESS);
}