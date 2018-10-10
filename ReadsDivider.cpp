#include "utils/helper.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "readsset/iterator/DivisionReadsSetDecorators.h"
#include <stdlib.h>    /* for exit */
#include <unistd.h>

using namespace PgSAHelpers;
using namespace PgTools;

void divideReads(string srcFastqFile, string pairFastqFile, string outputFile, double error_limit) {
    clock_checkpoint();
    ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile);
    QualityDividingReadsSetIterator<uint_read_len_max> *badReadsIterator = new QualityDividingReadsSetIterator<uint_read_len_max>(readsIterator, error_limit, false);
    std::ofstream filteredIndexesDest(outputFile, std::ios::out | std::ios::binary);
    if (filteredIndexesDest.fail()) {
        fprintf(stderr, "cannot write to filtered indexes file %s\n", outputFile.c_str());
        exit(EXIT_FAILURE);
    }
    cout << "Starting division... " << endl;
    uint64_t hitCounter = 0;
    while (badReadsIterator->moveNextVirtual()) {
        hitCounter++;
        writeValue(filteredIndexesDest, badReadsIterator->getReadOriginalIndex());
    }
    writeValue(filteredIndexesDest, UINT64_MAX);
    filteredIndexesDest.close();
    cout << "Filtered " << hitCounter << " reads (out of " << (badReadsIterator->getReadOriginalIndex()) << ") in " << clock_millis() << " msec." << endl;
    delete(badReadsIterator);
    delete(readsIterator);
}

int main(int argc, char *argv[]) {

    int opt; // current option
    bool tFlag = false;

    while ((opt = getopt(argc, argv, "t?")) != -1) {
        switch (opt) {
            case 't':
                plainTextWriteMode = true;
                plainTextReadMode = true;
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-t] error_probability readssrcfile.fastq [pairsrcfile.fastq] divisionfile\n\n",
                        argv[0]);
                fprintf(stderr, "-t write numbers in text mode\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 or 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    double err_limit = atof(argv[optind++]);
    string srcFile(argv[optind++]);
    string pairFile = "";
    if (optind == argc - 2)
        pairFile = argv[optind++];
    string divisionFile(argv[optind++]);

    divideReads(srcFile, pairFile, divisionFile, err_limit);

    exit(EXIT_SUCCESS);
}
