#include "division.h"

#include "../iterator/DivisionReadsSetDecorators.h"

using namespace PgSAHelpers;

void PgTools::divideReads(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, string outputFile, double error_limit) {
    clock_checkpoint();
    QualityDividingReadsSetIterator<uint_read_len_max> *badReadsIterator = new QualityDividingReadsSetIterator<uint_read_len_max>(readsIterator, error_limit, false);
    std::ofstream filteredIndexesDest(outputFile, std::ios::out | std::ios::binary);
    if (filteredIndexesDest.fail()) {
        fprintf(stderr, "cannot write to filtered indexes file %s\n", outputFile.c_str());
        exit(EXIT_FAILURE);
    }
    cout << "Starting division... " << endl;
    writeReadMode(filteredIndexesDest, PgSAHelpers::plainTextWriteMode);
    uint64_t hitCounter = 0;
    while (badReadsIterator->moveNextVirtual()) {
        hitCounter++;
        writeValue(filteredIndexesDest, badReadsIterator->getReadOriginalIndex());
    }
    writeValue(filteredIndexesDest, UINT64_MAX);
    filteredIndexesDest.close();
    cout << "Filtered " << hitCounter << " reads (out of " << (badReadsIterator->getReadOriginalIndex()) << ") in " << clock_millis() << " msec." << endl << endl;
    delete(badReadsIterator);
}