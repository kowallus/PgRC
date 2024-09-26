#include "ReadsSetAnalyzer.h"

namespace PgReadsSet {

    template<class ReadsSourceIterator>
    ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet(ReadsSourceIterator *readsIterator) {
        ReadsSetProperties* properties = new ReadsSetProperties();
        bool symbolOccured[UCHAR_MAX] = {0};

        while (readsIterator->moveNext()) {

            properties->readsCount++;

            // analyze read length
            uint_read_len_max currLength = readsIterator->getReadLength();
            if (properties->maxReadLength == 0) {
                properties->maxReadLength = currLength;
                properties->minReadLength = currLength;
            } else if (properties->maxReadLength != currLength) {
                properties->constantReadLength = false;
                if (properties->maxReadLength < currLength)
                    properties->maxReadLength = currLength;
                if (properties->minReadLength > currLength)
                    properties->minReadLength = currLength;
            }
            properties->allReadsLength += currLength;

            //analyze symbols
            string read(readsIterator->getRead());

            for (uint_read_len_max i = 0; i < currLength; i++) {
                read[i] = toupper(read[i]);
                if (!symbolOccured[(unsigned char) read[i]]) {
                    symbolOccured[(unsigned char) read[i]] = true;
                    properties->symbolsCount++;
                }
            }
        }

        // order symbols

        for (uint_symbols_cnt i = 0, j = 0; i < properties->symbolsCount; j++)
            if (symbolOccured[(unsigned char) j])
                properties->symbolsList[(unsigned char) (i++)] = j;

        properties->generateSymbolOrder();

        return properties;
    }

    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<ReadsSourceIteratorTemplate<uint_read_len_min>>(ReadsSourceIteratorTemplate<uint_read_len_min>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<ReadsSourceIteratorTemplate<uint_read_len_std>>(ReadsSourceIteratorTemplate<uint_read_len_std>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_min>>(ConcatenatedReadsSourceIterator<uint_read_len_min>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_std>>(ConcatenatedReadsSourceIterator<uint_read_len_std>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<FASTAReadsSourceIterator<uint_read_len_min>>(FASTAReadsSourceIterator<uint_read_len_min>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<FASTAReadsSourceIterator<uint_read_len_std>>(FASTAReadsSourceIterator<uint_read_len_std>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<FASTQReadsSourceIterator<uint_read_len_min>>(FASTQReadsSourceIterator<uint_read_len_min>* readsIterator);
    template ReadsSetProperties* ReadsSetAnalyzer::analyzeReadsSet<FASTQReadsSourceIterator<uint_read_len_std>>(FASTQReadsSourceIterator<uint_read_len_std>* readsIterator);

}