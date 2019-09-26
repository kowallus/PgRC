#include "DefaultReadsSet.h"

namespace PgSAReadsSet {
        
    template<class ReadsSourceIterator>
    DefaultReadsSet::DefaultReadsSet(ReadsSourceIterator* readsIterator) {

        bool symbolOccured[UCHAR_MAX] = {0};
        uint_read_len_max minReadLength = 0;
        
        while (readsIterator->moveNext()) {

            properties->readsCount++;

            // analyze read length
            uint_read_len_max length = readsIterator->getReadLength();
            if (properties->maxReadLength == 0) {
                properties->maxReadLength = length;
                minReadLength = length;
            } else if (properties->maxReadLength != length) {
                properties->constantReadLength = false;
                if (properties->maxReadLength < length)
                    properties->maxReadLength = length;
                if (minReadLength > length)
                    minReadLength = length;
            }

            properties->allReadsLength += length;

            //analyze symbols
            string read(readsIterator->getRead());

            for (uint_read_len_max i = 0; i < length; i++) {
                read[i] = toupper(read[i]);
                if (!symbolOccured[(unsigned char) read[i]]) {
                    symbolOccured[(unsigned char) read[i]] = true;
                    properties->symbolsCount++;
                }
            }

            reads.push_back(read);
        }

        // order symbols

        for (uint_symbols_cnt i = 0, j = 0; i < properties->symbolsCount; j++)
            if (symbolOccured[(unsigned char) j])
                properties->symbolsList[(unsigned char) (i++)] = j;

        properties->generateSymbolOrder();

    };

    template DefaultReadsSet::DefaultReadsSet<ReadsSourceIteratorTemplate<uint_read_len_min>>(ReadsSourceIteratorTemplate<uint_read_len_min>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<ReadsSourceIteratorTemplate<uint_read_len_std>>(ReadsSourceIteratorTemplate<uint_read_len_std>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_min>>(ConcatenatedReadsSourceIterator<uint_read_len_min>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_std>>(ConcatenatedReadsSourceIterator<uint_read_len_std>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<FASTAReadsSourceIterator<uint_read_len_min>>(FASTAReadsSourceIterator<uint_read_len_min>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<FASTAReadsSourceIterator<uint_read_len_std>>(FASTAReadsSourceIterator<uint_read_len_std>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<FASTQReadsSourceIterator<uint_read_len_min>>(FASTQReadsSourceIterator<uint_read_len_min>* readsIterator);
    template DefaultReadsSet::DefaultReadsSet<FASTQReadsSourceIterator<uint_read_len_std>>(FASTQReadsSourceIterator<uint_read_len_std>* readsIterator);
    
}
