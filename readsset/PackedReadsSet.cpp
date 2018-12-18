#include "PackedReadsSet.h"

namespace PgSAReadsSet {

    template<class ReadsSourceIterator>
    PackedReadsSet::PackedReadsSet(ReadsSourceIterator* readsIterator) {

        bool symbolOccured[UCHAR_MAX] = {0};
        uint_read_len_max minReadLength = 0;

        while (readsIterator->moveNextVirtual()) {

            properties->readsCount++;

            // analize read length
            uint_read_len_max length = readsIterator->getReadLengthVirtual();
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

            //analize symbols
            string read(readsIterator->getReadVirtual());

            for (uint_read_len_max i = 0; i < length; i++) {
                read[i] = toupper(read[i]);
                if (!symbolOccured[(unsigned char) read[i]]) {
                    symbolOccured[(unsigned char) read[i]] = true;
                    properties->symbolsCount++;
                }
            }

            lengths.push_back((uint_read_len_max) read.length());
        }

        // TODO: support variable length reads
        if (!properties->constantReadLength) {
            cout << "Unsupported: variable length reads :(" << endl;
            cout << "Trimming reads to " << minReadLength << " symbols." << endl;
            properties->constantReadLength = true;
            properties->maxReadLength = minReadLength;
            properties->allReadsLength = minReadLength * properties->readsCount;
            for(uint_reads_cnt_max i = 0; i < properties->readsCount; i++)
                lengths[i] = minReadLength;
        }

        // order symbols

        for (uint_symbols_cnt i = 0, j = 0; i < properties->symbolsCount; j++)
            if (symbolOccured[(unsigned char) j])
                properties->symbolsList[(unsigned char) (i++)] = j;

        properties->generateSymbolOrder();

        uchar symbolsPerElement = SymbolsPackingFacility<uint_ps_element_min>::maxSymbolsPerElement(properties->symbolsCount);
        sPacker = new SymbolsPackingFacility<uint_ps_element_min>(properties, symbolsPerElement);

        packedLength = (properties->maxReadLength + symbolsPerElement - 1) / symbolsPerElement;

        packedReads = new uint_ps_element_min[(size_t) packedLength * properties->readsCount];

        uint_ps_element_min* packedReadsPtr = packedReads;

        readsIterator->rewindVirtual();
        while (readsIterator->moveNextVirtual()) {
            sPacker->packSequence(readsIterator->getReadVirtual().data(), readsIterator->getReadLengthVirtual() > properties->maxReadLength?properties->maxReadLength:readsIterator->getReadLengthVirtual(), packedReadsPtr);
            packedReadsPtr += packedLength;
        }
        
    };

    PackedReadsSet::~PackedReadsSet() {
        delete[] packedReads;
        delete sPacker;
    }

    int PackedReadsSet::comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx){
        return sPacker->compareSequences(packedReads + lIdx * (size_t) packedLength, packedReads + rIdx * (size_t) packedLength, properties->maxReadLength);
    }

    int PackedReadsSet::comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx, uint_read_len_max offset) {
        return sPacker->compareSequences(packedReads + lIdx * (size_t) packedLength, packedReads + rIdx * (size_t) packedLength, offset, properties->maxReadLength - offset);
    }

    int PackedReadsSet::compareSuffixWithPrefix(uint_reads_cnt_max sufIdx, uint_reads_cnt_max preIdx, uint_read_len_max sufOffset) {
        return sPacker->compareSuffixWithPrefix(packedReads + sufIdx * (size_t) packedLength, packedReads + preIdx * (size_t) packedLength, sufOffset, properties->maxReadLength - sufOffset);
    }

    template PackedReadsSet::PackedReadsSet<ReadsSourceIteratorTemplate<uint_read_len_min>>(ReadsSourceIteratorTemplate<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<ReadsSourceIteratorTemplate<uint_read_len_std>>(ReadsSourceIteratorTemplate<uint_read_len_std>* readsIterator);
    template PackedReadsSet::PackedReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_min>>(ConcatenatedReadsSourceIterator<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_std>>(ConcatenatedReadsSourceIterator<uint_read_len_std>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTAReadsSourceIterator<uint_read_len_min>>(FASTAReadsSourceIterator<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTAReadsSourceIterator<uint_read_len_std>>(FASTAReadsSourceIterator<uint_read_len_std>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTQReadsSourceIterator<uint_read_len_min>>(FASTQReadsSourceIterator<uint_read_len_min>* readsIterator);
    template PackedReadsSet::PackedReadsSet<FASTQReadsSourceIterator<uint_read_len_std>>(FASTQReadsSourceIterator<uint_read_len_std>* readsIterator);
    
}
