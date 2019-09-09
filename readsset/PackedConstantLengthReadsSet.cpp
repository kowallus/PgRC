#include "PackedConstantLengthReadsSet.h"

#include "tools/ReadsSetAnalizer.h"

namespace PgSAReadsSet {

    PackedConstantLengthReadsSet::PackedConstantLengthReadsSet(uint_read_len_max readLength, const char *symbolsList,
                                                               uint8_t symbolsCount) {
        properties->maxReadLength = readLength;
        properties->minReadLength = readLength;
        properties->constantReadLength = true;
        properties->symbolsCount = symbolsCount;
        strcpy(properties->symbolsList, symbolsList);
        properties->generateSymbolOrder();
        properties->readsCount = 0;
        properties->allReadsLength = 0;

        uchar symbolsPerElement = SymbolsPackingFacility::maxSymbolsPerElement(
                properties->symbolsCount);
        sPacker = SymbolsPackingFacility::getInstance(properties, symbolsPerElement);

        packedLength = (properties->maxReadLength + symbolsPerElement - 1) / symbolsPerElement;
    }

    void PackedConstantLengthReadsSet::reserve(uint_reads_cnt_max readsCount) {
        packedReads.reserve((size_t) packedLength * readsCount);
    }

    void PackedConstantLengthReadsSet::resize(uint_reads_cnt_max readsCount) {
        properties->readsCount = readsCount;
        properties->allReadsLength = readsCount * properties->minReadLength;
        packedReads.resize((size_t) packedLength * readsCount);
//        packedReads.shrink_to_fit();
    }

    void PackedConstantLengthReadsSet::addRead(const char* read, uint_read_len_max readLength) {
        if (readLength != properties->minReadLength) {
            fprintf(stderr, "Unsupported variable length reads.\n");
            exit(EXIT_FAILURE);
        }
        packedReads.resize((size_t) packedLength * ++properties->readsCount);
        properties->allReadsLength += properties->minReadLength;
        uint_ps_element_min *packedReadsPtr = packedReads.data() + (size_t) packedLength * (properties->readsCount - 1);
        sPacker->packSequence(read, readLength, packedReadsPtr);
    }

    void PackedConstantLengthReadsSet::copyRead(uint_reads_cnt_max srcIdx, uint_reads_cnt_max destIdx,
            uint_reads_cnt_max n) {
        std::copy(packedReads.begin() + (size_t) packedLength * srcIdx, packedReads.begin() + (size_t) packedLength * (srcIdx + n),
                  packedReads.begin() + (size_t) packedLength * destIdx);
    }

    void PackedConstantLengthReadsSet::copyPackedRead(const uint_ps_element_min *packedSequence,
                                                      uint_reads_cnt_max destIdx, uint_reads_cnt_max n) {
        std::copy(packedSequence, packedSequence + (size_t) packedLength * n,
                  packedReads.begin() + (size_t) packedLength * destIdx);
    }

    template<class ReadsSourceIterator>
    PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet(ReadsSourceIterator* readsIterator,
                                                                             ReadsSetProperties *properties) {
        bool ownProperties = properties == 0;
        if (ownProperties)
            properties = ReadsSetAnalizer::analizeReadsSet(readsIterator);

        if (!properties->constantReadLength) {
            cout << "Unsupported: variable length reads :(" << endl;
            cout << "Trimming reads to " << properties->minReadLength << " symbols." << endl;
            properties->constantReadLength = true;
            properties->maxReadLength = properties->minReadLength;
            properties->allReadsLength = properties->minReadLength * properties->readsCount;
        }

        PackedConstantLengthReadsSet* readsSet = new PackedConstantLengthReadsSet(properties->maxReadLength, properties->symbolsList, properties->symbolsCount);
        readsSet->reserve(properties->readsCount);

        readsIterator->rewind();
        while (readsIterator->moveNext())
            readsSet->addRead(readsIterator->getRead().data(), properties->minReadLength);

        if (ownProperties)
            delete(properties);

        return readsSet;
    }

    PackedConstantLengthReadsSet::~PackedConstantLengthReadsSet() {
        if (!sPacker->isGloballyManaged())
            delete sPacker;
    }

    int PackedConstantLengthReadsSet::comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx){
        return sPacker->compareSequences(packedReads.data() + lIdx * (size_t) packedLength, packedReads.data() + rIdx * (size_t) packedLength, properties->maxReadLength);
    }

    int PackedConstantLengthReadsSet::comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx, uint_read_len_max offset) {
        return sPacker->compareSequences(packedReads.data() + lIdx * (size_t) packedLength, packedReads.data() + rIdx * (size_t) packedLength, offset, properties->maxReadLength - offset);
    }

    int PackedConstantLengthReadsSet::compareSuffixWithPrefix(uint_reads_cnt_max sufIdx, uint_reads_cnt_max preIdx, uint_read_len_max sufOffset) {
        return sPacker->compareSuffixWithPrefix(packedReads.data() + sufIdx * (size_t) packedLength, packedReads.data() + preIdx * (size_t) packedLength, sufOffset, properties->maxReadLength - sufOffset);
    }

    int PackedConstantLengthReadsSet::compareReadWithPattern(const uint_reads_cnt_max i, const char *pattern) {
        return sPacker->compareSequenceWithUnpacked(packedReads.data() + i * (size_t) packedLength, pattern, properties->maxReadLength);
    }

    int PackedConstantLengthReadsSet::compareReadWithPattern(const uint_reads_cnt_max i, const char *pattern, int length) {
        return sPacker->compareSequenceWithUnpacked(packedReads.data() + i * (size_t) packedLength, pattern, length);
    }

    uint8_t PackedConstantLengthReadsSet::countMismatchesVsPattern(uint_reads_cnt_max i, const char *pattern, uint_read_len_max length,
                                                           uint8_t maxMismatches) {
        return sPacker->countSequenceMismatchesVsUnpacked(packedReads.data() + i * (size_t) packedLength, pattern, length, maxMismatches);
    }

    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<ReadsSourceIteratorTemplate<uint_read_len_min>>(ReadsSourceIteratorTemplate<uint_read_len_min>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<ReadsSourceIteratorTemplate<uint_read_len_std>>(ReadsSourceIteratorTemplate<uint_read_len_std>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_min>>(ConcatenatedReadsSourceIterator<uint_read_len_min>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_std>>(ConcatenatedReadsSourceIterator<uint_read_len_std>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<FASTAReadsSourceIterator<uint_read_len_min>>(FASTAReadsSourceIterator<uint_read_len_min>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<FASTAReadsSourceIterator<uint_read_len_std>>(FASTAReadsSourceIterator<uint_read_len_std>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<FASTQReadsSourceIterator<uint_read_len_min>>(FASTQReadsSourceIterator<uint_read_len_min>* readsIterator, ReadsSetProperties* properties);
    template PackedConstantLengthReadsSet* PackedConstantLengthReadsSet::loadReadsSet<FASTQReadsSourceIterator<uint_read_len_std>>(FASTQReadsSourceIterator<uint_read_len_std>* readsIterator, ReadsSetProperties* properties);

}
