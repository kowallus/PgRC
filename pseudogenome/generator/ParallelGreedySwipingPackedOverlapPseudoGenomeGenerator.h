#ifndef PARALLELGREEDYSWIPINGPACKEDOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
#define PARALLELGREEDYSWIPINGPACKEDOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include "AbstractOverlapPseudoGenomeGenerator.h"
#include "../../readsset/PackedConstantLengthReadsSet.h"
#include <algorithm>
#include <deque>
#include <functional>

#define MAX_BLOCK_PREFIX_LENGTH 4
#define MAX_SYMBOLS_COUNT 5
#define MAX_BLOCKS_COUNT 625 // MAX_BLOCK_PREFIX_LENGTH * |"ACGTN"|

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt >
    class ParallelGreedySwipingPackedOverlapGeneratorTemplate: public AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>
    {
    private:

        bool ownReadsSet = false;
        PackedConstantLengthReadsSet* packedReadsSet = 0;

        vector<uint_reads_cnt> sortedReadsIdxs;
        vector<uint_reads_cnt> sortedSuffixIdxs;
        const uint_reads_cnt* sortedSuffixIdxsPtr;

        const uint8_t blockPrefixLength = 3;
        char blockPrefixes[MAX_BLOCKS_COUNT][MAX_SYMBOLS_COUNT] = { 0 };
        uint16_t blocksCount;

        uint_reads_cnt sortedReadsBlockPos[MAX_BLOCKS_COUNT + 1];
        uint_reads_cnt sortedReadsCount[MAX_BLOCKS_COUNT] = { 0 };

        uint_reads_cnt sortedSuffixBlockPlusSymbolPos[MAX_BLOCKS_COUNT + 1][MAX_SYMBOLS_COUNT + 1];
        uint_reads_cnt sortedSuffixBlockPos[MAX_BLOCKS_COUNT + 1];

        uint16_t threadStartBlock[UINT8_MAX] = { 0 };

        struct PackedReadsComparator {
            ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>* myGenerator;
            PackedReadsComparator(ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>* generator) : myGenerator(generator) {};

            bool operator() (uint_reads_cnt lIncIdx, uint_reads_cnt rIncIdx) const {
                return myGenerator->compareReads(lIncIdx, rIncIdx) < 0;
            }
        };

        struct PackedReadVsPatternComparator {
            const uint_reads_cnt PATTERN_INDEX = -1;

            ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>* myGenerator;
            char* pattern;
            int patLength;
            PackedReadVsPatternComparator(ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>* generator,
                    char* pattern) : myGenerator(generator), pattern(pattern), patLength(strlen(pattern)) {};

            bool operator() (uint_reads_cnt lIncIdx, uint_reads_cnt rIncIdx) const {
                if (lIncIdx == PATTERN_INDEX)
                    return myGenerator->packedReadsSet->compareReadWithPattern(rIncIdx - 1, pattern, patLength) >= 0;
                else if (rIncIdx == PATTERN_INDEX)
                    return myGenerator->packedReadsSet->compareReadWithPattern(lIncIdx - 1, pattern, patLength) < 0;
                return myGenerator->compareReads(lIncIdx, rIncIdx) < 0;
            }
        };

        int compareReads(uint_reads_cnt lIncIdx, uint_reads_cnt rIncIdx);

        uchar getSymbolOrderFromRead(uint_reads_cnt incIdx, uint_read_len offset);

        void updateSuffixQueue(uchar suffixGroup, uint_read_len suffixOffset,
                uint_reads_cnt* ssiSymbolIdx, uint_reads_cnt* ssiSymbolEnd, deque<uchar> &ssiOrder);

        int compareSuffixes(uchar lSymOrder, uchar rSymOrder, uint_read_len offset, uint_reads_cnt* ssiSymbolIdx);
        int compareSuffixWithPrefix(uint_reads_cnt sufIncIdx, uint_reads_cnt preIncIdx, uint_read_len sufOffset);

        virtual uint_read_len readLength(uint_reads_cnt incIdx) override;
        virtual string getReadUpToOverlap(uint_reads_cnt incIdx) override;
        virtual uint_reads_cnt readsTotal() override;

        virtual ReadsSetProperties* getReadsSetProperties() override;

        template<bool pgGenerationMode>
        void initAndFindDuplicates();
        void prepareSortedReadsBlocks();
        void mergeSortOfLeftSuffixes(uint8_t offset, const uint_reads_cnt *sortedSuffixesLeftCount,
                uint_reads_cnt *sortedSuffixLeftIdxsPtr, const uint_reads_cnt *sortedSuffixIdxsPtr);

        template<bool pgGenerationMode>
        void overlapSortedReadsAndSuffixes(uint_read_len suffixesOffset, uint_reads_cnt *sortedSuffixesLeftCount);

        void validateSortedSuffixes(uint8_t offset) const;

        void findOverlappingReads(double overlappedReadsCountStopCoef, bool pgGenerationMode) override;

    protected:
        template<bool pgGenerationMode>
        void setDuplicateSuccessor(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght);

    public:

        ParallelGreedySwipingPackedOverlapGeneratorTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet = false);
        virtual ~ParallelGreedySwipingPackedOverlapGeneratorTemplate();

        bool isPseudoGenomeLengthStandardVirtual();
        bool isPseudoGenomeLengthMaximalVirtual();
    };

    class ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory: public PseudoGenomeGeneratorFactory
    {
    private:

        template<typename uint_read_len, typename uint_reads_cnt>
        PseudoGenomeGeneratorBase* getGeneratorFullTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet);

        template<typename uint_read_len>
        PseudoGenomeGeneratorBase* getGeneratorPartialTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet);

    public:

        ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory() {};

        PseudoGenomeGeneratorBase* getGenerator(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator);
        PseudoGenomeGeneratorBase* getGenerator(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet);

        static PseudoGenomeBase* generatePg(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator);
        static PseudoGenomeBase* generatePg(PackedConstantLengthReadsSet *readsSet);
        static SeparatedPseudoGenome* generateSeparatedPg(PackedConstantLengthReadsSet *readsSet);
        static const vector<bool> getHQReads(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator,
                                             double qualityCoef);
        static const vector<bool> getHQReads(PackedConstantLengthReadsSet *readsSet,
                                             double qualityCoef);
    };

}


#endif // PARALLELGREEDYSWIPINGPACKEDOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
