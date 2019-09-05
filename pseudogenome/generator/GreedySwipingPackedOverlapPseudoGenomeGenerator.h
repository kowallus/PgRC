#ifndef GREEDYSWIPINGPACKEDOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
#define GREEDYSWIPINGPACKEDOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include "AbstractOverlapPseudoGenomeGenerator.h"
#include "../../readsset/PackedConstantLengthReadsSet.h"
#include <algorithm>
#include <deque>

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt >
    class GreedySwipingPackedOverlapGeneratorTemplate: public AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>
    {
    private:

        bool ownReadsSet = false;
        PackedConstantLengthReadsSet* packedReadsSet = 0;

        vector<uint_reads_cnt> sortedReadsIdxs;

        struct PackedReadsComparator {
            GreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>* myGenerator;
    PackedReadsComparator(GreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>* generator) : myGenerator(generator) {};

            bool operator() (uint_reads_cnt lIdx, uint_reads_cnt rIdx) const {
                return myGenerator->compareReads(lIdx, rIdx) < 0;
            }
        };

        int compareReads(uint_reads_cnt lIncIdx, uint_reads_cnt rIncIdx);

        uchar getSymbolOrderFromRead(uint_reads_cnt incIdx, uint_read_len offset);
        vector<uint_reads_cnt> sortedSuffixIdxs;
        vector<uint_reads_cnt> ssiSymbolIdx {vector<uint_reads_cnt>(UCHAR_MAX, 0)};
        vector<uint_reads_cnt> ssiSymbolEnd {vector<uint_reads_cnt>(UCHAR_MAX, 0)};
        deque<uchar> ssiOrder;
        void updateSuffixQueue(uchar suffixGroup, uint_read_len suffixOffset);

        int compareSuffixes(uchar lSymOrder, uchar rSymOrder, uint_read_len offset);
        int compareSuffixWithPrefix(uint_reads_cnt sufIncIdx, uint_reads_cnt preIncIdx, uint_read_len sufOffset);

        virtual uint_read_len readLength(uint_reads_cnt incIdx) override;
        virtual string getReadUpToOverlap(uint_reads_cnt incIdx) override;
        virtual uint_reads_cnt readsTotal() override;

        virtual ReadsSetProperties* getReadsSetProperties() override;

        bool isGenerationCyclesAware(bool pgGenerationMode) { return false; };

        template<bool pgGenerationMode>
        void initAndFindDuplicates();
        template<bool pgGenerationMode>
        void overlapSortedReadsAndMergeSortSuffixes(uint_read_len suffixesOffset);

        void findOverlappingReads(double overlappedReadsCountStopCoef, bool pgGenerationMode) override;

    public:

        GreedySwipingPackedOverlapGeneratorTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet = false);
        virtual ~GreedySwipingPackedOverlapGeneratorTemplate();

        bool isPseudoGenomeLengthStandardVirtual();
        bool isPseudoGenomeLengthMaximalVirtual();

    };

    class GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory: public PseudoGenomeGeneratorFactory
    {
    private:

        template<typename uint_read_len, typename uint_reads_cnt>
        PseudoGenomeGeneratorBase* getGeneratorFullTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet);

        template<typename uint_read_len>
        PseudoGenomeGeneratorBase* getGeneratorPartialTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet);

    public:

        GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory() {};

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


#endif // GREEDYSWIPINGPACKEDOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
