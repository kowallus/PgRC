#ifndef ABSTRACTOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
#define ABSTRACTOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include <algorithm>
#include <set>
#include "../persistence/SeparatedPseudoGenomePersistence.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt >
    class AbstractOverlapPseudoGenomeGeneratorTemplate: public PseudoGenomeGeneratorBase
    {

    protected:

        // auxiliary structures
        uint_reads_cnt* nextRead = 0;
        uint_read_len* overlap = 0;
        uint_reads_cnt* headRead = 0;
        uint_reads_cnt readsLeft;

        bool hasPredecessor(uint_reads_cnt incIdx);
        bool hasSuccessor(uint_reads_cnt incIdx);
        void unionOverlappedReads(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght);
        void setReadSuccessor(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght);

        uint_reads_cnt getHead(uint_reads_cnt idx);
        bool isHeadOf(uint_reads_cnt head, uint_reads_cnt idx);

        uint_pg_len_max pseudoGenomeLength;

        virtual uint_read_len readLength(uint_reads_cnt incIdx) = 0;
        virtual string getReadUpToOverlap(uint_reads_cnt incIdx) = 0;
        virtual uint_reads_cnt readsTotal() = 0;

        virtual ReadsSetProperties* getReadsSetProperties() = 0;

        void performOverlapping(double overlappedReadsCountStopCoef = 1, bool pgGenerationMode = true);

        template<class GeneratedPseudoGenome>
        GeneratedPseudoGenome* assemblePseudoGenomeTemplate();

        virtual void findOverlappingReads(double overlappedReadsCountStopCoef, bool pgGenerationMode) = 0;

        uint_pg_len_max countPseudoGenomeLength();
        uint_reads_cnt countSingles();
        uint_reads_cnt countComponents();
        void quick_stats();

        void init(bool pgGenerationMode = true);
        void dispose();

        virtual bool isGenerationCyclesAware(bool pgGenerationMode) = 0;
        void removeCyclesAndPrepareComponents();

    public:

        AbstractOverlapPseudoGenomeGeneratorTemplate() {}
        virtual ~AbstractOverlapPseudoGenomeGeneratorTemplate() {}

        SeparatedPseudoGenome *generateSeparatedPseudoGenome() override;

        PseudoGenomeBase *generatePseudoGenomeBase() override;

        const vector<bool> getBothSidesOverlappedReads(double overlappedReadsCountStopCoef) override;

        bool isPseudoGenomeLengthStandardVirtual();
        bool isPseudoGenomeLengthMaximalVirtual();
    };

}


#endif // ABSTRACTOVERLAPPSEUDOGENOMEGENERATOR_H_INCLUDED
