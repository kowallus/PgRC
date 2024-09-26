#ifndef PSEUDOGENOMEGENERATORBASE_H_INCLUDED
#define PSEUDOGENOMEGENERATORBASE_H_INCLUDED

#include "../../pseudogenome/DefaultPseudoGenome.h"
#include "../../pseudogenome/PseudoGenomeInterface.h"
#include "../../readsset/DefaultReadsSet.h"
#include "../SeparatedPseudoGenome.h"

using namespace PgHelpers;
using namespace PgReadsSet;
using namespace PgTools;

namespace PgIndex {

    class PseudoGenomeGeneratorBase
    {
        public:

            virtual SeparatedPseudoGenome* generateSeparatedPseudoGenome() = 0;
            virtual PseudoGenomeBase* generatePseudoGenomeBase() = 0;
            virtual const vector<bool> getBothSidesOverlappedReads(double overlappedReadsCountStopCoef) = 0;


            virtual ~PseudoGenomeGeneratorBase() {} ;

    };

    class PseudoGenomeGeneratorFactory
    {
        public:
            virtual ~PseudoGenomeGeneratorFactory() {};
            
            virtual PseudoGenomeGeneratorBase* getGenerator(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator) = 0;

    };
    
}


#endif // PSEUDOGENOMEGENERATORBASE_H_INCLUDED
