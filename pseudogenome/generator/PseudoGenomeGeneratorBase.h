#ifndef PSEUDOGENOMEGENERATORBASE_H_INCLUDED
#define PSEUDOGENOMEGENERATORBASE_H_INCLUDED

#include "../../pseudogenome/DefaultPseudoGenome.h"
#include "../../pseudogenome/PseudoGenomeInterface.h"
#include "../../readsset/DefaultReadsSet.h"

using namespace PgSAHelpers;
using namespace PgSAReadsSet;

namespace PgSAIndex {

    class PseudoGenomeGeneratorBase
    {
        public:

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
