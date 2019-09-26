/* 
 * Author: Tomek
 *
 * Created on 13 luty 2014, 10:27
 */

#ifndef PACKEDPSEUDOGENOMEGENERATOR_H_INCLUDED
#define	PACKEDPSEUDOGENOMEGENERATOR_H_INCLUDED

#include "PseudoGenomeGeneratorBase.h"
#include "../PackedPseudoGenome.h"

namespace PgSAIndex {

    class PackedPseudoGenomeGenerator: public PseudoGenomeGeneratorBase {
        private:
            uint_pg_len_max pseudoGenomeLength;
            PseudoGenomeBase* pgb; 
            uchar symbolsPerElement;
            
        public:
            
            PackedPseudoGenomeGenerator(PseudoGenomeBase* pgb, uchar symbolsPerElement)
            : pseudoGenomeLength(pgb->getPseudoGenomeLength()),
                    pgb(pgb), 
                    symbolsPerElement(symbolsPerElement) 
            { };
            
            ~PackedPseudoGenomeGenerator() {} ;
                        
            bool isPseudoGenomeLengthStandardVirtual() { return isPGLengthStd(pseudoGenomeLength); };
            bool isPseudoGenomeLengthMaximalVirtual() { return isPGLengthMax(pseudoGenomeLength); };
            
            PseudoGenomeBase* generatePseudoGenomeBase() {
                if (pgb->getTypeID() == PGTYPE_DEFAULT) {
                    if (pgb->isReadLengthMin()) {
                        if (pgb->isReadsCountStd()) {
                            if (pgb->isPGLengthStd()) {
                                return generatePseudoGenomeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>();
                            }
                            if (pgb->isPGLengthMax())
                                return generatePseudoGenomeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>();
                        }
                    }
                    if (pgb->isReadLengthStd()) {
                        if (pgb->isReadsCountStd()) {
                            if (pgb->isPGLengthStd())
                                return generatePseudoGenomeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>();
                            if (pgb->isPGLengthMax())
                                return generatePseudoGenomeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>();
                        }                            
                    }
                }
                cout << "ERROR: wrong source PGSATYPE " << pgb->getTypeID();
                return 0;
            };
     
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            PseudoGenomeBase* generatePseudoGenomeTemplate() {
                uchar symbolsCount = pgb->getReadsSetProperties()->symbolsCount;
                
                if (pgb->isReadLengthConstant()) {
                    typedef DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len> DefaultPseudoGenomeClass;
                    DefaultPseudoGenomeClass* pg = DefaultPseudoGenomeClass::castBase(pgb);
                    if(SymbolsPackingFacility::isCompatible(symbolsPerElement, symbolsCount)) {
                        typedef PackedPseudoGenomeOfConstantLengthReadsType <uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min> PackedPseudoGenomeClass;
                        return new PackedPseudoGenomeClass(pg, symbolsPerElement);
                    }
                }
                cout << "ERROR: wrong source PGSATYPE " << pgb->getTypeID();
                return 0;
            }
            
    };

}

#endif	/* PACKEDPSEUDOGENOMEGENERATOR_H_INCLUDED */

