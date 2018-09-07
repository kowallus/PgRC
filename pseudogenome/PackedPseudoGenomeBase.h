#ifndef PACKEDPSEUDOGENOMEBASE_H_INCLUDED
#define PACKEDPSEUDOGENOMEBASE_H_INCLUDED

#include "PseudoGenomeBase.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    class PackedPseudoGenomeBase: public PseudoGenomeBase
    {
        protected:
            uchar symbolsPerElement; 
            uchar bytesPerElement;
            
            PackedPseudoGenomeBase(uint_pg_len_max length, ReadsSetProperties* properties, uchar symbolsPerElement, uchar bytesPerElement)
            : PseudoGenomeBase(length, properties),
              symbolsPerElement(symbolsPerElement),
              bytesPerElement(bytesPerElement)
            { };

            PackedPseudoGenomeBase(uint_pg_len_max length, std::istream& src, uchar bytesPerElement)
            : PseudoGenomeBase(length, src){ 
                int srchelper;
                src >> srchelper;
                this->symbolsPerElement = srchelper;
                
                this->bytesPerElement = bytesPerElement;
            };

            virtual ~PackedPseudoGenomeBase() { };

        public:

            const uchar getSymbolsPerElement() { return this->symbolsPerElement; };
            const uchar getBytesPerElement() { return this->bytesPerElement; };
      
            bool isPgElementMinimal() { return bytesPerElement == 1; };
            bool isPgElementStandard() { return bytesPerElement == 2; };
    };

    class PackedPseudoGenomeHeaderExtension {
        private:
            uchar symbolsPerElement;
            uchar bytesPerElement;

        public:

            PackedPseudoGenomeHeaderExtension(PackedPseudoGenomeBase* ppgb) {
                this->symbolsPerElement = ppgb->getSymbolsPerElement();
                this->bytesPerElement = ppgb->getBytesPerElement();
            }

            PackedPseudoGenomeHeaderExtension(std::istream& src) {
                int srchelper;
                src >> srchelper;
                symbolsPerElement = srchelper;
                src >> srchelper;
                bytesPerElement = srchelper;                
                src.get();
            }

            ~PackedPseudoGenomeHeaderExtension() { };
            
            void write(std::ostream& dest) {
                dest << (int) symbolsPerElement << "\n";
                dest << (int) bytesPerElement << "\n";
            }

            bool isPgElementMinimal() { return bytesPerElement == 1; };
            bool isPgElementStandard() { return bytesPerElement == 2; };
    };
}

#endif // PACKEDPSEUDOGENOMEBASE_H_INCLUDED
