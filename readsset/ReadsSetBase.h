#ifndef READSSETBASE_H_INCLUDED
#define READSSETBASE_H_INCLUDED

#include "../utils/helper.h"
#include "../pseudogenome/readslist/ReadsListInterface.h"
#include "iterator/ReadsSetIterator.h"

namespace PgSAReadsSet {

    typedef unsigned long long uint_reads_total_len;

    class ReadsSetProperties {

        public:

            uint_reads_cnt_max readsCount = 0;
            uint_reads_total_len allReadsLength = 0;

            bool constantReadLength = true;
            uint_read_len_max maxReadLength = 0;

            uint_symbols_cnt symbolsCount = 0;
            char symbolsList[UCHAR_MAX] = {0};
            int symbolOrder[UCHAR_MAX] = {-1};

            ReadsSetProperties() {};

            ReadsSetProperties(std::istream& source) {
                source >> readsCount;
                source >> allReadsLength;
                source >> constantReadLength;
                source >> maxReadLength;
                int srchelper;
                source >> srchelper;
                symbolsCount = srchelper;
                source >> symbolsList;
                generateSymbolOrder();
            }

            void copy(ReadsSetProperties* properties) {
                readsCount = properties->readsCount;
                allReadsLength = properties->allReadsLength;

                constantReadLength = properties->constantReadLength;
                maxReadLength = properties->maxReadLength;

                symbolsCount = properties->symbolsCount;
                std::copy(std::begin(properties->symbolsList), std::end(properties->symbolsList), std::begin(symbolsList));
                std::copy(std::begin(properties->symbolOrder), std::end(properties->symbolOrder), std::begin(symbolOrder));
            }
            
            bool compareWith(ReadsSetProperties* properties) {
                return (readsCount == properties->readsCount &&
                        allReadsLength == properties->allReadsLength &&
                        constantReadLength == properties->constantReadLength &&
                        maxReadLength == properties->maxReadLength &&
                        symbolsCount == properties->symbolsCount &&
                        std::equal(std::begin(properties->symbolsList), std::end(properties->symbolsList), std::begin(symbolsList)) &&
                        std::equal(std::begin(properties->symbolOrder), std::end(properties->symbolOrder), std::begin(symbolOrder)));
            }

            void write(std::ostream& dest) {
                dest << readsCount << "\n"
                        << allReadsLength << "\n"
                        << constantReadLength << "\n"
                        << maxReadLength << "\n"
                        << (int) symbolsCount << "\n"
                        << symbolsList << "\n";
            }

            void generateSymbolOrder() {
                for (uint_symbols_cnt i = 0; i < symbolsCount; i++)
                    symbolOrder[(unsigned char) symbolsList[(unsigned char) i]] = i;
            }

            void printout() {
                std::cout << "reads count: " << readsCount << "\n";
                std::cout << "all reads length: " << allReadsLength << "\n";
                std::cout << "reads length is " << (constantReadLength?"constant":"variable") << "\n";
                std::cout << "maxReadLength: " << maxReadLength << "\n";
                std::cout << "symbolsCount: " << (int) symbolsCount << "\n";
                std::cout << "symbols: ";
                for(uint_symbols_cnt i = 0; i < symbolsCount; i++)
                    std::cout << symbolsList[(unsigned char) i];
                std::cout << "\n\n";
            }
    };

    class ReadsSetBase
    {
        protected:

            ReadsSetProperties* properties = 0;

            ReadsSetBase() {
                this->properties = new ReadsSetProperties();
            }

            ReadsSetBase(ReadsSetProperties* properties)
            : ReadsSetBase() {
                this->properties->copy(properties);
            };

            ReadsSetBase(std::istream& src) {
                this->properties = new ReadsSetProperties(src);
            };

            virtual ~ReadsSetBase() { 
                delete(properties); 
            };

        public:

            ReadsSetProperties* getReadsSetProperties() { return this->properties; };

            bool isReadLengthConstant() { return properties->constantReadLength; };

            bool isReadLengthMin() { return PgSAReadsSet::isReadLengthMin(properties->maxReadLength); };
            bool isReadLengthStd() { return PgSAReadsSet::isReadLengthStd(properties->maxReadLength); };

            bool isReadsCountStd() { return PgSAReadsSet::isReadsCountStd(properties->readsCount); };
            bool isReadsCountMax() { return PgSAReadsSet::isReadsCountMax(properties->readsCount); };

    };

}

#endif // READSSETBASE_H_INCLUDED
