#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <set>
#include <algorithm>
#include <climits>
#include <cmath>
#include <sstream>

using namespace std;

namespace PgSAHelpers {

    void clock_checkpoint();
    unsigned long long int clock_millis();

    void* readArray(std::istream&, size_t arraySizeInBytes);
    void writeArray(std::ostream&, void* srcArray, size_t arraySize);

    void* readWholeArray(std::istream&, size_t& arraySizeInBytes);
    void* readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes);
    void writeArrayToFile(string destFile, void* srcArray, size_t arraySize);

    string toString(unsigned long long value);
    string toMB(unsigned long long value, unsigned char decimalPlaces);
    string toString(long double value, unsigned char decimalPlaces);
   
    unsigned long long int powuint(unsigned long long int base, int exp);
    
    string reverseComplement(string kmer);
    
    template<typename uint>
    inline uint divideBySmallInteger(const uint dividend, const unsigned char divisor) {
        switch (divisor) {
            case 1: return dividend / 1;
            case 2: return dividend / 2;
            case 3: return dividend / 3;
            case 4: return dividend / 4;
            case 5: return dividend / 5;
            case 6: return dividend / 6;
            case 7: return dividend / 7;
            case 8: return dividend / 8;
        };
        cout << "Unsupported denominator " << divisor << "!\n";
        return 0;
    }
    
    template<typename uint>
    inline uint moduloBySmallInteger(const uint dividend, const unsigned char divisor, const uint resultOfDivision) {
        return dividend - resultOfDivision * divisor;
    }
    
    template<typename uint>
    inline uint moduloBySmallInteger(const uint dividend, const unsigned char divisor) {
        return moduloBySmallInteger(dividend, divisor, divideBySmallInteger(dividend, divisor));
    }
    
    class BufferedFileIStream : public istream {

        private:

            struct membuf : std::streambuf
            {
                membuf(char* begin, char* end) {
                    this->setg(begin, begin, end);
                }
            };

            membuf* sbuf = 0;
            char* readsArray = 0;

            BufferedFileIStream(membuf* sbuf, char* readsArray) : istream(sbuf) {
                this->sbuf = sbuf;
                this->readsArray = readsArray;
            }

        public:

            static BufferedFileIStream* getIStream(string filename) {
                size_t readsArraySize;
                char* readsArray = (char*) PgSAHelpers::readWholeArrayFromFile(filename, readsArraySize);

                membuf* sbuf = new membuf(readsArray, readsArray + readsArraySize);

                BufferedFileIStream* source = new BufferedFileIStream(sbuf, readsArray);
                if (!*source)
                    std::cout << "Whoops";

                return source;
            }

            virtual ~BufferedFileIStream() {
                delete(sbuf);
                delete[] readsArray;
            }
    };

    int readsSufPreCmp(const char* suffixPart, const char* prefixRead);
    
    int strcmplcp(const char* lStrPtr, const char* rStrPtr, int length);
    
}

#endif // HELPER_H_INCLUDED
