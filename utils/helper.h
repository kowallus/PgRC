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

    // bioinformatical routines

    string reverseComplement(string kmer);
    void reverseComplementInPlace(string &kmer);
    double qualityScore2approxCorrectProb(string quality);
    double qualityScore2correctProb(string quality);

    inline uint8_t symbol2value(char symbol);
    inline char value2symbol(uint8_t value);
    uint8_t mismatch2code(char actual, char mismatch);
    char code2mismatch(char actual, uint8_t code);

    void convertMisOffsets2RevOffsets(uint16_t* mismatchOffsets, uint8_t mismatchesCount, uint16_t readLength);

    // time routines

    void clock_checkpoint();
    unsigned long long int clock_millis();

    // string conversion routines

    string toString(unsigned long long value);
    string toMB(unsigned long long value, unsigned char decimalPlaces);
    string toString(long double value, unsigned char decimalPlaces);


    // mathematical routines

    unsigned long long int powuint(unsigned long long int base, int exp);

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

    // string comparison routines

    int readsSufPreCmp(const char* suffixPart, const char* prefixRead);
    
    int strcmplcp(const char* lStrPtr, const char* rStrPtr, int length);

    // input output routines

    void* readArray(std::istream&, size_t arraySizeInBytes);
    void readArray(std::istream&, void* destArray, size_t arraySizeInBytes);
    void writeArray(std::ostream&, void* srcArray, size_t arraySize);

    void* readWholeArray(std::istream&, size_t& arraySizeInBytes);
    void* readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes);
    void writeArrayToFile(string destFile, void* srcArray, size_t arraySize);

    extern bool plainTextWriteMode;
    const static string TEXT_MODE_ID = "TXT";
    const static string BINARY_MODE_ID = "BIN";

    void writeReadMode(std::ostream &dest, bool plainTextWriteMode);
    bool readReadMode(std::istream &src);

    template<typename t_val>
    void writeValue(std::ostream &dest, const t_val value, bool plainTextWriteMode) {
        if (plainTextWriteMode)
            dest << value << endl;
        else
            dest.write((char *) &value, sizeof(t_val));
    }
    template<typename t_val>
    void readValue(std::istream &src, t_val& value, bool plainTextReadMode) {
        if (plainTextReadMode)
            src >> value;
        else
            src.read((char *) &value, sizeof(t_val));
    }

    template<>
    void writeValue(std::ostream &dest, const uint8_t value, bool plainTextWriteMode);
    template<>
    void readValue(std::istream &src, uint8_t& value, bool plainTextReadMode);

    template<typename t_val>
    void writeValue(std::ostream &dest, const t_val value) {
        writeValue(dest, value, plainTextWriteMode);
    }

    extern bool bytePerReadLengthMode;

    void readReadLengthValue(std::istream &src, uint16_t& value, bool plainTextReadMode);

    void writeReadLengthValue(std::ostream &dest, const uint16_t value);

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

}

#endif // HELPER_H_INCLUDED
