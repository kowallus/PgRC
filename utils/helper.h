#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

#include <ctime>
#include <chrono>
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

class NullBuffer : public std::streambuf
{
public:
    int overflow(int c) { return c; }
};

extern std::ostream null_stream;

namespace PgSAHelpers {

    extern int numberOfThreads;

    // bioinformatical routines

    char reverseComplement(char symbol);
    string reverseComplement(string kmer);
    void reverseComplementInPlace(char* start, const std::size_t N);
    void reverseComplementInPlace(string &kmer);
    double qualityScore2approxCorrectProb(string quality);
    double qualityScore2correctProb(string quality);

    inline uint8_t symbol2value(char symbol);
    inline char value2symbol(uint8_t value);
    uint8_t mismatch2code(char actual, char mismatch);
    char code2mismatch(char actual, uint8_t code);

    template<typename uint_read_len>
    void convertMisRevOffsets2Offsets(uint_read_len *mismatchOffsets, uint8_t mismatchesCount, uint_read_len readLength) {
        for(uint8_t i = 0; i < mismatchesCount / 2; i++) {
            uint_read_len tmp = mismatchOffsets[mismatchesCount - i - 1];
            mismatchOffsets[mismatchesCount - i - 1] = mismatchOffsets[i];
            mismatchOffsets[i] = tmp;
        }
        uint_read_len pos = readLength;
        for(uint8_t i = mismatchesCount; i-- > 0;) {
            pos -= mismatchOffsets[i] + 1;
            mismatchOffsets[i] = pos;
        }
    }

    // time routines

    clock_t clock_checkpoint();
    unsigned long long int clock_millis();
    unsigned long long int clock_millis(clock_t checkpoint);

    chrono::steady_clock::time_point time_checkpoint();
    unsigned long long int time_millis();
    unsigned long long int time_millis(chrono::steady_clock::time_point checkpoint);


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
    inline uint multiplyBySmallInteger(const uint value, const unsigned char smallInteger) {
        switch (smallInteger) {
            case 0: return 0;
            case 1: return value * 1;
            case 2: return value * 2;
            case 3: return value * 3;
            case 4: return value * 4;
            case 5: return value * 5;
            case 6: return value * 6;
            case 7: return value * 7;
            case 8: return value * 8;
        };
        cout << "Unsupported multiplied " << smallInteger << "!\n";
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

    template<typename t_val>
    string transpose(string matrix, uint64_t rows, uint64_t cols) {
        string tRes;
        tRes.resize(matrix.size());
        t_val* nMatrix = (t_val*) matrix.data();
        t_val* tMatrix = (t_val*) tRes.data();
        for(uint64_t i = 0; i < rows; i++) {
            for(uint8_t m = 0; m < cols; m++)
                tMatrix[m * rows + i] = nMatrix[i * cols + m];
        }
        return tRes;
    }

    // string comparison routines

    int readsSufPreCmp(const char* suffixPart, const char* prefixRead);
    
    int strcmplcp(const char* lStrPtr, const char* rStrPtr, int length);

    // input output routines

    extern std::ostream *logout;

    void* readArray(std::istream&, size_t arraySizeInBytes);
    void readArray(std::istream&, void* destArray, size_t arraySizeInBytes);
    void writeArray(std::ostream&, void* srcArray, size_t arraySize, bool verbose = false);

    void* readWholeArray(std::istream&, size_t& arraySizeInBytes);
    void* readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes);
    void writeArrayToFile(string destFile, void* srcArray, size_t arraySize);
    void writeStringToFile(string destFile, const string &src);

    extern bool plainTextWriteMode;
    const static string TEXT_MODE_ID = "TXT";
    const static string BINARY_MODE_ID = "BIN";

    void writeReadMode(std::ostream &dest, bool plainTextWriteMode);
    bool confirmTextReadMode(std::istream &src);

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

    void writeUIntByteFrugal(std::ostream &dest, uint64_t value);

    template<typename t_val>
    void readUIntByteFrugal(std::istream &src, t_val& value) {
        value = 0;
        uint8_t yByte = 0;
        t_val base = 1;
        do {
            src.read((char *) &yByte, sizeof(uint8_t));
            value += base * (yByte % 128);
            base *= 128;
        } while (yByte >= 128);
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
