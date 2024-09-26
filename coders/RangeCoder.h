#ifndef PGTOOLS_RANGECODER_H
#define PGTOOLS_RANGECODER_H

#include "LzmaCoder.h"

class RangeCoderProps: public CoderProps {
private:

    uint16_t normalizeNoOfSymbols(uint16_t no_of_symbols) {
        if (no_of_symbols > 256) {
            cout << "Unsupported number of symbols " << no_of_symbols << " in range coder (only up to 256)." << endl;
            exit(EXIT_FAILURE);
        }
        if (no_of_symbols > 128)
            return 256;
        else if (no_of_symbols > 64)
            return 128;
        else if (no_of_symbols > 8)
            return 64;
        else if (no_of_symbols > 6)
            return 8;
        else if (no_of_symbols > 4)
            return 6;
        else if (no_of_symbols > 2)
            return 4;
        else
            return 2;
    }

public:
    uint16_t no_of_symbols;
    uint8_t schema;
    uint8_t bytes_period;

    const static uint8_t PERIODIC_SCHEMA = 0;

    RangeCoderProps(uint16_t no_of_symbols, uint8_t bytes_period = 1)
            : CoderProps(RANGE_CODER, 1), no_of_symbols(normalizeNoOfSymbols(no_of_symbols)),
            schema(PERIODIC_SCHEMA), bytes_period(bytes_period)  { }

    string log() {
        return " range_coder (sigma = " + to_string(no_of_symbols) + "; period = " + to_string(bytes_period) + ")";
    }

    virtual ~RangeCoderProps() { };
};

MY_STDAPI RangeCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                        RangeCoderProps* coder_props, double estimated_compression);

MY_STDAPI RangeUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen,
                          ostream* logout = PgHelpers::devout);

#endif //PGTOOLS_RANGECODER_H
