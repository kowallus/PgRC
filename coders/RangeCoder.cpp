/* RangeCoder library wrapper */

#include "RangeCoder.h"

#include "rangecoder/clr.cdr"
#include "rangecoder/simple_model.h"
#include "rangecoder/base_model.h"
#include "lzma/CpuArch.h"

#include "CodersLib.h"
#include <assert.h>

const static uint8_t RANGE_CODER_PROPS_SIZE = 4;

template<typename MODEL_T>
MY_STDAPI RangeCoderCompressTemplate(unsigned char *dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                                     int period) {

    char * output = (char*) dest;

    RangeCoder rc{};
    rc.output(output);
    rc.StartEncode();

    MODEL_T value;

    if (period < 2) {
        for (size_t i = 0; i < srcLen; ++i) {
            value.encodeSymbol(&rc, src[i]);
        }
    } else {
        MODEL_T first_value;
        size_t i = 0;
        while (i < srcLen) {
            first_value.encodeSymbol(&rc, src[i++]);
            for (size_t m = 1; m < period; ++m) {
                value.encodeSymbol(&rc, src[i++]);
            }
        }
    }
    rc.FinishEncode();
    destLen += rc.size_out();

    return SZ_OK;
}

MY_STDAPI RangeCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                        RangeCoderProps* props, double estimated_compression) {

    if (props->getNumThreads() != 1) {
        cout << "Unsupported number of threads " << props->getNumThreads() << " in range coder." << endl;
        exit(EXIT_FAILURE);
    }
    size_t propsSize = RANGE_CODER_PROPS_SIZE;
    size_t maxDestSize = propsSize + (srcLen + srcLen / 3) * estimated_compression + 1024;
    try {
        dest = new unsigned char[maxDestSize];
    } catch (const std::bad_alloc& e) {
        cout << "Allocation failed: " << e.what() << endl;
        maxDestSize -= srcLen / 3 * estimated_compression;
        dest = new unsigned char[maxDestSize];
    }
    *dest = (unsigned char) props->schema;
    *(dest + 1) = (unsigned char) props->bytes_period;
    SetUi16(dest + 2, props->no_of_symbols);
    destLen = propsSize;

    MY_STDAPI res;
    switch (props->no_of_symbols) {
        case 256: res = RangeCoderCompressTemplate<SIMPLE_MODEL<256>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        case 128: res = RangeCoderCompressTemplate<SIMPLE_MODEL<128>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        case 64: res = RangeCoderCompressTemplate<SIMPLE_MODEL<64>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        case 8: res = RangeCoderCompressTemplate<SIMPLE_MODEL<8>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        case 6: res = RangeCoderCompressTemplate<BASE_MODEL<uint8_t>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        case 4: res = RangeCoderCompressTemplate<SIMPLE_MODEL<4>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        case 2: res = RangeCoderCompressTemplate<SIMPLE_MODEL<2>>(dest + propsSize, destLen, src, srcLen, props->bytes_period);
            break;
        default:
            cout << "Unsupported number of symbols " << props->no_of_symbols << " in range coder." << endl;
            exit(EXIT_FAILURE);
    }

    return res;
}

template<typename MODEL_T>
MY_STDAPI RangeCoderDecompressTemplate(unsigned char *dest, size_t destLen, const unsigned char *src, size_t srcLen,
                                     int period) {

    RangeCoder rc{};
    rc.input((char*) src);
    rc.StartDecode();

    MODEL_T value;

    if (period < 2) {
        for (size_t i = 0; i < destLen; ++i) {
            dest[i] = value.decodeSymbol(&rc);
        }
    } else {
        MODEL_T first_value;
        size_t i = 0;
        while (i < destLen) {
            dest[i++] = first_value.decodeSymbol(&rc);
            for (size_t m = 1; m < period; ++m) {
                dest[i++] = value.decodeSymbol(&rc);
            }
        }
        assert(i == destLen);
    }
    rc.FinishDecode();

    return SZ_OK;
}

MY_STDAPI RangeUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen, ostream* logout) {
    size_t propsSize = RANGE_CODER_PROPS_SIZE;
    unsigned char propsBuf[RANGE_CODER_PROPS_SIZE];
    memcpy(propsBuf, src, propsSize);
    src += propsSize;
    uint8_t schema = propsBuf[0];
    if (schema != RangeCoderProps::PERIODIC_SCHEMA) {
        cout << "Unsupported schema " << schema << " in range coder." << endl;
        exit(EXIT_FAILURE);
    }
    uint8_t bytes_period = propsBuf[1];
    uint16_t no_of_symbols = GetUi16(propsBuf + 2);
    srcLen -= propsSize;

    *logout << "... range_coder (sigma = " + to_string(no_of_symbols) + "; period = " + to_string(bytes_period) + ") ... ";

    MY_STDAPI res;
    switch (no_of_symbols) {
        case 256: res = RangeCoderDecompressTemplate<SIMPLE_MODEL<256>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        case 128: res = RangeCoderDecompressTemplate<SIMPLE_MODEL<128>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        case 64: res = RangeCoderDecompressTemplate<SIMPLE_MODEL<64>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        case 8: res = RangeCoderDecompressTemplate<SIMPLE_MODEL<8>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        case 6: res = RangeCoderDecompressTemplate<BASE_MODEL<uint8_t>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        case 4: res = RangeCoderDecompressTemplate<SIMPLE_MODEL<4>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        case 2: res = RangeCoderDecompressTemplate<SIMPLE_MODEL<2>>(dest, *destLen, src, srcLen, bytes_period);
            break;
        default:
            cout << "Unsupported number of symbols " << no_of_symbols << " in range coder decompressor." << endl;
            exit(EXIT_FAILURE);
    }

    return res;
}
