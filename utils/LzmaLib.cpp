/* LzmaLib.c -- LZMA library wrapper
2015-06-13 : Igor Pavlov : Public domain */

#include <assert.h>
#include "../lzma/Alloc.h"
#include "../lzma/LzmaDec.h"
#include "../lzma/LzmaEnc.h"
#include "LzmaLib.h"

/*
RAM requirements for LZMA:
  for compression:   (dictSize * 11.5 + 6 MB) + state_size
  for decompression: dictSize + state_size
    state_size = (4 + (1.5 << (lc + lp))) KB
    by default (lc=3, lp=0), state_size = 16 KB.

LZMA properties (5 bytes) format
    Offset Size  Description
      0     1    lc, lp and pb in encoded form.
      1     4    dictSize (little endian).
*/

/*
LzmaProperties
------------

  LZMA Encoder will use default values for any parameter, if it is
  -1  for any from: level, loc, lp, pb, fb, numThreads
   0  for dictSize

level - compression level: 0 <= level <= 9;

  level dictSize algo  fb
    0:    16 KB   0    32
    1:    64 KB   0    32
    2:   256 KB   0    32
    3:     1 MB   0    32
    4:     4 MB   0    32
    5:    16 MB   1    32
    6:    32 MB   1    32
    7+:   64 MB   1    64

  The default value for "level" is 5.

  algo = 0 means fast method
  algo = 1 means normal method

dictSize - The dictionary size in bytes. The maximum value is
        128 MB = (1 << 27) bytes for 32-bit version
          1 GB = (1 << 30) bytes for 64-bit version
     The default value is 16 MB = (1 << 24) bytes.
     It's recommended to use the dictionary that is larger than 4 KB and
     that can be calculated as (1 << N) or (3 << N) sizes.

lc - The number of literal context bits (high bits of previous literal).
     It can be in the range from 0 to 8. The default value is 3.
     Sometimes lc=4 gives the gain for big files.

lp - The number of literal pos bits (low bits of current position for literals).
     It can be in the range from 0 to 4. The default value is 0.
     The lp switch is intended for periodical data when the period is equal to 2^lp.
     For example, for 32-bit (4 bytes) periodical data you can use lp=2. Often it's
     better to set lc=0, if you change lp switch.

pb - The number of pos bits (low bits of current position).
     It can be in the range from 0 to 4. The default value is 2.
     The pb switch is intended for periodical data when the period is equal 2^pb.

fb - Word size (the number of fast bytes).
     It can be in the range from 5 to 273. The default value is 32.
     Usually, a big number gives a little bit better compression ratio and
     slower compression process.
*/

void LzmaEncProps_Set(CLzmaEncProps *p, int coder_level, size_t dataLength, int numThreads, 
        int dataPeriodCode = -1) {
    switch(coder_level) {
        case PGRC_CODER_LEVEL_NORMAL:
            p->level = -1;
            p->dictSize = 0;
            p->lc = -1;
            p->lp = -1;
            p->pb = -1;
            p->fb = -1;
            break;
        case PGRC_CODER_LEVEL_MAXIMUM:
            p->level = 9;
            p->dictSize = 3 << 29;
            p->lc = 3; // test 4 for big files
            p->lp = dataPeriodCode;
            p->pb = dataPeriodCode;
            p->fb = 273;
            break;
        default:
            fprintf(stderr, "Unsupported %d PgRC coding level for LZMA compression.\n", coder_level);
            exit(EXIT_FAILURE);
    }
    p->numThreads = numThreads;
    p->reduceSize = dataLength;
}

/*
LzmaCompress
------------

outPropsSize -
In:  the pointer to the size of outProps buffer; *outPropsSize = LZMA_PROPS_SIZE = 5.
Out: the pointer to the size of written properties in outProps buffer; *outPropsSize = LZMA_PROPS_SIZE = 5.

Out:
destLen  - processed output size
        Returns:
SZ_OK               - OK
SZ_ERROR_MEM        - Memory allocation error
SZ_ERROR_PARAM      - Incorrect paramater
SZ_ERROR_OUTPUT_EOF - output buffer overflow
SZ_ERROR_THREAD     - errors in multithreading functions (only for Mt version)
*/

MY_STDAPI LzmaCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                       uint8_t coder_level, int numThreads, int coder_param) {
    CLzmaEncProps props;
    LzmaEncProps_Init(&props);
    LzmaEncProps_Set(&props, coder_level, srcLen, numThreads, coder_param);
    size_t propsSize = LZMA_PROPS_SIZE;
    destLen = srcLen + srcLen / 3 + 128;
    dest = new unsigned char[propsSize + destLen];
    int res = LzmaEncode(dest + LZMA_PROPS_SIZE, &destLen, src, srcLen, &props, dest, &propsSize, 0,
                      NULL, &g_Alloc, &g_Alloc);
    assert(propsSize == LZMA_PROPS_SIZE);
    destLen += propsSize;
    return res;
}

/*
LzmaUncompress
--------------
In:
  dest     - output data
  destLen  - output data size
  src      - input data
  srcLen   - input data size
Out:
  destLen  - processed output size
  srcLen   - processed input size
Returns:
  SZ_OK                - OK
  SZ_ERROR_DATA        - Data error
  SZ_ERROR_MEM         - Memory allocation arror
  SZ_ERROR_UNSUPPORTED - Unsupported properties
  SZ_ERROR_INPUT_EOF   - it needs more bytes in input buffer (src)
*/

MY_STDAPI LzmaUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen,
                         const unsigned char *props, size_t propsSize) {
    ELzmaStatus status;
    return LzmaDecode(dest, destLen, src, srcLen, props, (unsigned) propsSize, LZMA_FINISH_ANY, &status, &g_Alloc);
}

char* Compress(size_t &destLen, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
               int coder_param) {
    clock_t start_t = clock();
    unsigned char* dest = 0;
    int res = 0;
    switch (coder_type) {
        case LZMA_CODER:
            res = LzmaCompress(dest, destLen, (const unsigned char*) src, srcLen, coder_level, 1, coder_param);
            break;
        case LZMA2_CODER:
        case PPMD7_CODER:
        default:
            fprintf(stderr, "Unsupported coder type: %d.\n", coder_type);
            exit(EXIT_FAILURE);
    }

    if (res != SZ_OK) {
        fprintf(stderr, "Error during compression (code: %d).\n", res);
        exit(EXIT_FAILURE);
    }
    cout << "compressed " << srcLen << " bytes to " << destLen << " bytes (ratio "
        << PgSAHelpers::toString(((double) destLen)/srcLen, 3) << ") in "
        << PgSAHelpers::clock_millis(start_t) << " msec." << endl;

    return (char*) dest;
}

void Uncompress(char* dest, size_t destLen, const char *src, size_t srcLen, uint8_t coder_type) {
    int res = 0;

    size_t outLen = destLen;
    size_t propsSize = LZMA_PROPS_SIZE;
    size_t srcBufSize = srcLen - LZMA_PROPS_SIZE;

    switch (coder_type) {
    case LZMA_CODER:
        res = LzmaUncompress((unsigned char*) dest, &outLen,
                             (unsigned char*) src + LZMA_PROPS_SIZE, &srcBufSize, (unsigned char*) src, propsSize);
    break;
    case LZMA2_CODER:
    case PPMD7_CODER:
    default:
    fprintf(stderr, "Unsupported coder type: %d.\n", coder_type);
    exit(EXIT_FAILURE);
    }
    assert(outLen == destLen);
    assert(srcBufSize == srcLen - LZMA_PROPS_SIZE);

    if (res != SZ_OK) {
        fprintf(stderr, "Error during decompression (code: %d).\n", res);
        exit(EXIT_FAILURE);
    }
}

void writeCompressed(ostream &dest, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
                     int coder_param) {
    PgSAHelpers::writeValue<uint64_t>(dest, srcLen, false);
    if (srcLen == 0) {
        cout << "skipped compression (0 bytes)." << endl;
        return;
    }
    size_t compLen = 0;
    char* compSeq = Compress(compLen, src, srcLen, coder_type, coder_level, coder_param);
    PgSAHelpers::writeValue<uint64_t>(dest, compLen, false);
    PgSAHelpers::writeValue<uint8_t>(dest, coder_type, false);
    PgSAHelpers::writeArray(dest, (void*) compSeq, compLen);
    delete(compSeq);
}

void writeCompressed(ostream &dest, const string srcStr, uint8_t coder_type, uint8_t coder_level, int coder_param) {
    writeCompressed(dest, srcStr.data(), srcStr.length(), coder_type, coder_level, coder_param);
}
void readCompressed(istream &src, string& dest) {
    size_t destLen = 0;
    size_t srcLen = 0;
    uint8_t coder_type = 0;
    PgSAHelpers::readValue<uint64_t>(src, destLen, false);
    dest.resize(destLen);
    if (destLen == 0)
        return;
    PgSAHelpers::readValue<uint64_t>(src, srcLen, false);
    PgSAHelpers::readValue<uint8_t>(src, coder_type, false);
    const char* srcArray = (const char*) PgSAHelpers::readArray(src, srcLen);
    Uncompress((char*) dest.data(), destLen, srcArray, srcLen, coder_type);
    delete(srcArray);
}