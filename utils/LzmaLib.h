/* Based on LzmaLib.h -- LZMA library interface
2013-01-18 : Igor Pavlov : Public domain */

#ifndef __LZMA_LIB_H
#define __LZMA_LIB_H

#include "../lzma/7zTypes.h"
#include "helper.h"
#include <vector>

using namespace std;

#define MY_STDAPI int MY_STD_CALL

#define LZMA_PROPS_SIZE 5

const static uint8_t LZMA_CODER = 1;
const static uint8_t LZMA2_CODER = 2;
const static uint8_t PPMD7_CODER = 3;

const static int PGRC_DATAPERIODCODE_8_t = 0;
const static int PGRC_DATAPERIODCODE_16_t = 1;
const static int PGRC_DATAPERIODCODE_32_t = 2;
const static int PGRC_DATAPERIODCODE_64_t = 3;
const static int PGRC_DATAPERIODCODE_128_t = 4;

const static uint8_t PGRC_CODER_LEVEL_FAST = 1;
const static uint8_t PGRC_CODER_LEVEL_NORMAL = 2;
const static uint8_t PGRC_CODER_LEVEL_MAX = 3;

const static double COMPRESSION_ESTIMATION_UINT8_BITMAP = 0.125;
const static double COMPRESSION_ESTIMATION_BASIC_DNA = 0.250;
const static double COMPRESSION_ESTIMATION_MIS_CNT = 0.5;
const static double COMPRESSION_ESTIMATION_MIS_SYM = 0.250;

double simpleUintCompressionEstimate(uint64_t dataMaxValue, uint64_t typeMaxValue);

char* Compress(size_t &destLen, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
        int coder_param = -1, double estimated_compression = 1);
void writeCompressed(ostream &dest, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
                     int coder_param = -1, double estimated_compression = 1);
void writeCompressed(ostream &dest, const string srcStr, uint8_t coder_type, uint8_t coder_level,
                     int coder_param = -1, double estimated_compression = 1);

void Uncompress(char* dest, size_t destLen, istream &src, size_t srcLen, uint8_t coder_type);
void readCompressed(istream &src, string& dest);

template<typename T>
void readCompressed(istream &src, vector<T>& dest) {
    size_t destLen = 0;
    size_t srcLen = 0;
    uint8_t coder_type = 0;
    PgSAHelpers::readValue<uint64_t>(src, destLen, false);
    if (destLen % sizeof(T)) {
        fprintf(stderr, "Invalid output size %zu for decompressing to the vector of %zu-byte elements",
                destLen, sizeof(T));
    }
    dest.resize(destLen / sizeof(T));
    if (destLen == 0)
        return;
    PgSAHelpers::readValue<uint64_t>(src, srcLen, false);
    PgSAHelpers::readValue<uint8_t>(src, coder_type, false);
    Uncompress((char*) dest.data(), destLen, src, srcLen, coder_type);
}

#endif