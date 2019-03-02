/* Based on LzmaLib.h -- LZMA library interface
2013-01-18 : Igor Pavlov : Public domain */

#ifndef __LZMA_LIB_H
#define __LZMA_LIB_H

#include "../lzma/7zTypes.h"
#include "helper.h"

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

const static uint8_t PGRC_CODER_LEVEL_NORMAL = 3;
const static uint8_t PGRC_CODER_LEVEL_MAXIMUM = 4;

char* Compress(size_t &destLen, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
        int coder_param = -1);
void writeCompressed(ostream &dest, const char *src, size_t srcLen, uint8_t coder_type, uint8_t coder_level,
                     int coder_param = -1);
void writeCompressed(ostream &dest, const string srcStr, uint8_t coder_type, uint8_t coder_level,
                     int coder_param = -1);

void Uncompress(char* dest, size_t destLen, const char *src, size_t srcLen, uint8_t coder_type);
void readCompressed(istream &src, string& dest);

#endif