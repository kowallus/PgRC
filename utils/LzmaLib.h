/* Based on LzmaLib.h -- LZMA library interface
2013-01-18 : Igor Pavlov : Public domain */

#ifndef __LZMA_LIB_H
#define __LZMA_LIB_H

#include "../lzma/7zTypes.h"

EXTERN_C_BEGIN

#define MY_STDAPI int MY_STD_CALL

#define LZMA_PROPS_SIZE 5

const static int LZMA_CODER = 1;
const static int LZMA2_CODER = 2;
const static int PPMD7_CODER = 3;

const static int PGRC_DATAPERIODCODE_8_t = 0;
const static int PGRC_DATAPERIODCODE_16_t = 1;
const static int PGRC_DATAPERIODCODE_32_t = 2;
const static int PGRC_DATAPERIODCODE_64_t = 3;

const static int PGRC_CODER_LEVEL_NORMAL = 3;
const static int PGRC_CODER_LEVEL_MAXIMUM = 4;

char* Compress(size_t &destLen, const char *src, size_t srcLen, int coder_type, int coder_level);

void Uncompress(char* dest, size_t destLen, const char *src, size_t srcLen, int coder_type);

EXTERN_C_END

#endif