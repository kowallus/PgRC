/* LZMA library wrapper
2015-06-13 : Igor Pavlov : Public domain */

#include "LzmaCoder.h"

#include <assert.h>
#include "lzma/LzmaDec.h"
#include "lzma/Alloc.h"

#include "CodersLib.h"

MY_STDAPI LzmaCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                       CLzmaEncProps* props, double estimated_compression) {
    props->reduceSize = srcLen;
    if (srcLen >= UINT32_MAX)
        props->numThreads = 1;
    size_t propsSize = LZMA_PROPS_SIZE;
    size_t maxDestSize = propsSize + (srcLen + srcLen / 3) * estimated_compression + 128;
    try {
        dest = new unsigned char[maxDestSize];
    } catch (const std::bad_alloc& e) {
        cerr << "WARNING: Allocation failed: " << e.what() << endl;
        maxDestSize -= srcLen / 3 * estimated_compression;
        dest = new unsigned char[maxDestSize];
    }
    destLen = maxDestSize - propsSize;
    int res = LzmaEncode(dest + LZMA_PROPS_SIZE, &destLen, src, srcLen, props, dest, &propsSize, 0,
                         NULL, &g_Alloc, &g_Alloc);
    LzmaEncProps_Normalize(props);
    assert(propsSize == LZMA_PROPS_SIZE);
    destLen += propsSize;
    return res;
}

#define RC_INIT_SIZE 5

MY_STDAPI LzmaUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen, ostream* logout) {
    size_t propsSize = LZMA_PROPS_SIZE;
    *logout << "... lzma ... ";
    unsigned char propsBuf[LZMA_PROPS_SIZE];
    memcpy(propsBuf, src, propsSize);

    CLzmaDec p;
    SRes res;
    SizeT outSize = *destLen, inSize = srcLen - propsSize;
    *destLen = srcLen = 0;
    ELzmaStatus status = LZMA_STATUS_NOT_SPECIFIED;
    if (inSize < RC_INIT_SIZE)
        return SZ_ERROR_INPUT_EOF;
    LzmaDec_Construct(&p);
    RINOK(LzmaDec_AllocateProbs(&p, propsBuf, propsSize, &g_Alloc));
    p.dic = dest;
    p.dicBufSize = outSize;
    LzmaDec_Init(&p);
    Int64 srcLeftCount = inSize;
    do {
        size_t srcBufSize = srcLeftCount;
        srcLen = srcBufSize;
        srcLeftCount -= srcBufSize;
        res = LzmaDec_DecodeToDic(&p, outSize, src + propsSize, &srcLen, LZMA_FINISH_ANY, &status);
    } while (srcLeftCount && status == LZMA_STATUS_NEEDS_MORE_INPUT);
    *destLen = p.dicPos;
    if (res == SZ_OK && status == LZMA_STATUS_NEEDS_MORE_INPUT)
        res = SZ_ERROR_INPUT_EOF;
    LzmaDec_FreeProbs(&p, &g_Alloc);
    return res;
}
