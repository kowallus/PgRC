/* FSECoder library wrapper */

#include "FSECoder.h"

#include "fse/fse.h"
#include "fse/huf.h"
#include "lzma/CpuArch.h"

#include "CodersLib.h"

#include <assert.h>

const static uint8_t FSE_CODER_PROPS_SIZE = 4;

MY_STDAPI FSECompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                        FSECoderProps* props, double estimated_compression) {

    if (props->getNumThreads() != 1) {
        cout << "Unsupported number of threads " << props->getNumThreads() << " in FSE coder." << endl;
        exit(EXIT_FAILURE);
    }
    size_t propsSize = FSE_CODER_PROPS_SIZE;
    size_t maxDestSize = propsSize + (srcLen + srcLen / 3) * estimated_compression + 1024;
    size_t fseBounds = props->mode == FSECoderProps::FSE_MODE ? FSE_compressBound(srcLen) : HUF_compressBound(srcLen);
    if (maxDestSize < fseBounds)
        maxDestSize = fseBounds;
    try {
        dest = new unsigned char[maxDestSize];
    } catch (const std::bad_alloc& e) {
        cout << "Allocation failed: " << e.what() << endl;
        maxDestSize -= srcLen / 3 * estimated_compression;
        dest = new unsigned char[maxDestSize];
    }
    *dest = props->mode;
    *(dest + 1) = props->blockSizeOrder;
    *(dest + 2) = props->maxSymbolValue;
    *(dest + 3) = props->tableLog;
    destLen = propsSize;

    size_t chunkSize = ((size_t) 1) << props->blockSizeOrder;
    size_t remainingSize = srcLen;
    unsigned char* destPtr = dest + propsSize;
    const unsigned char* srcPtr = src;
    while (remainingSize > 0) {
        size_t chunkLen = remainingSize > chunkSize ? chunkSize : remainingSize;
        uint8_t chunkSizeBytes = remainingSize == chunkLen ? 0 : (chunkLen > UINT32_MAX ? 8 : 4);
        destLen += chunkSizeBytes;
        size_t fse_res;
        if (props->mode == FSECoderProps::FSE_MODE)
            fse_res = FSE_compress2(destPtr + chunkSizeBytes, maxDestSize - destLen,
                                    srcPtr, chunkLen, props->maxSymbolValue, props->tableLog);
        else
            fse_res = HUF_compress2(destPtr + chunkSizeBytes, maxDestSize - destLen,
                                    srcPtr, chunkLen, props->maxSymbolValue, props->tableLog);
        if (FSE_isError(fse_res)) {
            cout << "Warning: fse_coder error " << FSE_getErrorName(fse_res) << "" << endl;
            return SZ_ERROR_FAIL;
        }
        if (fse_res == 0) {
            if (srcLen == 1) {
                fse_res = 1;
            } else {
                destLen = SIZE_MAX;
                return SZ_OK;
            };
        }
        destLen += fse_res;
        if (fse_res == 1) {
            destPtr[chunkSizeBytes] = src[0];
        }
        if (chunkSizeBytes == 8)
            SetUi64(destPtr, fse_res)
        else if (chunkSizeBytes == 4)
            SetUi32(destPtr, (uint32_t) fse_res);
        srcPtr += chunkLen;
        destPtr += chunkSizeBytes + fse_res;
        remainingSize -= chunkLen;
    }
    return SZ_OK;
}

MY_STDAPI FSEUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen, ostream* logout) {
    size_t propsSize = FSE_CODER_PROPS_SIZE;
    unsigned char propsBuf[FSE_CODER_PROPS_SIZE];
    memcpy(propsBuf, src, propsSize);
    uint8_t mode = propsBuf[0];
    if (mode != FSECoderProps::FSE_MODE && mode != FSECoderProps::HUF_MODE) {
        cout << "Unsupported mode " << mode << " in FSE coder." << endl;
        exit(EXIT_FAILURE);
    }
    uint8_t blockSizeOrder = propsBuf[1];
    uint8_t maxSymbolValue = propsBuf[2];
    uint8_t tableLog = propsBuf[3];
    srcLen -= propsSize;

    *logout << string("... fse_coder (") + (mode == FSECoderProps::HUF_MODE ? "huffman mode; ": "")
        + "blockSizeOrder = " + to_string(blockSizeOrder) + "; max_symbol = "
        + to_string(maxSymbolValue) + "; tableLog = " + to_string(tableLog) + ") ... ";

    size_t chunkSize = ((size_t) 1) << blockSizeOrder;
    size_t remainingSize = *destLen;
    unsigned char* destPtr = dest;
    const unsigned char* srcPtr = src + propsSize;
    while (remainingSize > 0) {
        size_t chunkLen = remainingSize > chunkSize ? chunkSize : remainingSize;
        uint8_t chunkSizeBytes = remainingSize == chunkLen ? 0 : (chunkLen > UINT32_MAX ? 8 : 4);
        size_t srcChunkLen = srcLen;
        if (chunkSizeBytes == 8)
            srcChunkLen = GetUi64(srcPtr);
        else if (chunkSizeBytes == 4)
            srcChunkLen = GetUi32(srcPtr);
        srcPtr += chunkSizeBytes;
        srcLen -= chunkSizeBytes;

        size_t fse_res;
        if (srcChunkLen == 1) {
            memset((void *) destPtr, srcPtr[0], chunkLen);
            fse_res = *destLen;
        } else if (mode == FSECoderProps::FSE_MODE)
            fse_res = FSE_decompress(destPtr, chunkLen, srcPtr, srcChunkLen);
        else
            fse_res = HUF_decompress(destPtr, chunkLen, srcPtr, srcChunkLen);
        if (FSE_isError(fse_res)) {
            cout << "Warning: fse_coder error " << FSE_getErrorName(fse_res) << "" << endl;
            return SZ_ERROR_FAIL;
        }
        assert(chunkLen == fse_res);
        srcPtr += srcChunkLen;
        srcLen -= srcChunkLen;
        destPtr += chunkLen;
        remainingSize -= chunkLen;
    }

    return SZ_OK;
}
