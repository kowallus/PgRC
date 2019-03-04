/* LzmaLib.c -- LZMA library wrapper
2015-06-13 : Igor Pavlov : Public domain */

#include <assert.h>
#include "../lzma/Alloc.h"
#include "../lzma/LzmaDec.h"
#include "../lzma/LzmaEnc.h"
#include "../lzma/Ppmd7.h"
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
    cout << " lzma (level = " << p->level << "; dictSize = " << (p->dictSize >> 20) << "MB; lp,pb = " << p->lp << ") ...";
}

void Ppmd7_SetProps(uint32_t &memSize, uint8_t coder_level, size_t dataLength, int& order_param) {
    switch(coder_level) {
        case PGRC_CODER_LEVEL_NORMAL:
            memSize = (uint32_t) 128 << 20;
            if (order_param > 2)
                order_param--;
            break;
        case PGRC_CODER_LEVEL_MAXIMUM:
            memSize = (uint32_t) 192 << 20;
            break;
        default:
            fprintf(stderr, "Unsupported %d PgRC coding level for LZMA compression.\n", coder_level);
            exit(EXIT_FAILURE);
    }
    const unsigned kMult = 16;
    if (memSize / kMult > dataLength)
    {
        for (unsigned i = 16; i <= 31; i++)
        {
            UInt32 m = (UInt32)1 << i;
            if (dataLength <= m / kMult)
            {
                if (memSize > m)
                    memSize = m;
                break;
            }
        }
    }
}


/*
Compress
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

struct CByteOutBufWrap
{
    IByteOut vt;
    Byte *Cur;
    const Byte *Lim;
    Byte *Buf;
    size_t Size;
    int Res;

    UInt64 GetProcessed() const { return Cur - Buf; }
    CByteOutBufWrap(unsigned char *buf, size_t bufLen) throw();
    ~CByteOutBufWrap() { }
};

static void Wrap_WriteByte(const IByteOut *pp, Byte b) throw()
{
    CByteOutBufWrap *p = CONTAINER_FROM_VTBL_CLS(pp, CByteOutBufWrap, vt);
    Byte *dest = p->Cur;
    if (dest == p->Lim) {
        p->Res == SZ_ERROR_OUTPUT_EOF;
        return;
    }
    *dest = b;
    p->Cur = ++dest;
}

CByteOutBufWrap::CByteOutBufWrap(unsigned char *buf, size_t bufLen) throw(): Buf(buf), Size(bufLen)
{
    Cur = Buf;
    Lim = Buf + Size;
    Res = SZ_OK;
    vt.Write = Wrap_WriteByte;
}

MY_STDAPI Ppmd7Compress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
              uint8_t coder_level, int numThreads, int coder_param) {
    if (numThreads != 1) {
        cout << "Unsupported number of threads " << numThreads << " in ppmd compressor." << endl;
        exit(EXIT_FAILURE);
    }
    CPpmd7 ppmd;
    uint32_t memSize = 0;
    Ppmd7_SetProps(memSize, coder_level, srcLen, coder_param);
    Ppmd7_Construct(&ppmd);
    if (!Ppmd7_Alloc(&ppmd, memSize, &g_Alloc))
        return SZ_ERROR_MEM;
    Ppmd7_Init(&ppmd, coder_param);
    cout << " ppmd (mem = " << (memSize >> 20) << "MB; ord = " << coder_param << ") ...";
    CPpmd7z_RangeEnc rEnc;
    Ppmd7z_RangeEnc_Init(&rEnc);
    size_t propsSize = 5;
    size_t maxDestSize = propsSize + srcLen + srcLen / 3 + 128; //TODO: better estimation of max size
    dest = new unsigned char[maxDestSize];
    *dest = (unsigned char) coder_param;
    SetUi32(dest + 1, memSize);
    CByteOutBufWrap _outStream(dest + propsSize, maxDestSize - propsSize);
    rEnc.Stream = &_outStream.vt;

    for(size_t i = 0; i < srcLen; i++)
        Ppmd7_EncodeSymbol(&ppmd, &rEnc, src[i]);

    Ppmd7z_RangeEnc_FlushData(&rEnc);
    destLen = propsSize + _outStream.GetProcessed();
    return SZ_OK;
}

/*
Uncompress
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

MY_STDAPI LzmaUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen) {
    size_t propsSize = LZMA_PROPS_SIZE;
    size_t srcBufSize = *srcLen - LZMA_PROPS_SIZE;
    ELzmaStatus status;
    return LzmaDecode(dest, destLen, src + propsSize, &srcBufSize, src, (unsigned) propsSize, LZMA_FINISH_ANY, &status, &g_Alloc);
}

struct CByteInBufWrap
{
    IByteIn vt;
    const Byte *Cur;
    const Byte *Lim;
    Byte *Buf;
    UInt32 Size;
    bool Extra;
    int Res;

    CByteInBufWrap(unsigned char *buf, size_t bufLen);
    ~CByteInBufWrap() { }

    UInt64 GetProcessed() const { return (Cur - Buf); }
};

static Byte Wrap_ReadByte(const IByteIn *pp) throw()
{
    CByteInBufWrap *p = CONTAINER_FROM_VTBL_CLS(pp, CByteInBufWrap, vt);
    if (p->Cur != p->Lim)
        return *p->Cur++;
    p->Res = SZ_ERROR_INPUT_EOF;
    return 0;
}

CByteInBufWrap::CByteInBufWrap(unsigned char *buf, size_t bufLen): Buf(buf), Size(bufLen)
{
    Cur = Buf;
    Lim = Buf + Size;
    Extra = false;
    Res = SZ_OK;
    vt.Read = Wrap_ReadByte;
}


MY_STDAPI PpmdUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen) {
    size_t propsSize = LZMA_PROPS_SIZE;
    CPpmd7 ppmd;
    //LzmaDecProps_Set(&props, coder_level, srcLen, numThreads, coder_param); // extract method, support numThreads?
    Ppmd7_Construct(&ppmd);
    uint32_t memSize = GetUi32(src + 1);
    if (!Ppmd7_Alloc(&ppmd, memSize, &g_Alloc))
        return SZ_ERROR_MEM;
    unsigned int order = src[0];
    Ppmd7_Init(&ppmd, order);
    CPpmd7z_RangeDec rDec;
    Ppmd7z_RangeDec_CreateVTable(&rDec);
    CByteInBufWrap _inStream((unsigned char *) src + propsSize, *srcLen - propsSize);
    rDec.Stream = &_inStream.vt;

    int Res = SZ_OK;

    if (!Ppmd7z_RangeDec_Init(&rDec))
        Res = SZ_ERROR_DATA;
    else if (_inStream.Extra)
        Res = (_inStream.Res != SZ_OK ? _inStream.Res : SZ_ERROR_DATA);
    else {
        size_t i;
        for (i = 0; i < *destLen; i++) {
            int sym = Ppmd7_DecodeSymbol(&ppmd, &rDec.vt);
            if (_inStream.Extra|| sym < 0)
                break;
            dest[i] = sym;
        }
        if (i != *destLen)
            Res = (_inStream.Res != SZ_OK ? _inStream.Res : SZ_ERROR_DATA);
        else if (_inStream.GetProcessed() != *srcLen - propsSize || !Ppmd7z_RangeDec_IsFinishedOK(&rDec))
            Res = SZ_ERROR_DATA;
    }

    Ppmd7_Free(&ppmd, &g_Alloc);
    return Res;
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
        case PPMD7_CODER:
            res = Ppmd7Compress(dest, destLen, (const unsigned char*) src, srcLen, coder_level, 1, coder_param);
            break;
        case LZMA2_CODER:
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

    switch (coder_type) {
    case LZMA_CODER:
        res = LzmaUncompress((unsigned char*) dest, &outLen,
                             (unsigned char*) src, &srcLen);
    break;
    case PPMD7_CODER:
        res = PpmdUncompress((unsigned char*) dest, &outLen,
                             (unsigned char*) src, &srcLen);
    break;
    case LZMA2_CODER:
    default:
    fprintf(stderr, "Unsupported coder type: %d.\n", coder_type);
    exit(EXIT_FAILURE);
    }
    assert(outLen == destLen);

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