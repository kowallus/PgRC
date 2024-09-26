/* LZMA library wrapper
2015-06-13 : Igor Pavlov : Public domain */

#include "PpmdCoder.h"

#include "LzmaCoder.h"
#include "lzma/Ppmd7.h"
#include "lzma/Alloc.h"

#include "CodersLib.h"

void Ppmd7_SetMemSize(uint32_t &memSize, size_t dataLength) {
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
    CByteOutBufWrap *p = Z7_CONTAINER_FROM_VTBL_CLS(pp, CByteOutBufWrap, vt);
    // Z7_CONTAINER_FROM_VTBL_TO_DECL_VAR_pp_vt_p(CByteOutBufWrap)

    Byte *dest = p->Cur;
    if (dest == p->Lim) {
        p->Res = SZ_ERROR_OUTPUT_EOF;
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
                        PpmdCoderProps* props, double estimated_compression) {

    if (props->getNumThreads() != 1) {
        cerr << "Unsupported number of threads " << props->getNumThreads() << " in ppmd compressor." << endl;
        exit(EXIT_FAILURE);
    }
    CPpmd7 ppmd;
    Ppmd7_SetMemSize(props->memSize, srcLen);
    Ppmd7_Construct(&ppmd);
    if (!Ppmd7_Alloc(&ppmd, props->memSize, &g_Alloc))
        return SZ_ERROR_MEM;
    Ppmd7z_Init_RangeEnc(&ppmd);
    Ppmd7_Init(&ppmd, props->order);
    size_t propsSize = 5;
    size_t maxDestSize = propsSize + (srcLen + srcLen / 3) * estimated_compression + 128;
    try {
        dest = new unsigned char[maxDestSize];
    } catch (const std::bad_alloc& e) {
        cerr << "Allocation failed: " << e.what() << endl;
        maxDestSize -= srcLen / 3 * estimated_compression;
        dest = new unsigned char[maxDestSize];
    }
    *dest = (unsigned char) props->order;
    SetUi32(dest + 1, props->memSize);
    CByteOutBufWrap _outStream(dest + propsSize, maxDestSize - propsSize);
    ppmd.rc.enc.Stream = &_outStream.vt;

    Ppmd7z_EncodeSymbols(&ppmd, src, src + srcLen);

    Ppmd7z_Flush_RangeEnc(&ppmd);
    destLen = propsSize + _outStream.GetProcessed();
    Ppmd7_Free(&ppmd, &g_Alloc);
    if (_outStream.Res == SZ_ERROR_OUTPUT_EOF)
        return SZ_ERROR_OUTPUT_EOF;
    return SZ_OK;
}

struct CByteInBufWrap
{
    IByteIn vt;
    const Byte *Cur = 0;
    const Byte *Lim = 0;
    const Byte *Buf;
    UInt64 Size;
    bool Extra;
    int Res;

    istream *src;
    unsigned char* const srcBuf = 0;
    Int64 srcLeftCount;

    bool reloadSrcBuf();

    CByteInBufWrap(const unsigned char *buf, size_t bufLen);
    CByteInBufWrap(istream &src, size_t bufLen);
    ~CByteInBufWrap() {
        delete[] srcBuf;
    }

    UInt64 GetProcessed() const { return srcBuf?(Size - srcLeftCount):(Cur - Buf); }
};

static Byte Wrap_ReadByte(const IByteIn *pp) throw()
{
    CByteInBufWrap *p = Z7_CONTAINER_FROM_VTBL_CLS(pp, CByteInBufWrap, vt);
    if (p->Cur != p->Lim)
        return *p->Cur++;
    p->Res = SZ_ERROR_INPUT_EOF;
    return 0;
}

CByteInBufWrap::CByteInBufWrap(const unsigned char *buf, size_t bufLen): Buf(buf), Size(bufLen)
{
    Cur = Buf;
    Lim = Buf + Size;
    Extra = false;
    Res = SZ_OK;
    vt.Read = Wrap_ReadByte;
}

static Byte Wrap_Stream_ReadByte(const IByteIn *pp) throw()
{
    CByteInBufWrap *p = Z7_CONTAINER_FROM_VTBL_CLS(pp, CByteInBufWrap, vt);
    if (p->Cur != p->Lim)
        return *p->Cur++;
    if (p->reloadSrcBuf())
        return *p->Cur++;
    p->Res = SZ_ERROR_INPUT_EOF;
    return 0;
}

bool CByteInBufWrap::reloadSrcBuf() {
    if (!srcLeftCount)
        return false;
    size_t srcBufSize = UNCOMPRESS_BUFFER_SIZE > srcLeftCount?srcLeftCount:UNCOMPRESS_BUFFER_SIZE;
    PgHelpers::readArray(*src, srcBuf, srcBufSize);
    Cur = srcBuf;
    Lim = srcBuf + srcBufSize;
    srcLeftCount -= srcBufSize;
    return true;
}

CByteInBufWrap::CByteInBufWrap(istream &src, size_t bufLen): src(&src),
                                                             srcBuf(new unsigned char[UNCOMPRESS_BUFFER_SIZE]), Size(bufLen)
{
    srcLeftCount = Size;
    reloadSrcBuf();
    Extra = false;
    Res = SZ_OK;
    vt.Read = Wrap_Stream_ReadByte;
}

MY_STDAPI PpmdUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen, ostream* logout) {
    size_t propsSize = LZMA_PROPS_SIZE;
    CPpmd7 ppmd;
    //LzmaDecProps_Set(&props, coder_level, srcLen, numThreads, coder_param); // extract method, support numThreads?
    Ppmd7_Construct(&ppmd);
    unsigned char propsBuf[LZMA_PROPS_SIZE];
    memcpy(propsBuf, src, propsSize);
    uint32_t memSize = GetUi32(propsBuf + 1);
    if (!Ppmd7_Alloc(&ppmd, memSize, &g_Alloc))
        return SZ_ERROR_MEM;
    unsigned int order = propsBuf[0];
    Ppmd7_Init(&ppmd, order);
    *logout << "... ppmd (mem = " << (memSize >> 20) << "MB; ord = " << order << ") ... ";
    CByteInBufWrap _inStream(src + propsSize, srcLen - propsSize);
    ppmd.rc.dec.Stream = &_inStream.vt;

    int Res = SZ_OK;

    if (!Ppmd7z_RangeDec_Init(&ppmd.rc.dec))
        Res = SZ_ERROR_DATA;
    else if (_inStream.Extra)
        Res = (_inStream.Res != SZ_OK ? _inStream.Res : SZ_ERROR_DATA);
    else {
        size_t i;
        for (i = 0; i < *destLen; i++) {
            int sym = Ppmd7z_DecodeSymbol(&ppmd);
            if (_inStream.Extra|| sym < 0)
                break;
            dest[i] = sym;
        }
        if (i != *destLen)
            Res = (_inStream.Res != SZ_OK ? _inStream.Res : SZ_ERROR_DATA);
        else if (_inStream.GetProcessed() != srcLen - propsSize || !Ppmd7z_RangeDec_IsFinishedOK(&ppmd.rc.dec))
            Res = SZ_ERROR_DATA;
    }

    Ppmd7_Free(&ppmd, &g_Alloc);
    return Res;
}
