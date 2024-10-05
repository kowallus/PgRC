#ifndef PGTOOLS_LZMACODER_H
#define PGTOOLS_LZMACODER_H

#include "lzma/7zTypes.h"
#include "../utils/helper.h"

#define MY_STDAPI int

#define UNCOMPRESS_BUFFER_SIZE 8192

#define LZMA_PROPS_SIZE 5

#include "CodersLib.h"
#include "lzma/LzmaEnc.h"

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

class LzmaCoderProps: public CoderProps {
private:
    CLzmaEncProps *p;

public:
    LzmaCoderProps(uint8_t level = 5, uint32_t dictSize = 0, int lc_highBits = -1, int lp_dataPeriodCode = -1,
                   int pb_dataPeriodCode = -1, int fb_wordSize = -1, int algo = -1, int numThreads = -1)
            : CoderProps(LZMA_CODER) {
        p = new CLzmaEncProps();
        LzmaEncProps_Init(p);
        p->level = level;
        p->dictSize = dictSize;
        p->lc = lc_highBits;
        p->lp = lp_dataPeriodCode;
        p->pb = pb_dataPeriodCode;
        p->fb = fb_wordSize;
        p->algo = algo;
        p->numThreads = numThreads;
    }

    int getNumThreads() const override { return p->numThreads; };

    string log() override {
        return " lzma (algo = " + to_string(p->algo) + "; dictSize = " + to_string(p->dictSize >> 20)
                             + "MB; lc = " + to_string(p->lc) + "; lp = " + to_string(p->lp) + "; pb = "
                             + to_string(p->pb) + "; fb = " + to_string(p->fb)
                             + "; th = " + to_string(p->numThreads) + ")";
    }

    ~LzmaCoderProps() override { delete p; };

    CLzmaEncProps* getProps() const { return p; }
};

MY_STDAPI LzmaCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                       CLzmaEncProps* props, double estimated_compression);

MY_STDAPI LzmaUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen,
        ostream* logout = PgHelpers::devout);

#endif //PGTOOLS_LZMACODER_H
