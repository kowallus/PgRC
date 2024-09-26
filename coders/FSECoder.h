#ifndef PGTOOLS_FSECODER_H
#define PGTOOLS_FSECODER_H

#include "LzmaCoder.h"

class FSECoderProps: public CoderProps {
private:
    static const int FSE_BLOCKSIZE_ORDER_DEFAULT = 40;
    static const int HUF_BLOCKSIZE_ORDER_MAX = 17;
    static const int FSE_MAX_SYMBOL_VALUE = 255;
    static const int FSE_DEFAULT_TABLELOG = 11;

public:
    uint8_t mode;
    uint8_t blockSizeOrder;
    uint8_t maxSymbolValue;
    uint8_t tableLog;

    const static uint8_t FSE_MODE = 1;
    const static uint8_t HUF_MODE = 0;

    FSECoderProps(bool huffman, uint8_t maxSymbolValue = 0, uint8_t tableLog = 0, uint8_t blockSizeOrder = 0)
            : CoderProps(FSE_HUF_CODER, 1), mode(huffman ? HUF_MODE : FSE_MODE),
            maxSymbolValue(maxSymbolValue ? maxSymbolValue : FSE_MAX_SYMBOL_VALUE),
            tableLog(tableLog ? tableLog : FSE_DEFAULT_TABLELOG),
            blockSizeOrder(blockSizeOrder ? blockSizeOrder : (huffman ? HUF_BLOCKSIZE_ORDER_MAX : FSE_BLOCKSIZE_ORDER_DEFAULT)) { }

    string log() {
        return string(" fse_coder (") + (mode == HUF_MODE ? "huffman mode; ": "") + "blockSizeOrder = " + to_string(blockSizeOrder) + "; max_symbol = " + to_string(maxSymbolValue) + "; tableLog = " + to_string(tableLog) + ")";
    }

    virtual ~FSECoderProps() { };
};

MY_STDAPI FSECompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                        FSECoderProps* coder_props, double estimated_compression);

MY_STDAPI FSEUncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t srcLen,
                          ostream* logout = PgHelpers::devout);

#endif //PGTOOLS_FSECODER_H
