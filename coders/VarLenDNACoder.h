#ifndef PGTOOLS_VARLENDNACODER_H
#define PGTOOLS_VARLENDNACODER_H

#include "../utils/helper.h"
#include "CodersLib.h"

namespace PgHelpers {

    const static int STATIC_CODES_CODER_PARAM = 0;
    const static int DYNAMIC_CODES_CODER_PARAM = 1;

    class VarLenDNACoderProps: public CoderProps {
    public:
        const uint8_t mode;
        const uint8_t staticCodeBookID;

        VarLenDNACoderProps(uint8_t staticCodeBookID = 0)
                : CoderProps(VARLEN_DNA_CODER), mode(STATIC_CODES_CODER_PARAM), staticCodeBookID(staticCodeBookID) { }

        string log() {
            return " var-len (mode = " + to_string(mode) + "; cbId = " + to_string(staticCodeBookID) + ")";
        }

        virtual ~VarLenDNACoderProps() { };
    };

    class VarLenDNACoder {
    private:

        const static uint32_t MAX_NUMBER_OF_CODES = UINT8_MAX + 1;
        const static uint8_t MAX_CODE_LENGTH = 4;
        const static uint32_t CODE_LUT_SIZE = 1 << 27;
        const static uint32_t CODE_LUT_MASK = CODE_LUT_SIZE - 1;
        const static uint8_t NOT_FOUND_CODE = 0;

        const static size_t VAR_LEN_PROPS_SIZE = 2;
        const static size_t MAX_CODEBOOK_SIZE = 2 + (MAX_CODE_LENGTH + 1) * MAX_NUMBER_OF_CODES;

        uint8_t* codeLUT;
        char codeBook[MAX_NUMBER_OF_CODES][MAX_CODE_LENGTH + 1] = {};

        const static string AG_EXTENDED_CODES;
        const static string SYNC_ON_A_CODES;
        const static string AG_SHORT_EXTENDED_CODES;

        void initUsing(const string &codes);

    public:
        VarLenDNACoder(const string &codes);
        VarLenDNACoder(uint8_t staticCodeBookID);
        ~VarLenDNACoder();

        static int Compress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                            VarLenDNACoderProps* props);
        static int Uncompress(unsigned char *dest, size_t *destLen, unsigned char *src, size_t srcLen);

        void writeBook(unsigned char *dest, size_t &destLen);
        void encode(unsigned char *dest, size_t &destLen, const unsigned char *src, size_t srcLen);
        static string readBook(unsigned char *src);
        int decode(unsigned char *dest, size_t expDestLen, unsigned char *src, size_t srcLen);

        enum CODEBOOK_ID {
            AG_EXTENDED_CODES_ID,
            SYNC_ON_A_CODES_ID,
            AG_SHORT_EXTENDED_CODES_ID,
            VARLEN_CODEBOOKS_COUNT
        };

        constexpr static double COMPRESSION_ESTIMATION = 0.4;
    };

}

#endif //PGTOOLS_VARLENDNACODER_H

