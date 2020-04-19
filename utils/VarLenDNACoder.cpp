#include "VarLenDNACoder.h"
#include "../lzma/7zTypes.h"
#include <cassert>
#include <cstring>

PgSAHelpers::VarLenDNACoder::VarLenDNACoder(const string &codes) {
    initUsing(codes);
}

void PgSAHelpers::VarLenDNACoder::initUsing(const string &codes) {
    int sPos = 0;
    int pos = 0;
    int counter = 0;
    while ((sPos = codes.find('\n', pos)) != string::npos && counter <= UINT8_MAX) {
        strncpy(codeBook[counter++], codes.data() + pos, sPos - pos);
        pos = sPos + 1;
    };
    if (counter > UINT8_MAX) {
        fprintf(stderr, "Too many codes in the var-len coder book!");
        exit(EXIT_FAILURE);
    }
    strncpy(codeBook[counter++], codes.data() + pos, codes.length() - pos);
    if (strlen(codeBook[NOT_FOUND_CODE]) != 1) {
        fprintf(stderr, "First code in the var-len coder book should be a single symbol!");
        exit(EXIT_FAILURE);
    }

    codeLUT = new uint8_t[CODE_LUT_SIZE];
    memset(codeLUT, NOT_FOUND_CODE, CODE_LUT_SIZE);
    for (int i = 0; i < counter; i++) {
        uint32_t temp = 0;
        memcpy(&temp, codeBook[i], MAX_CODE_LENGTH);
        codeLUT[temp & CODE_LUT_MASK] = i;
    }
}

PgSAHelpers::VarLenDNACoder::~VarLenDNACoder() {
    delete(codeLUT);
}

void
PgSAHelpers::VarLenDNACoder::writeBook(unsigned char *dest, size_t &destLen) {
    uint8_t* destPtr = dest;
    for (size_t i = 0; i < MAX_NUMBER_OF_CODES; i++) {
        char* ptr = codeBook[i];
        while(*ptr)
            *(destPtr++) = *(ptr++);
        *(destPtr++) = '\n';
    }
    *(destPtr - 1) = 0;
    destLen = destPtr - dest;
}

void
PgSAHelpers::VarLenDNACoder::encode(unsigned char *dest, size_t &destLen, const unsigned char *src, size_t srcLen) {
    assert(MAX_CODE_LENGTH == 4);
    destLen = 0;
    size_t pos = 0;
    while (pos <= srcLen - 4) {
        uint32_t temp = 0;
        memcpy(&temp, src + pos, 4);
        temp &= CODE_LUT_MASK;
        if (codeLUT[temp] != NOT_FOUND_CODE) {
            dest[destLen++] = codeLUT[temp];
            pos += 4;
            continue;
        }
        temp &= 0x00FFFFFF;
        if (codeLUT[temp] != NOT_FOUND_CODE) {
            dest[destLen++] = codeLUT[temp];
            pos += 3;
            continue;
        }
        temp &= 0x0000FFFF;
        if (codeLUT[temp] != NOT_FOUND_CODE) {
            dest[destLen++] = codeLUT[temp];
            pos += 2;
            continue;
        }
        temp &= 0x000000FF;
        dest[destLen++] = codeLUT[temp];
        pos++;
    }
    while (pos < srcLen) {
        uint32_t temp = 0;
        memcpy(&temp, src + pos, srcLen - pos);
        temp &= CODE_LUT_MASK;
        temp &= 0x00FFFFFF;
        if (srcLen - pos >= 3 && codeLUT[temp] != NOT_FOUND_CODE) {
            dest[destLen++] = codeLUT[temp];
            pos += 3;
            continue;
        }
        temp &= 0x0000FFFF;
        if (srcLen - pos >= 2 && codeLUT[temp] != NOT_FOUND_CODE) {
            dest[destLen++] = codeLUT[temp];
            pos += 2;
            continue;
        }
        temp &= 0x000000FF;
        dest[destLen++] = codeLUT[temp];
        pos++;
    }
}

int PgSAHelpers::VarLenDNACoder::decode(unsigned char *dest, size_t expDestLen, unsigned char *src, size_t srcLen) {
    uint8_t* destPtr = dest;
    for (size_t i = 0; i < srcLen; i++) {
        uint8_t code = src[i];
        char* ptr = codeBook[code];
        while(*ptr)
            *(destPtr++) = *(ptr++);
    }
    size_t destLen = destPtr - dest;
    if (expDestLen != destLen) {
        fprintf(stderr, "Unexpected decoded length: %d (expected %d).\n", destLen, expDestLen);
        exit(EXIT_FAILURE);
    }
    return 0;
}

int PgSAHelpers::VarLenDNACoder::getCoderParam(uint8_t coder_mode_param, uint8_t coder_value_param)  {
    return coder_mode_param * 256 + coder_value_param;
};

int
PgSAHelpers::VarLenDNACoder::Compress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                                      int coder_param) {
    size_t headerSize = VAR_LEN_PROPS_SIZE;
    size_t maxDestSize = VAR_LEN_PROPS_SIZE + MAX_CODEBOOK_SIZE + (srcLen + srcLen / 3) * COMPRESSION_ESTIMATION + 128;
    try {
        dest = new unsigned char[maxDestSize];
    } catch (const std::bad_alloc& e) {
        cout << "WARNING: Allocation failed: " << e.what() << endl;
        maxDestSize -= srcLen / 3 * COMPRESSION_ESTIMATION;
        dest = new unsigned char[maxDestSize];
    }
    VarLenDNACoder* coder;
    uint8_t coder_mode_param = coder_param / 256;
    uint8_t coder_value_param = coder_param % 256;
    dest[0] = coder_mode_param;
    dest[1] = coder_value_param;
    switch (dest[0]) {
        case STATIC_CODES_CODER_PARAM:
            coder = new VarLenDNACoder(coder_value_param);
            break;
        default:
            fprintf(stderr, "Unsupported %d PgRC var-len coder parameter.\n", coder_param);
            exit(EXIT_FAILURE);
    }
    coder->writeBook(dest + headerSize, destLen);
    headerSize += destLen;
    coder->encode(dest + headerSize, destLen, src, srcLen);
    delete(coder);
    destLen += headerSize;
    if (destLen > maxDestSize)
        fprintf(stderr, "WARNING: exceeded maximum destination size: %d .\n", maxDestSize);

    return SZ_OK;
}

string PgSAHelpers::VarLenDNACoder::readBook(unsigned char *src) {
    return string((char*) src);
}

int PgSAHelpers::VarLenDNACoder::Uncompress(unsigned char *dest, size_t *destLen, unsigned char *src, size_t *srcLen) {
    int coder_param = src[0];
    size_t headerSize = VAR_LEN_PROPS_SIZE;
    VarLenDNACoder* coder;
    switch (coder_param) {
        case STATIC_CODES_CODER_PARAM:
        case DYNAMIC_CODES_CODER_PARAM:
            {
                string codeBook = readBook(src + headerSize);
                coder = new VarLenDNACoder(codeBook);
                headerSize += codeBook.length() + 1;
            }
            break;
        default:
            fprintf(stderr, "Unsupported %d PgRC var-len coder parameter.\n", coder_param);
            exit(EXIT_FAILURE);
    }

    int res = coder->decode(dest, *destLen, src + headerSize, (*srcLen) - headerSize);
    delete(coder);
    return res;
}

PgSAHelpers::VarLenDNACoder::VarLenDNACoder(uint8_t staticCodeBookID) {
    switch(staticCodeBookID) {
        case AG_EXTENDED_CODES_ID:
            initUsing(AG_EXTENDED_CODES);
            break;
        case SYNC_ON_A_CODES_ID:
            initUsing(SYNC_ON_A_CODES);
            break;
        case AG_SHORT_EXTENDED_CODES_ID:
            initUsing(AG_SHORT_EXTENDED_CODES);
            break;
        default:
            fprintf(stderr, "Unknown var-len coder static codebook id: %d.\n", staticCodeBookID);
            exit(EXIT_FAILURE);
    }
}

const string PgSAHelpers::VarLenDNACoder::AG_EXTENDED_CODES =
        "A\nAAA\nAAAA\nGAAA\nCAA\nACAA\nGCAA\nGAA\nAGAA\nGGAA\nTAA\nATAA\nGTAA\nACA\nAACA\nGACA\nCCA\nACCA\nGCCA\nGCA\n"
        "AGCA\nGGCA\nTCA\nATCA\nGTCA\nAGA\nAAGA\nGAGA\nCGA\nACGA\nGCGA\nGGA\nAGGA\nGGGA\nTGA\nATGA\nGTGA\nATA\nAATA\n"
        "GATA\nCTA\nACTA\nGCTA\nGTA\nAGTA\nGGTA\nTTA\nATTA\nGTTA\nTG%\nT%\nAT%\nCT%\nGT%\nTT%\n"
        "N\nAN\nAAN\nCAN\nGAN\nNAN\nTAN\n"
        "\n\n"
        "C\nAAC\nAAAC\nGAAC\nCAC\nACAC\nGCAC\nGAC\nAGAC\nGGAC\nTAC\nATAC\nGTAC\nACC\nAACC\nGACC\nCCC\nACCC\nGCCC\nGCC\n"
        "AGCC\nGGCC\nTCC\nATCC\nGTCC\nAGC\nAAGC\nGAGC\nCGC\nACGC\nGCGC\nGGC\nAGGC\nGGGC\nTGC\nATGC\nGTGC\nATC\nAATC\n"
        "GATC\nCTC\nACTC\nGCTC\nGTC\nAGTC\nGGTC\nTTC\nATTC\nGTTC\n"
        "CN\nACN\nCCN\nGCN\nNCN\nTCN\nGN\nAGN\nCGN\nGGN\nNGN\nTGN\n"
        "\n\n\n"
        "G\nAAG\nAAAG\nGAAG\nCAG\nACAG\nGCAG\nGAG\nAGAG\nGGAG\nTAG\nATAG\nGTAG\nACG\nAACG\nGACG\nCCG\nACCG\nGCCG\nGCG\n"
        "AGCG\nGGCG\nTCG\nATCG\nGTCG\nAGG\nAAGG\nGAGG\nCGG\nACGG\nGCGG\nGGG\nAGGG\nGGGG\nTGG\nATGG\nGTGG\nATG\nAATG\n"
        "GATG\nCTG\nACTG\nGCTG\nGTG\nAGTG\nGGTG\nTTG\nATTG\nGTTG\n"
        "TN\nATN\nCTN\nGTN\nNTN\nTTN\nNN\n"
        "\n\n\n\n\n\n\n\n"
        "T\nAAT\nAAAT\nGAAT\nCAT\nACAT\nGCAT\nGAT\nAGAT\nGGAT\nTAT\nATAT\nGTAT\nACT\nAACT\nGACT\nCCT\nACCT\nGCCT\nGCT\n"
        "AGCT\nGGCT\nTCT\nATCT\nGTCT\nAGT\nAAGT\nGAGT\nCGT\nACGT\nGCGT\nGGT\nAGGT\nGGGT\nTGT\nATGT\nGTGT\nATT\nAATT\n"
        "GATT\nCTT\nACTT\nGCTT\nGTT\nAGTT\nGGTT\nTTT\nATTT\nGTTT\n"
        "%\nA%\nAA%\nCA%\nGA%\nTA%\nC%\nAC%\nCC%\nGC%\nTC%\nG%\nAG%\nCG%\nGG%";

const string PgSAHelpers::VarLenDNACoder::SYNC_ON_A_CODES =
        "A\nAA\nAAA\nAAAA\nACAA\nAGAA\nATAA\nACA\nAACA\nACCA\nAGCA\nATCA\nAGA\nAAGA\nACGA\nAGGA\nATGA\nATA\nAATA\nACTA\n"
        "AGTA\nATTA\n"
        "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
        "C\nAC\nAAC\nAAAC\nACAC\nAGAC\nATAC\nCC\nACC\nAACC\nCCC\nACCC\nGCC\nAGCC\nTCC\nATCC\nGC\nAGC\nAAGC\nCGC\nACGC\n"
        "GGC\nAGGC\nTGC\nATGC\nTC\nATC\nAATC\nCTC\nACTC\nGTC\nAGTC\nTTC\nATTC\n"
        "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
        "G\nAG\nAAG\nAAAG\nACAG\nAGAG\nATAG\nCG\nACG\nAACG\nCCG\nACCG\nGCG\nAGCG\nTCG\nATCG\nGG\nAGG\nAAGG\nCGG\nACGG\n"
        "GGG\nAGGG\nTGG\nATGG\nTG\nATG\nAATG\nCTG\nACTG\nGTG\nAGTG\nTTG\nATTG\n"
        "N\nAN\nAAN\nCAN\nGAN\nNAN\nTAN\n"
        "CN\nACN\nCCN\nGCN\nNCN\nTCN\nGN\nAGN\nCGN\nGGN\nNGN\nTGN\n"
        "TN\nATN\nCTN\nGTN\nNTN\nTTN\nNN\n"
        "\n\n\n\n"
        "T\nAT\nAAT\nAAAT\nACAT\nAGAT\nATAT\nCT\nACT\n"
        "AACT\nCCT\nACCT\nGCT\nAGCT\nTCT\nATCT\nGT\nAGT\nAAGT\nCGT\nACGT\nGGT\nAGGT\nTGT\nATGT\nTT\nATT\nAATT\nCTT\n"
        "ACTT\nGTT\nAGTT\nTTT\nATTT\n"
        "%\nA%\nAA%\nCA%\nGA%\nTA%\nC%\nAC%\nCC%\nGC%\nTC%\nG%\nAG%\nCG%\nGG%";

const string PgSAHelpers::VarLenDNACoder::AG_SHORT_EXTENDED_CODES =
        "A\nAA\n\n\nAAA\n\n\n\n\n\nGAA\n\n\nCA\n\n\nACA\n\n\n\n\n\nGCA\n\n\nGA\n\n\nAGA\n\n\n\n\n\nGGA\n\n\nTA\n\n\n"
        "ATA\n\n\n\n\n\nGTA\n\n\nTG%\nT%\nAT%\nCT%\nGT%\nTT%\n"
        "N\nAN\nAAN\nCAN\nGAN\nNAN\nTAN\n"
        "\n\n"
        "C\nAC\n\n\nAAC\n\n\n\n\n\nGAC\n\n\nCC\n\n\nACC\n\n\n\n\n\nGCC\n\n\nGC\n\n\nAGC\n\n\n\n\n\nGGC\n\n\nTC\n\n\n"
        "ATC\n\n\n\n\n\nGTC\n"
        "CN\nACN\nCCN\nGCN\nNCN\nTCN\nGN\nAGN\nCGN\nGGN\nNGN\nTGN\n"
        "\n\n\n\n\n"
        "G\nAG\n\n\nAAG\n\n\n\n\n\nGAG\n\n\nCG\n\n\nACG\n\n\n\n\n"
        "\nGCG\n\n\nGG\n\n\nAGG\n\n\n\n\n\nGGG\n\n\nTG\n\n\nATG\n\n\n\n\n\nGTG\n"
        "TN\nATN\nCTN\nGTN\nNTN\nTTN\nNN\n"
        "\n\n\n\n\n\n\n\n\n\n"
        "T\nAT\n\n\nAAT\n\n\n\n\n\nGAT\n\n\nCT\n\n\nACT\n\n\n\n\n\nGCT\n\n\nGT\n\n\nAGT\n\n\n\n\n\nGGT\n\n\nTT\n\n\n"
        "ATT\n\n\n\n\n\nGTT\n\n\n"
        "%\nA%\nAA%\nCA%\nGA%\nTA%\nC%\nAC%\nCC%\nGC%\nTC%\nG%\nAG%\nCG%\nGG%";