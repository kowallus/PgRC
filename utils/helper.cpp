#include "helper.h"

#include "byteswap.h"

std::ostream *PgHelpers::logout = &std::cout;
std::ostream *PgHelpers::appout = &std::cout;
std::ostream *PgHelpers::devout = &std::cout;

int PgHelpers::numberOfThreads = 8;

NullBuffer null_buffer;
std::ostream null_stream(&null_buffer);

// TIME

clock_t checkpoint;

clock_t PgHelpers::clock_checkpoint() {
    checkpoint = clock();
    return checkpoint;
//    cout << "Clock reset!\n";
}

unsigned long long int PgHelpers::clock_millis(clock_t checkpoint) {
    return (clock() - checkpoint) * (unsigned long long int) 1000 / CLOCKS_PER_SEC;
}

unsigned long long int PgHelpers::clock_millis() {
    return clock_millis(checkpoint);
}


chrono::steady_clock::time_point chronocheckpoint;

chrono::steady_clock::time_point PgHelpers::time_checkpoint() {
    chronocheckpoint = chrono::steady_clock::now();
    return chronocheckpoint;
}

unsigned long long int PgHelpers::time_millis(chrono::steady_clock::time_point checkpoint) {
    chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now() - checkpoint);
    return (double)time_span.count() / 1000000.0;
}

unsigned long long int PgHelpers::time_millis() {
    return time_millis(chronocheckpoint);
}




const size_t chunkSize = 10000000;

void* PgHelpers::readArray(std::istream& in, size_t arraySizeInBytes) {

    if (in) {
        char* destArray = new char[arraySizeInBytes];
        readArray(in, destArray, arraySizeInBytes);
        return (void*) destArray;
    } else
        throw(errno);
}

void PgHelpers::readArray(std::istream& in, void* destArray, size_t arraySizeInBytes) {

    if (in) {
        size_t length = arraySizeInBytes;
        char* ptr = (char*) destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);
    } else
        throw(errno);
}

void PgHelpers::writeArray(std::ostream& out, void* srcArray, size_t arraySize, bool verbose) {
    size_t bytesLeft = arraySize;
    if (out) {
        while (bytesLeft > chunkSize) {
            out.write((char*) srcArray, chunkSize);
            srcArray = (void*) (((char*) srcArray) + chunkSize);
            bytesLeft -= chunkSize;
        }
        out.write((char*) srcArray, bytesLeft);

    } else
        throw(errno);
    if (verbose)
        cout << "Written " << arraySize << " bytes\n";
}

void* PgHelpers::readWholeArray(std::istream& in, size_t& arraySizeInBytes) {

    if (in) {
        in.seekg(0, in.end);
        size_t length = in.tellg();
        char* destArray = new char[length];
        in.seekg(0, in.beg);
        char* ptr = destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);

        arraySizeInBytes = length;
        return (void*) destArray;
    } else
        throw(errno);
}

void* PgHelpers::readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes) {
    time_checkpoint();

    std::ifstream in(srcFile.c_str(), std::ifstream::binary);

    void* destArray = PgHelpers::readWholeArray(in, arraySizeInBytes);

    cout << "Read " << arraySizeInBytes << " bytes from " << srcFile << " in " << time_millis() << " msec \n";

    return destArray;
}

void PgHelpers::writeArrayToFile(string destFile, void* srcArray, size_t arraySize) {
    time_checkpoint();

    std::ofstream out(destFile.c_str(), std::ios::out | std::ios::binary);

    PgHelpers::writeArray(out, srcArray, arraySize);

    cout << "Write " << arraySize << " bytes to " << destFile << " in " << time_millis() << " msec \n";
}

void PgHelpers::writeStringToFile(string destFile, const string &src) {
    writeArrayToFile(destFile, (void*) src.data(), src.length());
}

bool PgHelpers::plainTextWriteMode = false;

void PgHelpers::writeReadMode(std::ostream &dest, bool plainTextWriteMode) {
    dest << (plainTextWriteMode?TEXT_MODE_ID:BINARY_MODE_ID) << "\n";
}

bool PgHelpers::confirmTextReadMode(std::istream &src) {
    string readMode;
    src >> readMode;
    if (readMode != TEXT_MODE_ID && readMode != BINARY_MODE_ID) {
        fprintf(stderr, "Expected READ MODE id (not: %s)\n", readMode.c_str());
        exit(EXIT_FAILURE);
    }
    char check = src.get();
    return readMode == TEXT_MODE_ID;
}

template<>
void PgHelpers::writeValue(std::ostream &dest, const uint8_t value, bool plainTextWriteMode) {
    if (plainTextWriteMode)
        dest << (uint16_t) value << endl;
    else
        dest.write((char *) &value, sizeof(uint8_t));
}

template<>
void PgHelpers::readValue(std::istream &src, uint8_t &value, bool plainTextReadMode) {
    if (plainTextReadMode) {
        uint16_t temp;
        src >> temp;
        value = (uint8_t) temp;
    } else
        src.read((char *) &value, sizeof(uint8_t));
}

void PgHelpers::writeUIntByteFrugal(std::ostream &dest, uint64_t value) {
    while (value >= 128) {
        uint8_t yByte = 128 + (value % 128);
        dest.write((char *) &yByte, sizeof(uint8_t));
        value = value / 128;
    }
    uint8_t yByte = value;
    dest.write((char *) &yByte, sizeof(uint8_t));
}

bool PgHelpers::bytePerReadLengthMode = false;

void PgHelpers::readReadLengthValue(std::istream &src, uint16_t &value, bool plainTextReadMode) {
    if (bytePerReadLengthMode)
        readValue<uint8_t>(src, (uint8_t&) value, plainTextReadMode);
    else
        readValue<uint16_t>(src, value, plainTextReadMode);
}

void PgHelpers::writeReadLengthValue(std::ostream &dest, const uint16_t value) {
    if (bytePerReadLengthMode)
        writeValue<uint8_t>(dest, (uint8_t) value, plainTextWriteMode);
    else
        writeValue<uint16_t>(dest, value, plainTextWriteMode);
}

string PgHelpers::toString(unsigned long long value) {
        std::ostringstream oss;
        oss << value;
        return oss.str();
};

string PgHelpers::toMB(unsigned long long value, unsigned char decimalPlaces) {
    std::ostringstream oss;
    int power = 1000000;
    oss << value / power;
    if (decimalPlaces > 0) {
        oss << ".";
        for(int i = 0; i < decimalPlaces; i++) 
            oss << ((value / (power /= 10)) % 10);
        
    }
    return oss.str();
}

string PgHelpers::toString(long double value, unsigned char decimalPlaces) {
    std::ostringstream oss;
    oss << (long long int) value;
    if (decimalPlaces > 0) {
        oss << ".";
        int power = 1000000;        
        unsigned long long decimals = (value - (long long int) value) * power;
        for(int i = 0; i < decimalPlaces; i++) 
            oss << ((decimals / (power /= 10)) % 10);
    }
    return oss.str();
}

unsigned long long int PgHelpers::powuint(unsigned long long int base, int exp)
{
    if (exp == 0) return 1;
    if (exp == 1) return base;

    unsigned long long int tmp = PgHelpers::powuint(base, exp/2);
    if (exp%2 == 0) return tmp * tmp;
        else return base * tmp * tmp;
}

struct LUT
{
    char complementsLut[256];
    char sym2val[256];
    char val2sym[5];
    float qualityLut[133];

    void reorderSymAndVal(const char* basesOrdered) {
        for(char i = 0; i < 5; i++)
            val2sym[i] = basesOrdered[i];
        for(char i = 0; i < 5; i++)
            sym2val[val2sym[i]] = i;
    }

    LUT() {
        memset(complementsLut, 0, 256);
        complementsLut['A'] = 'T'; complementsLut['a'] = 'T';
        complementsLut['C'] = 'G'; complementsLut['c'] = 'G';
        complementsLut['G'] = 'C'; complementsLut['g'] = 'C';
        complementsLut['T'] = 'A'; complementsLut['t'] = 'A';
        complementsLut['N'] = 'N'; complementsLut['n'] = 'N';
        complementsLut['U'] = 'A'; complementsLut['u'] = 'A';
        complementsLut['Y'] = 'R'; complementsLut['y'] = 'R';
        complementsLut['R'] = 'Y'; complementsLut['r'] = 'Y';
        complementsLut['K'] = 'M'; complementsLut['k'] = 'M';
        complementsLut['M'] = 'K'; complementsLut['m'] = 'K';
        complementsLut['B'] = 'V'; complementsLut['b'] = 'V';
        complementsLut['D'] = 'H'; complementsLut['d'] = 'H';
        complementsLut['H'] = 'D'; complementsLut['h'] = 'D';
        complementsLut['V'] = 'B'; complementsLut['v'] = 'B';
        val2sym[0] = 'A';
        val2sym[1] = 'C';
        val2sym[2] = 'G';
        val2sym[3] = 'T';
        val2sym[4] = 'N';
        memset(sym2val, -1, 256);
        reorderSymAndVal("ACGTN");
        qualityLut[33 + 0] = 0;
        qualityLut[33 + 1] = 0.2056717652757185;
        qualityLut[33 + 2] = 0.36904265551980675;
        qualityLut[33 + 3] = 0.49881276637272776;
        qualityLut[33 + 4] = 0.6018928294465028;
        qualityLut[33 + 5] = 0.683772233983162;
        qualityLut[33 + 6] = 0.748811356849042;
        qualityLut[33 + 7] = 0.800473768503112;
        qualityLut[33 + 8] = 0.8415106807538887;
        qualityLut[33 + 9] = 0.8741074588205833;
        qualityLut[33 + 10] = 0.9;
        qualityLut[33 + 11] = 0.9205671765275718;
        qualityLut[33 + 12] = 0.9369042655519807;
        qualityLut[33 + 13] = 0.9498812766372727;
        qualityLut[33 + 14] = 0.9601892829446502;
        qualityLut[33 + 15] = 0.9683772233983162;
        qualityLut[33 + 16] = 0.9748811356849042;
        qualityLut[33 + 17] = 0.9800473768503112;
        qualityLut[33 + 18] = 0.9841510680753889;
        qualityLut[33 + 19] = 0.9874107458820583;
        qualityLut[33 + 20] = 0.99;
        qualityLut[33 + 21] = 0.9920567176527572;
        qualityLut[33 + 22] = 0.993690426555198;
        qualityLut[33 + 23] = 0.9949881276637272;
        qualityLut[33 + 24] = 0.996018928294465;
        qualityLut[33 + 25] = 0.9968377223398316;
        qualityLut[33 + 26] = 0.9974881135684904;
        qualityLut[33 + 27] = 0.9980047376850312;
        qualityLut[33 + 28] = 0.9984151068075389;
        qualityLut[33 + 29] = 0.9987410745882058;
        qualityLut[33 + 30] = 0.999;
        qualityLut[33 + 31] = 0.9992056717652757;
        qualityLut[33 + 32] = 0.9993690426555198;
        qualityLut[33 + 33] = 0.9994988127663728;
        qualityLut[33 + 34] = 0.9996018928294464;
        qualityLut[33 + 35] = 0.9996837722339832;
        qualityLut[33 + 36] = 0.999748811356849;
        qualityLut[33 + 37] = 0.9998004737685031;
        qualityLut[33 + 38] = 0.9998415106807539;
        qualityLut[33 + 39] = 0.9998741074588205;
        qualityLut[33 + 40] = 0.9999;
        for(int i = 41; i < 100; i++)
            qualityLut[33 + i] = 1;
    }
} lutInstance;

char* complementsLUT = lutInstance.complementsLut;
char* sym2val = lutInstance.sym2val;
char* val2sym = lutInstance.val2sym;
float* qualityLut = lutInstance.qualityLut;

void PgHelpers::reorderSymAndVal(const char* basesOrdered) {
    lutInstance.reorderSymAndVal(basesOrdered);
}

uint8_t PgHelpers::symbol2value(char symbol) {
    return sym2val[symbol];
}

char PgHelpers::value2symbol(uint8_t value) {
    return val2sym[value];
}

uint8_t PgHelpers::mismatch2code(char actual, char mismatch) {
    uint8_t actualValue = symbol2value(actual);
    uint8_t mismatchValue = symbol2value(mismatch);
    return mismatchValue - (mismatchValue > actualValue?1:0);
}

char PgHelpers::code2mismatch(char actual, uint8_t code) {
    uint8_t actualValue = symbol2value(actual);
    return value2symbol(code < actualValue?code:(code+1));
}

uint8_t PgHelpers::mismatch2CxtCode(char actual, char mismatch) {
    uint8_t actualValue = symbol2value(actual);
    uint8_t mismatchValue = symbol2value(mismatch);
    return (actualValue << 4) + mismatchValue;
}

uint8_t PgHelpers::cxtCode2ActualValue(uint8_t code) {
    uint8_t actualValue = code >> 4;
    return actualValue;
}

uint8_t PgHelpers::cxtCode2MismatchValue(uint8_t code) {
    uint8_t mismatchValue = code & 0x0F;
    return mismatchValue;
}

char PgHelpers::cxtCode2Mismatch(uint8_t code) {
    uint8_t mismatchValue = cxtCode2MismatchValue(code);
    return value2symbol(mismatchValue);
}

char PgHelpers::reverseComplement(char symbol) {
    return complementsLUT[symbol];
}

void PgHelpers::reverseComplementInPlace(char* start, const std::size_t N) {
    char* left = start - 1;
    char* right = start + N;
    while (--right > ++left) {
        char tmp = complementsLUT[*left];
        *left = complementsLUT[*right];
        *right = tmp;
    }
    if (left == right)
        *left = complementsLUT[*left];
}

string PgHelpers::reverseComplement(string kmer) {
    size_t kmer_length = kmer.size();
    string revcomp;
    revcomp.resize(kmer_length);
    size_t j = kmer_length;
    for(size_t i = 0; i < kmer_length; i++)
        revcomp[--j] = complementsLUT[kmer[i]];
    return revcomp;
}

void PgHelpers::reverseComplementInPlace(string &kmer) {
    reverseComplementInPlace((char*) kmer.data(), kmer.length());
}

double PgHelpers::qualityScore2approxCorrectProb(string& quality) {
    double val = 1;
    for (char q : quality) {
        switch (q) {
            case 33:case 34:case 35:case 36: return 0;
            case 37: val *= 0.6018928294465028; break;
            case 38: val *= 0.683772233983162; break;
            case 39: val *= 0.748811356849042; break;
            case 40: val *= 0.800473768503112; break;
            case 41: val *= 0.8415106807538887; break;
            case 42: val *= 0.8741074588205833; break;
            case 43: val *= 0.9; break;
            case 44: val *= 0.9205671765275718; break;
            case 45: val *= 0.9369042655519807; break;
            case 46: val *= 0.9498812766372727; break;
            case 47: val *= 0.9601892829446502; break;
            case 48: val *= 0.9683772233983162; break;
            case 49: val *= 0.9748811356849042; break;
            case 50: val *= 0.9800473768503112; break;
            case 51: val *= 0.9841510680753889; break;
            case 52: val *= 0.9874107458820583; break;
            case 53: val *= 0.99; break;
            case 54: val *= 0.9920567176527572; break;
            case 55: val *= 0.993690426555198; break;
            case 56: val *= 0.9949881276637272; break;
            case 57: val *= 0.996018928294465; break;
            case 58: val *= 0.9968377223398316; break;
            case 59: val *= 0.9974881135684904; break;
            case 60: val *= 0.9980047376850312; break;
            case 61: val *= 0.9984151068075389; break;
            case 62: val *= 0.9987410745882058; break;
            case 63: val *= 0.999; break;
            case 64: val *= 0.9992056717652757; break;
            case 65: val *= 0.9993690426555198; break;
            case 66: val *= 0.9994988127663728; break;
            case 67: val *= 0.9996018928294464; break;
            case 68: val *= 0.9996837722339832; break;
            default: ;
        }
    }
    return pow(val, 1.0/quality.length());
}

double PgHelpers::qualityScore2correctProbArithAvg(string& quality, int fraction, bool onlyRightFraction) {
    double val1 = 0, val2 = 0;
    int fractionLength = quality.length() / fraction;
    int rightLength = fractionLength;
    if (!onlyRightFraction) {
        int leftLength = fractionLength / 2;
        rightLength = fractionLength - leftLength;
        int i = 0;
        for (; i < (leftLength / 2) * 2; i += 2) {
            val1 += qualityLut[quality[i]];
            val2 += qualityLut[quality[i + 1]];
        }
        for (; i < leftLength; i++)
            val1 += qualityLut[quality[i]];
    }
    int i = quality.length() - rightLength;
    for (; i < (quality.length() / 2) * 2; i += 2) {
        val1 += qualityLut[quality[i]];
        val2 += qualityLut[quality[i + 1]];
    }
    for (; i < quality.length(); i++)
        val1 += qualityLut[quality[i]];
    return (val1 + val2) / fractionLength;
}


double PgHelpers::qualityScore2correctProb(string& quality) {
    double val = 1;
    for (char q : quality) {
        switch (q) {
            case 33: return 0;
            case 34: val *= 0.2056717652757185; break;
            case 35: val *= 0.36904265551980675; break;
            case 36: val *= 0.49881276637272776; break;
            case 37: val *= 0.6018928294465028; break;
            case 38: val *= 0.683772233983162; break;
            case 39: val *= 0.748811356849042; break;
            case 40: val *= 0.800473768503112; break;
            case 41: val *= 0.8415106807538887; break;
            case 42: val *= 0.8741074588205833; break;
            case 43: val *= 0.9; break;
            case 44: val *= 0.9205671765275718; break;
            case 45: val *= 0.9369042655519807; break;
            case 46: val *= 0.9498812766372727; break;
            case 47: val *= 0.9601892829446502; break;
            case 48: val *= 0.9683772233983162; break;
            case 49: val *= 0.9748811356849042; break;
            case 50: val *= 0.9800473768503112; break;
            case 51: val *= 0.9841510680753889; break;
            case 52: val *= 0.9874107458820583; break;
            case 53: val *= 0.99; break;
            case 54: val *= 0.9920567176527572; break;
            case 55: val *= 0.993690426555198; break;
            case 56: val *= 0.9949881276637272; break;
            case 57: val *= 0.996018928294465; break;
            case 58: val *= 0.9968377223398316; break;
            case 59: val *= 0.9974881135684904; break;
            case 60: val *= 0.9980047376850312; break;
            case 61: val *= 0.9984151068075389; break;
            case 62: val *= 0.9987410745882058; break;
            case 63: val *= 0.999; break;
            case 64: val *= 0.9992056717652757; break;
            case 65: val *= 0.9993690426555198; break;
            case 66: val *= 0.9994988127663728; break;
            case 67: val *= 0.9996018928294464; break;
            case 68: val *= 0.9996837722339832; break;
            case 69: val *= 0.999748811356849; break;
            case 70: val *= 0.9998004737685031; break;
            case 71: val *= 0.9998415106807539; break;
            case 72: val *= 0.9998741074588205; break;
            case 73: val *= 0.9999; break;
            default: val *= 1;
        }
    }
    return pow(val, 1.0/quality.length());
}

int PgHelpers::readsSufPreCmp(const char* suffixPart, const char* prefixRead) {
    while (*suffixPart) {
        if (*suffixPart > *prefixRead)
            return 1;
        if (*suffixPart++ < *prefixRead++)
            return -1;
    }
    return 0;
}

int PgHelpers::strcmplcp(const char* lStrPtr, const char* rStrPtr, int length) {
    
    int i = 0;
    while (length - i >= 4) {
        int cmp = bswap_32(*(uint32_t*) lStrPtr) - bswap_32(*(uint32_t*) rStrPtr);
        if (cmp != 0)
            break;
        lStrPtr += 4;
        rStrPtr += 4;
        i += 4;
    }

    while (i < length) {
        i++;
        int cmp = *(unsigned char*)(lStrPtr++) - *(unsigned char*)(rStrPtr++);
        if (cmp > 0)
            return 1;
        if (cmp < 0)
            return -1;
    }

    return 0;

}
