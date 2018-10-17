#include "helper.h"

#include "byteswap.h"

// TIME

clock_t checkpoint;

bool PgSAHelpers::plainTextWriteMode = false;
bool PgSAHelpers::plainTextReadMode = false;

void PgSAHelpers::clock_checkpoint() {
    checkpoint = clock();
//    cout << "Clock reset!\n";
}

unsigned long long int PgSAHelpers::clock_millis() {
    return (clock() - checkpoint) * (unsigned long long int) 1000 / CLOCKS_PER_SEC;
}

const size_t chunkSize = 10000000;

void* PgSAHelpers::readArray(std::istream& in, size_t arraySizeInBytes) {

    if (in) {
        size_t length = arraySizeInBytes;
        char* destArray = new char[length];
        char* ptr = destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);

        return (void*) destArray;
    } else
        throw(errno);
}

void PgSAHelpers::writeArray(std::ostream& out, void* srcArray, size_t arraySize) {
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
    
    cout << "Written " << arraySize << " bytes\n";
}

void* PgSAHelpers::readWholeArray(std::istream& in, size_t& arraySizeInBytes) {

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

void* PgSAHelpers::readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes) {
    clock_checkpoint();

    std::ifstream in(srcFile.c_str(), std::ifstream::binary);

    void* destArray = PgSAHelpers::readWholeArray(in, arraySizeInBytes);

    cout << "Read " << arraySizeInBytes << " bytes from " << srcFile << " in " << clock_millis() << " msec \n";

    return destArray;
}

void PgSAHelpers::writeArrayToFile(string destFile, void* srcArray, size_t arraySize) {
    clock_checkpoint();

    std::ofstream out(destFile.c_str(), std::ios::out | std::ios::binary);

    PgSAHelpers::writeArray(out, srcArray, arraySize);

    cout << "Write " << arraySize << " bytes to " << destFile << " in " << clock_millis() << " msec \n";
}

string PgSAHelpers::toString(unsigned long long value) {
        std::ostringstream oss;
        oss << value;
        return oss.str();
};

string PgSAHelpers::toMB(unsigned long long value, unsigned char decimalPlaces) {
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

string PgSAHelpers::toString(long double value, unsigned char decimalPlaces) {
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

unsigned long long int PgSAHelpers::powuint(unsigned long long int base, int exp)
{
    if (exp == 0) return 1;
    if (exp == 1) return base;

    unsigned long long int tmp = PgSAHelpers::powuint(base, exp/2);
    if (exp%2 == 0) return tmp * tmp;
        else return base * tmp * tmp;
}

uint8_t PgSAHelpers::symbol2value(char symbol) {
    switch (symbol) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
    }
    return -1;
}

char PgSAHelpers::value2symbol(uint8_t value) {
    switch (value) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        case 4: return 'N';
    }
    return '?';
}

uint8_t PgSAHelpers::mismatch2code(char actual, char mismatch) {
    uint8_t actualValue = symbol2value(actual);
    uint8_t mismatchValue = symbol2value(mismatch);
    return mismatchValue > actualValue?(mismatchValue - 1):mismatchValue;
}

char PgSAHelpers::code2mismatch(char actual, uint8_t code) {
    uint8_t actualValue = symbol2value(actual);
    return value2symbol(code < actualValue?code:(code+1));
}

string PgSAHelpers::reverseComplement(string kmer) {
    size_t kmer_length = kmer.size();
    string revcomp;
    revcomp.resize(kmer_length);
    size_t j = kmer_length;
    for(size_t i = 0; i < kmer_length; i++) {
        switch(kmer[i]) {
            case 'A': revcomp[--j] = 'T'; break;
            case 'T': revcomp[--j] = 'A'; break;
            case 'U': revcomp[--j] = 'A'; break;
            case 'G': revcomp[--j] = 'C'; break;
            case 'C': revcomp[--j] = 'G'; break;
            case 'Y': revcomp[--j] = 'R'; break;
            case 'R': revcomp[--j] = 'Y'; break;
            case 'K': revcomp[--j] = 'M'; break;
            case 'M': revcomp[--j] = 'K'; break;
            case 'B': revcomp[--j] = 'V'; break;
            case 'D': revcomp[--j] = 'H'; break;
            case 'H': revcomp[--j] = 'D'; break;
            case 'V': revcomp[--j] = 'B'; break;
            case 'N': revcomp[--j] = 'N'; break;
            default: revcomp[--j] = kmer[i];
                cout << "WARNING: Unsupported reverse compliment: " << kmer[i] << "\n";
                break;
        }
    }
    return revcomp;
}

double PgSAHelpers::qualityScore2approxCorrectProb(string quality) {
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
    return val;
}

double PgSAHelpers::qualityScore2correctProb(string quality) {
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
    return val;
}

int PgSAHelpers::readsSufPreCmp(const char* suffixPart, const char* prefixRead) {
    while (*suffixPart) {
        if (*suffixPart > *prefixRead)
            return 1;
        if (*suffixPart++ < *prefixRead++)
            return -1;
    }
    return 0;
}

int PgSAHelpers::strcmplcp(const char* lStrPtr, const char* rStrPtr, int length) {
    
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
