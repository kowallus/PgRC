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

double PgSAHelpers::combinedQuality(string quality) {
    return 0;
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
