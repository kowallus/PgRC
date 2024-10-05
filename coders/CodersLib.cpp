#include <assert.h>
#include "CodersLib.h"
#include "LzmaCoder.h"
#include "PpmdCoder.h"
#include "RangeCoder.h"
#include "FSECoder.h"
#include "VarLenDNACoder.h"
#include <omp.h>
#include <memory>
#include "PropsLibrary.h"

#ifdef DEVELOPER_BUILD
bool dump_after_decompression = false;
int dump_after_decompression_counter = 1;
string dump_after_decompression_prefix;

struct CompressorAutoSelector {

    int level = 0;
    vector<unique_ptr<CoderProps>> coderPropsList;

    CompressorAutoSelector() {
        coderPropsList.push_back(getDefaultFSECoderProps());
        coderPropsList.push_back(getDefaultRangeCoderProps(256, 1));
        coderPropsList.push_back(getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 2));
        coderPropsList.push_back(getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t));
        coderPropsList.push_back(getDefaultRangeCoderProps(256, 2));
        coderPropsList.push_back(getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 3));
        coderPropsList.push_back(getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_16_t));
        coderPropsList.push_back(getDefaultRangeCoderProps(256, 4));
        coderPropsList.push_back(getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 5));
        coderPropsList.push_back(getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_32_t));
        coderPropsList.push_back(getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 8));
        coderPropsList.push_back(getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 13));
    }

    void setAutoSelectorLevel(int level) {
        this->level = level;
        if (level < 0)
            this->level = 0;
        if (level > coderPropsList.size())
            this->level = coderPropsList.size();
    }

    unique_ptr<CoderProps> getSelectorCoder(CoderProps* primaryCoder) {
        if (coderPropsList.size() <= 0 || level <= 0) {
            fprintf(stderr, "CompressorAutoSelector requires at least 2 coders, supplied 1.\n");
            exit(EXIT_FAILURE);
        }

        auto propsIter = coderPropsList.begin() + level;
        CoderProps* res = (--propsIter)->get();
        while (propsIter-- != coderPropsList.begin()) {
            res = new SelectorCoderProps(propsIter->get(), res, 1);
        }
        return unique_ptr<CoderProps>(new SelectorCoderProps(primaryCoder, res, 1));
    }

}  casInstance;

int& auto_selector_level = casInstance.level;
vector<unique_ptr<CoderProps>>& auto_selector_props_list = casInstance.coderPropsList;

void setAutoSelectorLevel(int level) {
    casInstance.setAutoSelectorLevel(level);
}

unique_ptr<CoderProps> getSelectorCoder(CoderProps* primaryCoder) {
    return casInstance.getSelectorCoder(primaryCoder);
}
#endif

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

using namespace PgHelpers;

unsigned char* Compress(size_t &destLen, const unsigned char *src, size_t srcLen, CoderProps* props,
        double estimated_compression, std::ostream* logout, bool disable_auto_selector) {
#if DEVELOPER_BUILD
    unique_ptr<CoderProps> scPropsPtr;
    if (!disable_auto_selector && auto_selector_level > 0 && props->getCoderType() != VARLEN_DNA_CODER
        && props->getCoderType() != PARALLEL_BLOCKS_CODER_TYPE && props->getCoderType() != COMPOUND_CODER_TYPE
            && props->getCoderType() != SELECTOR_CODER_TYPE) {
        scPropsPtr = getSelectorCoder(props);
        props = scPropsPtr.get();
    }
#endif
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    unsigned char* dest = nullptr;
    unsigned char* compSeq = nullptr;
    CompoundCoderProps *ccProps;
    SelectorCoderProps* scProps;
    int res = 0;
    switch (props->getCoderType()) {
        case LZMA_CODER:
            res = LzmaCompress(dest, destLen, src, srcLen,
                               ((LzmaCoderProps *) props)->getProps(), estimated_compression);
            break;
        case PPMD7_CODER:
            res = Ppmd7Compress(dest, destLen, src, srcLen, (PpmdCoderProps *) props,
                                estimated_compression);
            break;
        case RANGE_CODER:
            res = RangeCompress(dest, destLen, src, srcLen, (RangeCoderProps *) props,
                                estimated_compression);
            break;
        case FSE_HUF_CODER:
            res = FSECompress(dest, destLen, src, srcLen, (FSECoderProps *) props,
                              estimated_compression);
            break;
        case VARLEN_DNA_CODER:
            res = VarLenDNACoder::Compress(dest, destLen, src, srcLen,
                                           (VarLenDNACoderProps *) props);
            estimated_compression = VarLenDNACoder::COMPRESSION_ESTIMATION;
            break;
        case PARALLEL_BLOCKS_CODER_TYPE:
            res = parallelBlocksCompress(dest, destLen, src, srcLen,
                                         (ParallelBlocksCoderProps *) props, estimated_compression, logout);
            break;
        case COMPOUND_CODER_TYPE:
            *logout << "\n\t\t";
            ccProps = (CompoundCoderProps *) props;
            compSeq = Compress(destLen, (const unsigned char *) src, srcLen,
                               ccProps->primaryCoder, estimated_compression, logout);
            if (destLen >= srcLen) {
                destLen = srcLen;
                memcpy(compSeq, src, srcLen);
            }
            ccProps->compLen = destLen;
            *logout << "\t\t";
            dest = Compress(destLen, (const unsigned char *) compSeq, destLen,
                            ccProps->secondaryCoder, estimated_compression, logout);
            if (destLen >= ccProps->compLen) {
                destLen = ccProps->compLen;
                swap(dest, compSeq);
            }
            *logout << "\t\t";
            delete[] compSeq;
            res = SZ_OK;
            break;
        case SELECTOR_CODER_TYPE:
            scProps = (SelectorCoderProps*) props;
            if (scProps->selected) {
                dest = Compress(destLen, (const unsigned char *) src, srcLen,
                                scProps->getSelectedCoder(), estimated_compression, logout, true);
            } else {
                size_t compALen, compBLen;
                size_t srcProbeLen = scProps->getProbeLength(srcLen);
                stringstream aLogout, bLogout;
                unsigned char *compASeq = Compress(compALen, (const unsigned char *) src, srcProbeLen, scProps->coderA,
                                    estimated_compression, &aLogout, true);
                unsigned char *compBSeq = Compress(compBLen, (const unsigned char *) src, srcProbeLen, scProps->coderB,
                                    estimated_compression, &bLogout, true);
                bool selectA = compALen <= compBLen;
                scProps->selectCoder(selectA);
                if (srcProbeLen == srcLen) {
                    dest = selectA ? compASeq : compBSeq;
                    destLen = selectA ? compALen : compBLen;
                    delete[] (selectA ? compBSeq : compASeq);
                    *logout << (selectA ? aLogout : bLogout).str();
                } else {
                    delete[] compASeq;
                    delete[] compBSeq;
                    dest = Compress(destLen, (const unsigned char *) src, srcLen,
                                    scProps->getSelectedCoder(), estimated_compression, logout, true);
                }
            }
            return dest;
        case LZMA2_CODER:
        default:
            fprintf(stderr, "Unsupported coder type: %d.\n", props->getCoderType());
            exit(EXIT_FAILURE);
    }

    if (res != SZ_OK) {
        fprintf(stderr, "Error during compression (code: %d).\n", res);
        fprintf(stderr,"%s\tsrcLen: %zu\tdestLen: %zu\tcoder: %d\n",
                props->log().c_str(), srcLen, destLen, (int) props->getCoderType());
        exit(EXIT_FAILURE);
    }
    if (destLen < srcLen) {
        const double ratio = ((double) destLen) / srcLen;
        *logout << props->log() << " ... ";
        *logout << "compressed " << srcLen << " bytes to " << destLen << " bytes (ratio "
                << PgHelpers::toString(ratio, 3) << " vs estimated "
                << PgHelpers::toString(estimated_compression, 3) << ") in "
                << PgHelpers::time_millis(start_t) << " msec." << endl;
        if (ratio > estimated_compression)
            *logout << "WARNING: compression ratio " << PgHelpers::toString(ratio / estimated_compression, 5)
                    << " times greater than estimation." << endl;
    } else {
        *logout << " ignored" << props->log() << " due ratio >=1 ... " << srcLen << " bytes disregarding compression in "
                << PgHelpers::time_millis(start_t) << " msec." << endl;
    }
    return dest;
}

MY_STDAPI NoCoderUncompress(unsigned char *dest, size_t *destLen, unsigned char *src, size_t srcLen, ostream* logout) {
    *destLen = srcLen;
    *logout << "... raw ... ";
    memcpy(dest, (const void *) src, *destLen);
    return SZ_OK;
}

void Uncompress(unsigned char* dest, size_t destLen, istream &src, size_t srcLen, uint8_t coder_type,
        std::ostream* logout) {
    switch (coder_type) {
    default:
        string srcString;
        srcString.resize(srcLen);
        PgHelpers:readArray(src, (void*) srcString.data(), srcLen);
        Uncompress(dest, destLen, (unsigned char*) srcString.data(), srcLen, coder_type, logout);
    }
}

void Uncompress(unsigned char* dest, size_t destLen, unsigned char* src, size_t srcLen, uint8_t coder_type,
                std::ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int res = SZ_OK;
    size_t outLen = destLen;
    switch (coder_type) {
        case VARLEN_DNA_CODER:
            res = PgHelpers::VarLenDNACoder::Uncompress(dest, &outLen, src, srcLen);
            break;
        case NO_CODER:
            res = NoCoderUncompress(dest, &outLen, src, srcLen, logout);
            break;
        case LZMA_CODER:
            res = LzmaUncompress(dest, &outLen, src, srcLen, logout);
            break;
        case PPMD7_CODER:
            res = PpmdUncompress(dest, &outLen, src, srcLen, logout);
            break;
        case RANGE_CODER:
            res = RangeUncompress(dest, &outLen, src, srcLen, logout);
            break;
        case FSE_HUF_CODER:
            res = FSEUncompress(dest, &outLen, src, srcLen, logout);
            break;
        case PARALLEL_BLOCKS_CODER_TYPE:
            res = parallelBlocksDecompress(dest, &outLen, src, logout);
            break;
        case LZMA2_CODER:
        default:
            fprintf(stderr, "Unsupported coder type: %d.\n", coder_type);
            exit(EXIT_FAILURE);
    }
    assert(outLen == destLen);

    if (res != SZ_OK) {
        fprintf(stderr, "Error during decompression (code: %d).\n", res);
        fprintf(stderr, "srcLen: %zu\tdestLen: %zu\tcoder: %d\n", srcLen, destLen, (int) coder_type);
        exit(EXIT_FAILURE);
    }
    *logout << "uncompressed " << srcLen << " bytes to " << destLen << " bytes in "
                         << PgHelpers::time_millis(start_t) << " msec." << endl;
}

void writeCompressed(ostream &dest, const char *src, size_t srcLen, CoderProps* props, double estimated_compression) {
    if (srcLen == 0) {
        PgHelpers::writeValue<uint64_t>(dest, 0);
        *PgHelpers::devout << "skipped compression (0 bytes)." << endl;
        return;
    }
    size_t compLen = 0;
    unsigned char* compSeq = Compress(compLen, (const unsigned char*) src, srcLen, props, estimated_compression);
    writeHeader(dest, srcLen, compLen, props);
    if (compLen < srcLen) {
        PgHelpers::writeArray(dest, (void*) compSeq, compLen);
    } else {
        PgHelpers::writeArray(dest, (void *) src, srcLen);
    }
    delete[] compSeq;
}

void writeCompressed(ostream &dest, const string& srcStr, CoderProps* props, double estimated_compression) {
    writeCompressed(dest, srcStr.data(), srcStr.length(), props, estimated_compression);
}

void writeHeader(ostream &dest, size_t srcLen, size_t destLen, CoderProps *props) {
    PgHelpers::writeValue<uint64_t>(dest, srcLen);
    if (props->getCoderType() == COMPOUND_CODER_TYPE) {
        PgHelpers::writeValue<uint64_t>(dest,destLen + ((CompoundCoderProps*) props)->primaryCoder->getHeaderLen());
        PgHelpers::writeValue<uint8_t>(dest, COMPOUND_CODER_TYPE);
        PgHelpers::writeValue<uint64_t>(dest, ((CompoundCoderProps*) props)->compLen);
        CoderProps* primaryCoderProps = ((CompoundCoderProps*) props)->primaryCoder;
        primaryCoderProps = primaryCoderProps->getCoderType() == SELECTOR_CODER_TYPE ? ((SelectorCoderProps*) primaryCoderProps)->getSelectedCoder() : primaryCoderProps;
        PgHelpers::writeValue<uint8_t>(dest, ((CompoundCoderProps*) props)->compLen == srcLen ?
            NO_CODER : primaryCoderProps->getCoderType());
        writeHeader(dest, ((CompoundCoderProps*) props)->compLen, destLen, ((CompoundCoderProps*) props)->secondaryCoder);
    } else if (destLen >= srcLen) {
        PgHelpers::writeValue<uint64_t>(dest, srcLen);
        PgHelpers::writeValue<uint8_t>(dest, NO_CODER);
    } else {
        PgHelpers::writeValue<uint64_t>(dest, destLen);
        props = props->getCoderType() == SELECTOR_CODER_TYPE ? ((SelectorCoderProps*) props)->getSelectedCoder() : props;
        PgHelpers::writeValue<uint8_t>(dest, props->getCoderType());
    }
}

vector<string*> prefetchedStreams;
int prefetchedStreamsLimit = 0;
int prefetchedStreamsRead = 0;

bool areStreamsPrefetched() {
    return prefetchedStreamsRead < prefetchedStreamsLimit;
}


void readCompressed(istream &src, string& dest, ostream* logout) {
    if (areStreamsPrefetched()) {
        dest = std::move(*prefetchedStreams[prefetchedStreamsRead]);
        delete(prefetchedStreams[prefetchedStreamsRead++]);
        return;
    }
    uint64_t destLen = 0;
    uint64_t srcLen = 0;
    uint8_t coder_type = 0;
    PgHelpers::readValue<uint64_t>(src, destLen);
    dest.resize(destLen);
    if (destLen == 0)
        return;
    PgHelpers::readValue<uint64_t>(src, srcLen);
    PgHelpers::readValue<uint8_t>(src, coder_type);
    if (coder_type == COMPOUND_CODER_TYPE) {
        uint64_t compLen = 0;
        uint8_t primary_coder_type = 0;
        PgHelpers::readValue<uint64_t>(src, compLen);
        PgHelpers::readValue<uint8_t>(src, primary_coder_type);
        string component;
        *logout << "\t\t";
        readCompressed(src, component, logout);
        assert(compLen == component.length());
        *logout << "\t\t";
        Uncompress((unsigned char *) dest.data(), destLen, (unsigned char*) component.data(), compLen, primary_coder_type, logout);
    } else {
        Uncompress((unsigned char *) dest.data(), destLen, src, srcLen, coder_type, logout);
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10?"0":"");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                      (void*) dest.data(), destLen);
    }
#endif
}

void readCompressed(unsigned char *src, string& dest, ostream* logout) {
    uint64_t destLen = 0;
    uint64_t srcLen = 0;
    uint8_t coder_type = 0;
    PgHelpers::readValue<uint64_t>(src, destLen);
    dest.resize(destLen);
    if (destLen == 0)
        return;
    PgHelpers::readValue<uint64_t>(src, srcLen);
    PgHelpers::readValue<uint8_t>(src, coder_type);
    if (coder_type == COMPOUND_CODER_TYPE) {
        uint64_t compLen = 0;
        uint8_t primary_coder_type = 0;
        PgHelpers::readValue<uint64_t>(src, compLen);
        PgHelpers::readValue<uint8_t>(src, primary_coder_type);
        string component;
        *logout << "\t\t";
        readCompressed(src, component, logout);
        assert(compLen == component.length());
        *logout << "\t\t";
        Uncompress((unsigned char *) dest.data(), destLen, (unsigned char*) component.data(), compLen, primary_coder_type, logout);
    } else {
        Uncompress((unsigned char *) dest.data(), destLen, src, srcLen, coder_type, logout);
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10?"0":"");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                      (void*) dest.data(), destLen);
    }
#endif
}

double simpleUintCompressionEstimate(uint64_t dataMaxValue, uint64_t typeMaxValue) {
    const double dataBits = 64 - ((double) __builtin_clzl(dataMaxValue));
    const double typeBits = 64 - __builtin_clzl(typeMaxValue);
    return dataBits / typeBits;
}

int parallelBlocksCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
        ParallelBlocksCoderProps *props, double estimated_compression, ostream* logout) {
    props->prepare(srcLen);
    size_t blockSize = ((srcLen / props->numOfBlocks) / props->blockAlignment) * props->blockAlignment;
    stringstream destOut;
    PgHelpers::writeValue(destOut, props->numOfBlocks);
    vector<CompressionJob> cJobs;
    size_t offset = 0;
    for (int i = 0; i < props->numOfBlocks - 1; i++) {
        cJobs.push_back(CompressionJob("block " + toString(i + 1) + "... ", src + offset, blockSize,
                props->blocksCoder));
        offset += blockSize;
    }
    cJobs.push_back(CompressionJob("block " + toString(props->numOfBlocks) + "... ", src + offset, srcLen - offset,
                    props->blocksCoder));
    CompressionJob::writeCompressedCollectiveParallel(destOut, cJobs, &null_stream);
    destLen = destOut.tellp();
    dest = new unsigned char[destLen];
    destOut.seekg(0);
    destOut.read((char*) dest, destLen);

    return SZ_OK;
}

int parallelBlocksDecompress(unsigned char *dest, size_t *destLen, unsigned char* src, ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int numOfBlocks;
    PgHelpers::readValue<int>(src, numOfBlocks);
    *logout << "... parallel_blocks (no = " << numOfBlocks << ") of ";
    size_t offset = 0;
#pragma omp parallel
    {
#pragma omp single
        {
            for (int i = 0; i < numOfBlocks; i++) {
                unsigned char* destPtr = dest + offset;
                uint64_t blockLen = 0;
                uint64_t srcLen = 0;
                uint8_t coder_type = 0;
                unsigned char* srcString;
                PgHelpers::readValue<uint64_t>(src, blockLen);
                offset += blockLen;
                if (blockLen > 0) {
                    PgHelpers::readValue<uint64_t>(src, srcLen);
                    PgHelpers::readValue<uint8_t>(src, coder_type);
                    srcString = src;
                    src += srcLen;
                }
                if (blockLen == 0)
                    continue;
#pragma omp task
                {
                    ostringstream tmpout;
                    ostream* currentOut = i == 0?&tmpout:&null_stream;
                    Uncompress(destPtr, blockLen, srcString, srcLen, coder_type, currentOut);
                    if (i == 0) {
                        string log = tmpout.str();
                        *logout << log.substr(4, log.find("...", 4));
                    }
                }
            }
        }
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10 ? "0" : "");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                    (void *) dest, *destLen);
    }
#endif
    return SZ_OK;
}

CompressionJob::CompressionJob(string label, const unsigned char *src, size_t srcLen,
        CoderProps *props, double estimated_compression) : log(label), src(src), srcLen(srcLen),
        props(props), estimated_compression(estimated_compression) {}

CompressionJob::CompressionJob(string label, const string &src,
        CoderProps *props, double estimated_compression) : log(label), src((const unsigned char*) src.data()),
        srcLen(src.length()), props(props), estimated_compression(estimated_compression) { }

void CompressionJob::writeCompressedCollectiveParallel(ostream &dest, vector<CompressionJob> &cJobs, ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    *logout << "collective compression of streams..." << endl;
    vector<unsigned char*> compSeqs;
    vector<size_t> compLens;
    compSeqs.resize(cJobs.size());
    compLens.resize(cJobs.size(), 0);
#ifdef __APPLE__
    omp_set_max_active_levels(4);
#else
    omp_set_nested(true);
#endif
#pragma omp parallel for
    for(int i = 0; i < cJobs.size(); i++) {
        ostringstream localLogOut;
        compSeqs[i] = Compress(compLens[i], cJobs[i].src, cJobs[i].srcLen, cJobs[i].props,
                cJobs[i].estimated_compression, &localLogOut);
        cJobs[i].log.append(localLogOut.str());
    }
    for(int i = 0; i < cJobs.size(); i++) {
        if (cJobs[i].srcLen == 0) {
            PgHelpers::writeValue<uint64_t>(dest, 0);
            *PgHelpers::devout << "skipped compression (0 bytes)." << endl;
            continue;
        }
        writeHeader(dest, cJobs[i].srcLen, compLens[i], cJobs[i].props);
        *logout << "\t" << cJobs[i].log;
        if (compLens[i] < cJobs[i].srcLen)
            PgHelpers::writeArray(dest, (void *) compSeqs[i], compLens[i]);
        else
            PgHelpers::writeArray(dest, (void *) cJobs[i].src, cJobs[i].srcLen);

        delete[] compSeqs[i];
    }
    *logout << "collective compression finished in " << PgHelpers::time_millis(start_t) << " msec." << endl;
    if (logout != &null_stream && logout != PgHelpers::logout)
        *PgHelpers::logout << PgHelpers::time_millis(start_t) << "       \t";
}

void prefetchCompressedCollectiveParallel(istream &src, int streamsLimit) {
    for (int i = 0; i < streamsLimit; i++)
        prefetchedStreams.push_back(new string());
    readCompressedCollectiveParallel(src, prefetchedStreams);
    prefetchedStreamsLimit = prefetchedStreams.size();
}

void readCompressedCollectiveParallel(istream &src, vector<string*>& destStrings, int obsolete_pgrc_compound_coder) {
    int i = 0;
    if (areStreamsPrefetched()) {
        for (; i < destStrings.size(); i++) {
            if (areStreamsPrefetched()) {
                *destStrings[i] = std::move(*prefetchedStreams[prefetchedStreamsRead]);
                delete(prefetchedStreams[prefetchedStreamsRead++]);
                continue;
            }
        }
    }
    if (i == destStrings.size())
        return;
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    *PgHelpers::devout << "collective decompression of streams..." << endl;
    vector<ostringstream> logOuts;
    logOuts.resize(destStrings.size());
#ifdef __APPLE__
    omp_set_max_active_levels(4);
#else
    omp_set_nested(true);
#endif
#pragma omp parallel
    {
#pragma omp single
        {
            for (int i = 0; i < destStrings.size(); i++) {
                uint64_t destLen = 0;
                uint64_t srcLen = 0;
                uint64_t compLen = 0;
                uint8_t coder_type = 0;
                uint8_t primary_coder_type = 0;
                string srcString, component;
                if (prefetchedStreamsLimit == 0 && !prefetchedStreams.empty() && src.peek() == -1) {
                    for(int j = i; j < prefetchedStreams.size(); j++)
                        delete(prefetchedStreams[j]);
                    prefetchedStreams.resize(i);
                    assert(i == destStrings.size());
                    continue;
                }
                PgHelpers::readValue<uint64_t>(src, destLen);
                destStrings[i]->resize(destLen);
                if (destLen == 0)
                    continue;
                if (obsolete_pgrc_compound_coder == i) {
                    PgHelpers::readValue<uint64_t>(src, compLen);
                    PgHelpers::readValue<uint8_t>(src, coder_type);
                    assert(coder_type == COMPOUND_CODER_TYPE);
                    PgHelpers::readValue<uint8_t>(src, primary_coder_type);
                    readCompressed(src, component, &logOuts[i]);
                } else {
                    PgHelpers::readValue<uint64_t>(src, srcLen);
                    PgHelpers::readValue<uint8_t>(src, coder_type);
                    srcString.resize(srcLen);
                    if (coder_type == COMPOUND_CODER_TYPE) {
                        PgHelpers::readValue<uint64_t>(src, compLen);
                        PgHelpers::readValue<uint8_t>(src, primary_coder_type);
                    }
                    PgHelpers::readArray(src, (void *) srcString.data(), srcLen);
                }
#pragma omp task
                {
                    if (coder_type == COMPOUND_CODER_TYPE) {
                        logOuts[i] << "\t";
                        if (obsolete_pgrc_compound_coder != i)
                            readCompressed((unsigned char*) srcString.data(), component, &logOuts[i]);
                        assert(compLen == component.length());
                        logOuts[i] << "\t\t";
                        Uncompress((unsigned char *) destStrings[i]->data(), destLen, (unsigned char*) component.data(), compLen,
                                   primary_coder_type, &logOuts[i]);
                    } else
                        Uncompress((unsigned char *) destStrings[i]->data(), destLen, (unsigned char*) srcString.data(), srcLen, coder_type,
                                   &logOuts[i]);
                }
            }
        }
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression)
        for (int i = 0; i < destStrings.size(); i++) {
            string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10 ? "0" : "");
            PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                        (void *) destStrings[i]->data(), destStrings[i]->size());
        }
#endif
    for (int i = 0; i < destStrings.size(); i++) {
        *PgHelpers::devout << "\t" << logOuts[i].str();
    }
    *PgHelpers::devout << "collective decompression finished in " << PgHelpers::time_millis(start_t) << " msec."
                       << endl;
    if (PgHelpers::logout != PgHelpers::devout)
        *PgHelpers::logout << PgHelpers::time_millis(start_t) << "       \t";
}

