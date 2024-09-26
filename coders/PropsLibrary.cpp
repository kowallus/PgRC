#include "PropsLibrary.h"

#include "PpmdCoder.h"
#include "LzmaCoder.h"
#include "RangeCoder.h"
#include "FSECoder.h"
#include "VarLenDNACoder.h"


LzmaCoderProps* getDefaulLzmaCoderProps(int coder_level, int dataPeriodCode = -1) {
    int noOfThreads = PgHelpers::numberOfThreads > 1 ? 2 : 1;
    switch(coder_level) {
        case CODER_LEVEL_FAST:
            return new LzmaCoderProps(5, 1 << 24, 3, dataPeriodCode, dataPeriodCode, 6, -1, noOfThreads);
        case CODER_LEVEL_NORMAL:
            return new LzmaCoderProps(9, 3 << 29, 3, dataPeriodCode, dataPeriodCode, 128, -1, noOfThreads);
        case CODER_LEVEL_MAX:
            //            p->lc = 3; // test 4 for big files
            return new LzmaCoderProps(9, 3 << 29, 3, dataPeriodCode, dataPeriodCode, 273, -1, noOfThreads);
        default:
            fprintf(stderr, "Unsupported %d PgRC coding level for LZMA compression.\n", coder_level);
            exit(EXIT_FAILURE);
    }
}

PpmdCoderProps* getDefaultPpmdCoderProps(uint8_t coder_level, int order_param) {
    switch (coder_level) {
        case CODER_LEVEL_FAST:
            return new PpmdCoderProps((uint32_t) 16 << 20, 2);
        case CODER_LEVEL_NORMAL:
            return new PpmdCoderProps((uint32_t) 192 << 20, order_param > 2 ? (order_param / 2) + 1 : order_param);
        case CODER_LEVEL_MAX:
            return new PpmdCoderProps((uint32_t) 192 << 20, order_param);
        default:
            fprintf(stderr, "Unsupported %d PgRC coding level for Ppmd compression.\n", coder_level);
            exit(EXIT_FAILURE);
    }
}

RangeCoderProps* getRangeCoderProps(uint16_t no_of_symbols, uint8_t bytes_period) {
        return new RangeCoderProps(no_of_symbols, bytes_period);
}

unique_ptr<CoderProps> getDefaultFSECoderProps(uint8_t tableLog, uint8_t maxSymbolValue) {
    return unique_ptr<CoderProps>(new FSECoderProps(false, maxSymbolValue, tableLog));
}

unique_ptr<CoderProps> getDefaultHufCoderProps(uint8_t tableLog, uint8_t maxSymbolValue) {
    return unique_ptr<CoderProps>(new FSECoderProps(true, maxSymbolValue, tableLog));
}

unique_ptr<CoderProps> getDefaultCoderProps(uint8_t coder_type, uint8_t coder_level, int coder_param) {
    switch(coder_type) {
        case LZMA_CODER:
            return unique_ptr<CoderProps>(getDefaulLzmaCoderProps(coder_level, coder_param));
        case PPMD7_CODER:
            return unique_ptr<CoderProps>(getDefaultPpmdCoderProps(coder_level, coder_param));
        case RANGE_CODER:
            return unique_ptr<CoderProps>(getRangeCoderProps(coder_param, 1));
        case VARLEN_DNA_CODER:
            return unique_ptr<CoderProps>(
                    new PgHelpers::VarLenDNACoderProps(PgHelpers::VarLenDNACoder::AG_EXTENDED_CODES_ID));
        case LZMA2_CODER:
        default:
            fprintf(stderr, "Unsupported coder type: %d.\n", (int) coder_type);
            exit(EXIT_FAILURE);

    }
}

unique_ptr<CoderProps> getCompoundCoderProps(CoderProps* firstCoderProps, CoderProps* secondCoderProps) {
    return unique_ptr<CoderProps>(new CompoundCoderProps(firstCoderProps, secondCoderProps));
}

unique_ptr<CoderProps> getSelectorCoderProps(std::initializer_list<CoderProps*> coderPropsList, float probeFraction,
                                             size_t minProbeSize) {
    if (coderPropsList.size() <= 1) {
        fprintf(stderr, "Selector coder requires at least 2 coders, supplied %d.\n", (int) coderPropsList.size());
        exit(EXIT_FAILURE);
    }
    auto propsIter = coderPropsList.end();
    CoderProps* res = *(--propsIter);
    while (--propsIter != coderPropsList.begin()) {
        res = new SelectorCoderProps(*propsIter, res, 1, minProbeSize);
    }
    return unique_ptr<CoderProps>(new SelectorCoderProps(*propsIter, res, probeFraction, minProbeSize));
}


unique_ptr<CoderProps> getReadsPositionsCoderProps(uint8_t coder_level, uint8_t lzma_pos_dataperiod_param) {
    int noOfThreads = PgHelpers::numberOfThreads > 1 ? 2 : 1;
    int lc = 8;
    switch(coder_level) {
        case CODER_LEVEL_FAST:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 1 << 20, lc,
                    lzma_pos_dataperiod_param, lzma_pos_dataperiod_param, 5, 1, noOfThreads));
        case CODER_LEVEL_NORMAL:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 8 << 20, lc,
                    lzma_pos_dataperiod_param, lzma_pos_dataperiod_param, 16, 1, noOfThreads));
        case CODER_LEVEL_MAX:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 64 << 20, lc,
                    lzma_pos_dataperiod_param, lzma_pos_dataperiod_param, 16, 1, noOfThreads));
        default:
            fprintf(stderr, "Unsupported coder level: %d.\n", (int) coder_level);
            exit(EXIT_FAILURE);
    }
}

unique_ptr<CoderProps> getVarLenEncodedPgCoderProps(uint8_t coder_level) {
    int noOfThreads = PgHelpers::numberOfThreads > 1 ? 2 : 1;
    switch(coder_level) {
        case CODER_LEVEL_FAST:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 1 << 24, 8,
                    LZMA_DATAPERIODCODE_8_t, LZMA_DATAPERIODCODE_8_t, 16, 1, noOfThreads));
        case CODER_LEVEL_NORMAL:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 3 << 29, 8,
                    LZMA_DATAPERIODCODE_8_t, LZMA_DATAPERIODCODE_8_t, 32, 1, noOfThreads));
        case CODER_LEVEL_MAX:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 3 << 29, 8,
                    LZMA_DATAPERIODCODE_8_t, LZMA_DATAPERIODCODE_8_t, 128, 1, noOfThreads));
        default:
            fprintf(stderr, "Unsupported coder level: %d.\n", (int) coder_level);
            exit(EXIT_FAILURE);
    }
}

unique_ptr<CoderProps> getRelativeOffsetDeltasOfPairsValueCoderProps(uint8_t coder_level) {
    int noOfThreads = PgHelpers::numberOfThreads > 1 ? 2 : 1;
    switch(coder_level) {
        case CODER_LEVEL_FAST:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 1 << 24, 8,
                    LZMA_DATAPERIODCODE_8_t, LZMA_DATAPERIODCODE_8_t, 5, 1, noOfThreads));
        case CODER_LEVEL_NORMAL:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 3 << 29, 8,
                    LZMA_DATAPERIODCODE_8_t, LZMA_DATAPERIODCODE_8_t, 16, 1, noOfThreads));
        case CODER_LEVEL_MAX:
            return unique_ptr<CoderProps>(new LzmaCoderProps(9, 3 << 29, 8,
                    LZMA_DATAPERIODCODE_8_t, LZMA_DATAPERIODCODE_8_t, 16, 1, noOfThreads));
        default:
            fprintf(stderr, "Unsupported coder level: %d.\n", (int) coder_level);
            exit(EXIT_FAILURE);
    }
}

unique_ptr<CoderProps> getDefaultRangeCoderProps(uint16_t max_read_len, uint8_t mismatches_count) {
    return unique_ptr<CoderProps>(getRangeCoderProps(max_read_len, mismatches_count));
}
