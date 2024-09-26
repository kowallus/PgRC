#ifndef PGTOOLS_PARAMSLIBRARY_H
#define PGTOOLS_PARAMSLIBRARY_H

#include "../utils/helper.h"
#include "CodersLib.h"
#include <memory>

unique_ptr<CoderProps> getDefaultCoderProps(uint8_t coder_type, uint8_t coder_level, int coder_param = -1);

unique_ptr<CoderProps> getDefaultFSECoderProps(uint8_t tableLog = 0, uint8_t maxSymbolValue = 0);

unique_ptr<CoderProps> getDefaultHufCoderProps(uint8_t tableLog = 0, uint8_t maxSymbolValue = 0);

unique_ptr<CoderProps> getCompoundCoderProps(CoderProps* firstCoderProps, CoderProps* secondCoderProps);

unique_ptr<CoderProps> getSelectorCoderProps(std::initializer_list<CoderProps*> coderPropsList, float probeFraction = 1,
                                             size_t minProbeSize = SelectorCoderProps::DEFAULT_MIN_PROBE_SIZE);

unique_ptr<CoderProps> getReadsPositionsCoderProps(uint8_t coder_level, uint8_t lzma_pos_dataperiod_param);

unique_ptr<CoderProps> getVarLenEncodedPgCoderProps(uint8_t coder_level);

unique_ptr<CoderProps> getRelativeOffsetDeltasOfPairsValueCoderProps(uint8_t coder_level);

unique_ptr<CoderProps> getDefaultRangeCoderProps(uint16_t no_of_symbols = 256, uint8_t bytes_period = 1);

#endif //PGTOOLS_PARAMSLIBRARY_H
