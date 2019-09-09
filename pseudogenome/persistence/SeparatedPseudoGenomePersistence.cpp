#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"
#include "../SeparatedPseudoGenomeBase.h"

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, const string &pseudoGenomePrefix,
            IndexesMapping* orgIndexesMapping, bool revComplPairFile) {
        clock_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(pseudoGenomePrefix, !revComplPairFile, true);
        builder.writePseudoGenome(pgb, orgIndexesMapping, revComplPairFile);
        builder.build();
        cout << "Writing (" << pseudoGenomePrefix << ") pseudo genome files in " << clock_millis() << " msec." << endl << endl;
    }

    void SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(SeparatedPseudoGenome *sPg,
            const string &pseudoGenomePrefix, bool skipPgSequence) {
        clock_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(sPg->getReadsList()->isRevCompEnabled(),
                sPg->getReadsList()->areMismatchesEnabled());
        builder.feedSeparatedPseudoGenome(sPg, skipPgSequence);
        builder.build(pseudoGenomePrefix);
        cout << "Writing (" << pseudoGenomePrefix << ") pseudo genome files in " << clock_millis() << " msec." << endl << endl;
    }

    void SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(SeparatedPseudoGenome *sPg,
                                                                      ostream* pgrcOut, uint8_t coder_level,
                                                                      bool ignoreOffDest, bool skipPgSequence) {
        clock_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(!sPg->getReadsList()->isRevCompEnabled(),
                                                   !sPg->getReadsList()->areMismatchesEnabled());
        builder.feedSeparatedPseudoGenome(sPg, true);
        builder.compressedBuild(*pgrcOut, coder_level, ignoreOffDest);
        *logout << "Compressed Pg reads list in " << clock_millis() << " msec." << endl << endl;
        if (!skipPgSequence) {
            clock_checkpoint();
            writeCompressed(*pgrcOut, sPg->getPgSequence().data(), sPg->getPgSequence().size(), LZMA_CODER, coder_level,
                            PGRC_DATAPERIODCODE_8_t, COMPRESSION_ESTIMATION_BASIC_DNA);
            *logout << "Compressed Pg sequence in " << clock_millis() << " msec." << endl << endl;
        }
    }

    SeparatedPseudoGenome* SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(const string &pgPrefix,
            bool skipReadsList){
        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string pgSequence = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(pgPrefix);
        ExtendedReadsListWithConstantAccessOption* caeRl = skipReadsList?0:
                ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgPrefix, pgh->getPseudoGenomeLength());
        SeparatedPseudoGenome* spg = new SeparatedPseudoGenome(std::move(pgSequence), caeRl, prop);
        delete(pgh);
        delete(prop);
        return spg;
    }


    std::ifstream SeparatedPseudoGenomePersistence::getPseudoGenomeSrc(const string &pseudoGenomePrefix) {
        const string pgFile = pseudoGenomePrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX;
        std::ifstream pgSrc(pgFile, std::ios::in | std::ios::binary);
        if (pgSrc.fail()) {
            fprintf(stderr, "cannot open pseudogenome file %s\n", pgFile.c_str());
            exit(EXIT_FAILURE);
        }
        return pgSrc;
    }

    string SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(const string &pseudoGenomePrefix) {
        std::ifstream pgSrc = getPseudoGenomeSrc(pseudoGenomePrefix);
        pgSrc.seekg(0, std::ios::end);
        size_t size = pgSrc.tellg();
        std::string buffer(size, ' ');
        pgSrc.seekg(0);
        PgSAHelpers::readArray(pgSrc, &buffer[0], size);
        return buffer;
    }


    void SeparatedPseudoGenomePersistence::writePseudoGenomeSequence(string &pgSequence, string pgPrefix) {
        ofstream pgDest =
                getPseudoGenomeElementDest(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX, true);
        PgSAHelpers::writeArray(pgDest, (void*) pgSequence.data(), pgSequence.length());
        pgDest.close();
        acceptTemporaryPseudoGenomeElement(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX, true);
    }

    std::ifstream SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(const string &pseudoGenomePrefix,
                                                                              const string &fileSuffix) {
        const string pgElFile = pseudoGenomePrefix + fileSuffix;
        std::ifstream pgElSrc(pgElFile, std::ios::in | std::ios::binary);
        if (pgElSrc.fail())
            fprintf(stderr, "warning: pseudogenome element file %s does not exist (or cannot be opened for reading)\n",
                    pgElFile.c_str());

        return pgElSrc;
    }

    std::ofstream SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(const string &pseudoGenomePrefix,
                                                                               const string &fileSuffix,
                                                                               bool temporary) {
        string pgElFile = pseudoGenomePrefix + fileSuffix;
        if (temporary)
            pgElFile = pgElFile + TEMPORARY_FILE_SUFFIX;

        std::ofstream destPgEl(pgElFile, std::ios::out | std::ios::binary | std::ios::trunc);
        return destPgEl;
    }


    bool SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElement(const string &pseudoGenomePrefix,
                                                                              const string &fileSuffix,
                                                                              bool alwaysRemoveExisting) {
        string pgElFile = pseudoGenomePrefix + fileSuffix;
        string pgElTempFile = pgElFile + TEMPORARY_FILE_SUFFIX;
        if (std::ifstream(pgElFile) && (std::ifstream(pgElTempFile) || alwaysRemoveExisting))
            remove(pgElFile.c_str());
        if (std::ifstream(pgElTempFile))
            return rename(pgElTempFile.c_str(), pgElFile.c_str()) == 0;
        return false;
    }

    void SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(const string &pseudoGenomePrefix,
            bool clearReadsListDescriptionIfOverride) {
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX, false);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX, false);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_POSITIONS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_REVERSECOMPL_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_OFFSETS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX, clearReadsListDescriptionIfOverride);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX, false);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX, false);
    }

    const string SeparatedPseudoGenomePersistence::TEMPORARY_FILE_SUFFIX = ".temp";

    bool SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = false;
    bool SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = true;

    void SeparatedPseudoGenomePersistence::appendIndexesFromPg(string pgFilePrefix, vector<uint_reads_cnt_std> &idxs) {
        bool plainTextReadMode;
        PseudoGenomeHeader* pgh;
        ReadsSetProperties* rsProp;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgFilePrefix, pgh, rsProp, plainTextReadMode);
        ifstream orgIdxsSrc = getPseudoGenomeElementSrc(pgFilePrefix, SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        const uint_reads_cnt_max readsCount = pgh->getReadsCount();
        idxs.reserve(idxs.size() + readsCount);
        for(uint_reads_cnt_std i = 0; i < readsCount; i++) {
            uint_reads_cnt_std idx;
            readValue<uint_reads_cnt_std>(orgIdxsSrc, idx, plainTextReadMode);
            idxs.push_back(idx);
        }
        delete(pgh);
        delete(rsProp);
    }

    void SeparatedPseudoGenomePersistence::writePairMapping(string &pgFilePrefix,
                                                            vector<uint_reads_cnt_std> orgIdxs) {
        ofstream pair1OffsetsDest = getPseudoGenomeElementDest(pgFilePrefix, SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX, true);
        ofstream pair1SrcFlagDest = getPseudoGenomeElementDest(pgFilePrefix, SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX, true);
        ofstream pair1IndexesDest = getPseudoGenomeElementDest(pgFilePrefix, SeparatedPseudoGenomeBase::READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX, true);
        writeReadMode(pair1OffsetsDest, PgSAHelpers::plainTextWriteMode);
        writeReadMode(pair1IndexesDest, PgSAHelpers::plainTextWriteMode);
        uint_reads_cnt_std readsCount = orgIdxs.size();
        vector<uint_reads_cnt_std> rev(readsCount);
        for(uint_reads_cnt_std i = 0; i < readsCount; i++)
            rev[orgIdxs[i]] = i;
        vector<bool> isReadDone(readsCount, false);
        for(uint32_t i = 0; i < readsCount; i++) {
            if (isReadDone[i])
                continue;
            uint_reads_cnt_std idx = orgIdxs[i];
            uint_reads_cnt_std pairIdx = idx % 2?(idx-1):(idx+1);
            uint_reads_cnt_std pairI = rev[pairIdx];
            isReadDone[pairI] = true;
            writeValue<uint_reads_cnt_std>(pair1OffsetsDest, pairI > i?pairI - i: readsCount - (i - pairI));
            writeValue<uint8_t>(pair1SrcFlagDest, idx % 2);
            writeValue<uint_reads_cnt_std>(pair1IndexesDest, idx);
        }
        pair1IndexesDest.close();
        pair1OffsetsDest.close();
        acceptTemporaryPseudoGenomeElements(pgFilePrefix, false);
    }

    void SeparatedPseudoGenomePersistence::dumpPgPairs(vector<string> pgFilePrefixes) {
        clock_checkpoint();
        vector<uint_reads_cnt_std> orgIdxs;
        for(string pgFilePrefix: pgFilePrefixes)
            SeparatedPseudoGenomePersistence::appendIndexesFromPg(pgFilePrefix, orgIdxs);

        SeparatedPseudoGenomePersistence::writePairMapping(pgFilePrefixes[0], orgIdxs);
        *logout << "... dumping pairs completed in " << clock_millis() << " msec. " << endl;
    }

    void SeparatedPseudoGenomePersistence::compressReadsOrder(ostream &pgrcOut,
            const vector<uint_reads_cnt_std>& orgIdxs, uint8_t coder_level,
            bool completeOrderInfo, bool ignorePairOrderInformation, bool singleFileMode) {
        clock_checkpoint();
        uint_reads_cnt_std readsCount = orgIdxs.size();
        int lzma_reads_dataperiod_param = readsCount <= UINT32_MAX ? PGRC_DATAPERIODCODE_32_t : PGRC_DATAPERIODCODE_64_t;
        vector<uint_reads_cnt_std> rev(readsCount);
        for (uint_reads_cnt_std i = 0; i < readsCount; i++)
            rev[orgIdxs[i]] = i;
        if (completeOrderInfo && singleFileMode) {
            *logout << "Reverse index of original indexes... ";
            writeCompressed(pgrcOut, (char *) rev.data(), rev.size() * sizeof(uint_reads_cnt_std), LZMA_CODER,
                    coder_level, lzma_reads_dataperiod_param);
        } else {
            // absolute pair base index of original pair
            vector<uint_reads_cnt_std> revPairBaseOrgIdx;
            if (completeOrderInfo)
                revPairBaseOrgIdx.resize(readsCount / 2);
            // flag indicating a processed pair base file (0 - Second, 1 - First)
            vector<uint8_t> offsetPairBaseFileFlag;
            vector<uint8_t> nonOffsetPairBaseFileFlag;
            if (!ignorePairOrderInformation) {
                offsetPairBaseFileFlag.reserve(readsCount / 2);
                nonOffsetPairBaseFileFlag.reserve(readsCount / 4);
            }

            // coding reads list index relative offset of paired read
            vector<uint8_t> offsetInUint8Flag;
            offsetInUint8Flag.reserve(readsCount / 2);
            vector<uint8_t> offsetInUint8Value;
            offsetInUint8Value.reserve(readsCount / 2); // estimated
            vector<uint8_t> deltaInInt8Flag;
            deltaInInt8Flag.reserve(readsCount / 4); // estimated
            vector<int8_t> deltaInInt8Value;
            deltaInInt8Value.reserve(readsCount / 16); // estimated
            vector<uint_reads_cnt_std> fullOffset;
            deltaInInt8Value.reserve(readsCount / 8); // estimated

            vector<bool> isReadDone(readsCount, false);
            int64_t refPrev = 0;
            int64_t prev = 0;
            bool match = false;
            for (uint32_t i1 = 0; i1 < readsCount; i1++) {
                if (isReadDone[i1])
                    continue;
                uint_reads_cnt_std orgIdx = orgIdxs[i1];
                uint_reads_cnt_std pairOrgIdx = orgIdx % 2 ? (orgIdx - 1) : (orgIdx + 1);
                uint_reads_cnt_std i2 = rev[pairOrgIdx]; // i2 > i1
                isReadDone[i2] = true;
                if (completeOrderInfo)
                    revPairBaseOrgIdx[orgIdx / 2] = offsetInUint8Flag.size() * 2 + orgIdx % 2;
                int64_t pairRelativeOffset = i2 - i1;
                offsetInUint8Flag.push_back((uint8_t) (pairRelativeOffset <= UINT8_MAX));
                if (pairRelativeOffset <= UINT8_MAX) {
                    offsetInUint8Value.push_back((uint8_t) pairRelativeOffset);
                    if (!completeOrderInfo && !ignorePairOrderInformation)
                        offsetPairBaseFileFlag.push_back(orgIdx % 2);
                    continue;
                }
                if (!completeOrderInfo && !ignorePairOrderInformation)
                    nonOffsetPairBaseFileFlag.push_back(orgIdx % 2);
                const int64_t delta = pairRelativeOffset - refPrev;
                const bool isDeltaInInt8 = delta <= INT8_MAX && delta >= INT8_MIN;
                deltaInInt8Flag.push_back((uint8_t) isDeltaInInt8);
                if (isDeltaInInt8) {
                    match = true;
                    deltaInInt8Value.push_back((int8_t) delta);
                    refPrev = pairRelativeOffset;
                } else {
                    if (!match || refPrev != prev)
                        refPrev = pairRelativeOffset;
                    fullOffset.push_back(pairRelativeOffset);
                    match = false;
                }
                prev = pairRelativeOffset;
            }

            *logout << "Uint8 reads list relative offsets of pair reads (flag)... ";
            writeCompressed(pgrcOut, (char *) offsetInUint8Flag.data(), offsetInUint8Flag.size() * sizeof(uint8_t),
                            PPMD7_CODER, coder_level, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
            *logout << "Uint8 reads list relative offsets of pair reads (value)... ";
            writeCompressed(pgrcOut, (char *) offsetInUint8Value.data(), offsetInUint8Value.size() * sizeof(uint8_t),
                            PPMD7_CODER, coder_level, 2);
            *logout << "Relative offsets deltas of pair reads (flag)... ";
            writeCompressed(pgrcOut, (char *) deltaInInt8Flag.data(), deltaInInt8Flag.size() * sizeof(uint8_t),
                            PPMD7_CODER, coder_level, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
            *logout << "Relative offsets deltas of pair reads (value)... ";
            writeCompressed(pgrcOut, (char*) deltaInInt8Value.data(), deltaInInt8Value.size() * sizeof(uint8_t),
                            LZMA_CODER, coder_level, PGRC_DATAPERIODCODE_8_t);
//            writeCompressed(pgrcOut, (char *) deltaInInt8Value.data(), deltaInInt8Value.size() * sizeof(int8_t), PPMD7_CODER, coder_level, 2);
            *logout << "Full reads list relative offsets of pair reads ... ";
            double estimated_reads_ratio = simpleUintCompressionEstimate(readsCount, readsCount <= UINT32_MAX?UINT32_MAX:UINT64_MAX);
            writeCompressed(pgrcOut, (char *) fullOffset.data(), fullOffset.size() * sizeof(uint_reads_cnt_std),
                            LZMA_CODER, coder_level, lzma_reads_dataperiod_param, estimated_reads_ratio);
            if (completeOrderInfo) {
                *logout << "Original indexes of pair bases... ";
                writeCompressed(pgrcOut, (char *) revPairBaseOrgIdx.data(), revPairBaseOrgIdx.size() * sizeof(uint_reads_cnt_std),
                                LZMA_CODER, coder_level, lzma_reads_dataperiod_param, estimated_reads_ratio);
            } else if (!ignorePairOrderInformation) {
                *logout << "File flags of pair bases (for offsets)... ";
                writeCompressed(pgrcOut, (char *) offsetPairBaseFileFlag.data(), offsetPairBaseFileFlag.size() * sizeof(uint8_t),
                                PPMD7_CODER, coder_level, 2, COMPRESSION_ESTIMATION_UINT8_BITMAP);
                *logout << "File flags of pair bases (for non-offsets)... ";
                writeCompressed(pgrcOut, (char *) nonOffsetPairBaseFileFlag.data(), nonOffsetPairBaseFileFlag.size() * sizeof(uint8_t),
                                PPMD7_CODER, coder_level, 2, COMPRESSION_ESTIMATION_UINT8_BITMAP);
            }
        }
        *logout << "... compressing order information completed in " << clock_millis() << " msec. " << endl;
        *logout << endl;
    }

    void SeparatedPseudoGenomePersistence::decompressReadsOrder(istream &pgrcIn,
                                                                vector<uint_reads_cnt_std> &rlIdxOrder,
                                                                bool completeOrderInfo, bool ignorePairOrderInformation,
                                                                bool singleFileMode) {
        if (singleFileMode) {
            if (!completeOrderInfo)
                return;
            readCompressed<uint_reads_cnt_std>(pgrcIn, rlIdxOrder);
        } else {
            vector<uint8_t> offsetInUint8Flag;
            vector<uint8_t> offsetInUint8Value;
            vector<uint8_t> deltaInInt8Flag;
            vector<int8_t> deltaInInt8Value;
            vector<uint_reads_cnt_std> fullOffset;
            readCompressed<uint8_t>(pgrcIn, offsetInUint8Flag);
            readCompressed<uint8_t>(pgrcIn, offsetInUint8Value);
            readCompressed<uint8_t>(pgrcIn, deltaInInt8Flag);
            readCompressed<int8_t>(pgrcIn, deltaInInt8Value);
            readCompressed<uint_reads_cnt_std>(pgrcIn, fullOffset);

            uint_reads_cnt_std readsCount = offsetInUint8Flag.size() * 2;
            rlIdxOrder.resize(readsCount);
            vector<bool> isReadDone(readsCount, false);
            int64_t pairOffset = 0;
            int64_t pairCounter = -1;
            int64_t offIdx = -1;
            int64_t delFlagIdx = -1;
            int64_t delIdx = -1;
            int64_t fulIdx = -1;
            int64_t refPrev = 0;
            int64_t prev = 0;
            bool match = false;

            for(uint_reads_cnt_std i = 0; i < readsCount; i++) {
                if (isReadDone[i])
                    continue;
                if (offsetInUint8Flag[++pairCounter])
                    pairOffset = offsetInUint8Value[++offIdx];
                else if (deltaInInt8Flag[++delFlagIdx]) {
                    pairOffset = refPrev + deltaInInt8Value[++delIdx];
                    refPrev = pairOffset;
                    match = true;
                    prev = pairOffset;
                } else {
                    pairOffset = fullOffset[++fulIdx];
                    if (!match || refPrev != prev)
                        refPrev = pairOffset;
                    match = false;
                    prev = pairOffset;
                }
                rlIdxOrder[pairCounter * 2] = i;
                rlIdxOrder[pairCounter * 2 + 1] = i + pairOffset;
                isReadDone[i + pairOffset] = true;
            }
            if (completeOrderInfo) {
                vector<uint_reads_cnt_std> revPairBaseOrgIdx;
                revPairBaseOrgIdx.reserve(readsCount);
                readCompressed<uint_reads_cnt_std>(pgrcIn, revPairBaseOrgIdx);
                revPairBaseOrgIdx.resize(readsCount);
                vector<uint_reads_cnt_std> peRlIdxOrder = std::move(rlIdxOrder);
                for(uint_reads_cnt_std p = readsCount / 2; p-- > 0;) {
                    uint_reads_cnt_std rlIdx = revPairBaseOrgIdx[p];
                    revPairBaseOrgIdx[p * 2] = peRlIdxOrder[rlIdx];
                    revPairBaseOrgIdx[p * 2 + 1] = peRlIdxOrder[rlIdx % 2?rlIdx - 1:rlIdx + 1];
                }
                rlIdxOrder = std::move(revPairBaseOrgIdx);
            } else if (!ignorePairOrderInformation) {
                vector<uint8_t> offsetPairBaseFileFlag;
                vector<uint8_t> nonOffsetPairBaseFileFlag;
                readCompressed<uint8_t>(pgrcIn, offsetPairBaseFileFlag);
                readCompressed<uint8_t>(pgrcIn, nonOffsetPairBaseFileFlag);
                int64_t offIdx = -1;
                int64_t nonOffIdx = -1;
                for(uint_reads_cnt_std p = 0; p < readsCount / 2; p++) {
                    bool swapPair = offsetInUint8Flag[p]?offsetPairBaseFileFlag[++offIdx]:nonOffsetPairBaseFileFlag[++nonOffIdx];
                    if (swapPair) {
                        uint_reads_cnt_std tmpIdx = rlIdxOrder[p * 2];
                        rlIdxOrder[p * 2] = rlIdxOrder[p * 2 + 1];
                        rlIdxOrder[p * 2 + 1] = tmpIdx;
                    }
                }
            }
        }
    }

    template <typename uint_pg_len>
    void SeparatedPseudoGenomePersistence::compressReadsPgPositions(ostream &pgrcOut,
            vector<uint_pg_len_max> orgIdx2PgPos, uint_pg_len_max joinedPgLength, uint8_t coder_level,
            bool singleFileMode, bool deltaPairEncodingEnabled) {
        clock_checkpoint();
        uint_reads_cnt_std readsTotalCount = orgIdx2PgPos.size();
        int lzma_pos_dataperiod_param = sizeof(uint_pg_len) == 4 ? PGRC_DATAPERIODCODE_32_t : PGRC_DATAPERIODCODE_64_t;
        double estimated_pos_ratio = simpleUintCompressionEstimate(joinedPgLength, sizeof(uint_pg_len) == 4?UINT32_MAX:UINT64_MAX);
        if (singleFileMode) {
            uint_pg_len_max* const maxPgPosPtr = orgIdx2PgPos.data();
            if (sizeof(uint_pg_len) < sizeof(uint_pg_len_max)) {
                uint_pg_len* PgPosPtr = (uint_pg_len*) maxPgPosPtr;
                for (uint_reads_cnt_std i = 0; i < readsTotalCount; i++)
                    *(PgPosPtr++) = (uint_pg_len) orgIdx2PgPos[i];
            }
            writeCompressed(pgrcOut, (char*) maxPgPosPtr, readsTotalCount * sizeof(uint_pg_len),
                            LZMA_CODER, coder_level, lzma_pos_dataperiod_param, estimated_pos_ratio);
        } else {
            vector<uint_pg_len> basePairPos;
            // pair relative offset info
            vector<uint8_t> offsetInUint16Flag;
            vector<uint8_t> offsetIsBaseFirstFlag;
            vector<uint16_t> offsetInUint16Value;
            vector<uint8_t> deltaInInt16Flag;
            vector<uint8_t> deltaIsBaseFirstFlag;
            vector<int16_t> deltaInInt16Value;
            vector<uint_pg_len> notBasePairPos;
            const uint_reads_cnt_std pairsCount = readsTotalCount / 2;
            basePairPos.reserve(pairsCount);
            offsetInUint16Flag.reserve(pairsCount);
            offsetIsBaseFirstFlag.reserve(pairsCount);
            offsetInUint16Value.reserve(pairsCount);
            if (deltaPairEncodingEnabled) {
                deltaInInt16Flag.reserve(pairsCount / 2);
                deltaIsBaseFirstFlag.reserve(pairsCount / 8);
                deltaInInt16Value.reserve(pairsCount / 8);
            }
            notBasePairPos.reserve(pairsCount / 4);

            vector<uint_reads_cnt_std> bppRank;
            bppRank.reserve(pairsCount);
            for (uint_reads_cnt_std i = 0; i < readsTotalCount; i += 2) {
                basePairPos.push_back(orgIdx2PgPos[i]);
                bppRank.push_back(i >> 1);
            }
            std::stable_sort(bppRank.begin(), bppRank.end(),
                    [&](const uint_reads_cnt_std &idx1, const uint_reads_cnt_std &idx2) -> bool
                        { return basePairPos[idx1] < basePairPos[idx2]; });
            *logout << "... reordering bases checkpoint: " << clock_millis() << " msec. " << endl;
            int64_t refPrev = 0;
            int64_t prev = 0;
            bool match = false;
            for (uint_reads_cnt_std p = 0; p < pairsCount; p++) {
                uint_reads_cnt_std i = bppRank[p] * 2;

                bool isBaseBefore = orgIdx2PgPos[i] < orgIdx2PgPos[i + 1];
                uint_pg_len relativeAbsOffset = isBaseBefore?(orgIdx2PgPos[i + 1] - orgIdx2PgPos[i]):
                                    orgIdx2PgPos[i] - orgIdx2PgPos[i + 1];
                const bool isOffsetInUint16 = relativeAbsOffset <= UINT16_MAX;
                offsetInUint16Flag.push_back(isOffsetInUint16 ? 1 : 0);
                if (isOffsetInUint16) {
                    offsetIsBaseFirstFlag.push_back(isBaseBefore?1:0);
                    offsetInUint16Value.push_back((uint16_t) relativeAbsOffset);
                    continue;
                }
                if (deltaPairEncodingEnabled) {
                    const int64_t delta = relativeAbsOffset - refPrev;
                    const bool isDeltaInInt16 = delta <= INT16_MAX && delta >= INT16_MIN;
                    deltaInInt16Flag.push_back((uint8_t) isDeltaInInt16);
                    if (isDeltaInInt16) {
                        match = true;
                        deltaIsBaseFirstFlag.push_back(isBaseBefore ? 1 : 0);
                        deltaInInt16Value.push_back((int16_t) delta);
                        refPrev = relativeAbsOffset;
                    } else {
                        if (!match || refPrev != prev)
                            refPrev = relativeAbsOffset;
                        notBasePairPos.push_back((uint_pg_len) orgIdx2PgPos[i + 1]);
                        match = false;
                    }
                    prev = relativeAbsOffset;
                } else {
                    notBasePairPos.push_back((uint_pg_len) orgIdx2PgPos[i + 1]);
                }
            }
            pgrcOut.put(deltaPairEncodingEnabled);
            *logout << "Base pair position... ";
            writeCompressed(pgrcOut, (char *) basePairPos.data(), basePairPos.size() * sizeof(uint_pg_len),
                            LZMA_CODER, coder_level, lzma_pos_dataperiod_param, estimated_pos_ratio);
            *logout << "Uint16 relative offset of pair positions (flag)... ";
            writeCompressed(pgrcOut, (char *) offsetInUint16Flag.data(), offsetInUint16Flag.size() * sizeof(uint8_t),
                            PPMD7_CODER, coder_level, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
            *logout << "Is uint16 relative offset of pair positions positive (flag)... ";
            writeCompressed(pgrcOut, (char *) offsetIsBaseFirstFlag.data(), offsetIsBaseFirstFlag.size() * sizeof(uint8_t),
                            PPMD7_CODER, coder_level, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
            *logout << "Uint16 relative offset of pair positions (value)... ";
            writeCompressed(pgrcOut, (char*) offsetInUint16Value.data(), offsetInUint16Value.size() * sizeof(uint16_t),
                            PPMD7_CODER, coder_level, 3);
            if (deltaPairEncodingEnabled) {
                *logout << "Relative offset deltas of pair positions (flag)... ";
                writeCompressed(pgrcOut, (char *) deltaInInt16Flag.data(), deltaInInt16Flag.size() * sizeof(uint8_t),
                                PPMD7_CODER, coder_level, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
                *logout << "Is relative offset (for deltas stream) of pair positions positive (flag)... ";
                writeCompressed(pgrcOut, (char *) deltaIsBaseFirstFlag.data(),
                                deltaIsBaseFirstFlag.size() * sizeof(uint8_t),
                                PPMD7_CODER, coder_level, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
                *logout << "Relative offset deltas of pair positions (value)... ";
                writeCompressed(pgrcOut, (char *) deltaInInt16Value.data(), deltaInInt16Value.size() * sizeof(int16_t),
                                PPMD7_CODER, coder_level, 3);
            }
            *logout << "Not-base pair position... ";
            writeCompressed(pgrcOut, (char *) notBasePairPos.data(), notBasePairPos.size() * sizeof(uint_pg_len),
                            LZMA_CODER, coder_level, lzma_pos_dataperiod_param, estimated_pos_ratio);
        }
        *logout << "... compressing reads positions completed in " << clock_millis() << " msec. " << endl;
        *logout << endl;
    }
    template void SeparatedPseudoGenomePersistence::compressReadsPgPositions<uint_pg_len_std>(ostream &pgrcOut,
            vector<uint_pg_len_max> orgIdx2PgPos, uint_pg_len_max joinedPgLength, uint8_t coder_level,
            bool singleFileMode, bool deltaPairEncodingEnabled);
    template void SeparatedPseudoGenomePersistence::compressReadsPgPositions<uint_pg_len_max>(ostream &pgrcOut,
            vector<uint_pg_len_max> orgIdx2PgPos, uint_pg_len_max joinedPgLength, uint8_t coder_level,
            bool singleFileMode, bool deltaPairEncodingEnabled);

    template <typename uint_pg_len>
    void SeparatedPseudoGenomePersistence::decompressReadsPgPositions(istream &pgrcIn, vector<uint_pg_len> &pgPos,
                                                                          uint_reads_cnt_std readsTotalCount, bool singleFileMode) {
        if (singleFileMode)
            readCompressed(pgrcIn, pgPos);
        else {
            bool deltaPairEncodingEnabled = (bool) pgrcIn.get();
            vector<uint8_t> offsetInUint16Flag;
            vector<uint8_t> offsetIsBaseFirstFlag;
            vector<uint16_t> offsetInUint16Value;
            vector<uint8_t> deltaInInt16Flag;
            vector<uint8_t> deltaIsBaseFirstFlag;
            vector<int16_t> deltaInInt16Value;
            vector<uint_pg_len> notBasePairPos;
            pgPos.reserve(readsTotalCount);
            readCompressed(pgrcIn, pgPos);
            readCompressed(pgrcIn, offsetInUint16Flag);
            readCompressed(pgrcIn, offsetIsBaseFirstFlag);
            readCompressed(pgrcIn, offsetInUint16Value);
            if (deltaPairEncodingEnabled) {
                readCompressed(pgrcIn, deltaInInt16Flag);
                readCompressed(pgrcIn, deltaIsBaseFirstFlag);
                readCompressed(pgrcIn, deltaInInt16Value);
            }
            readCompressed(pgrcIn, notBasePairPos);
            const uint_reads_cnt_std pairsCount = readsTotalCount / 2;
            vector<uint_reads_cnt_std> bppRank;
            bppRank.reserve(pairsCount);
            for (uint_reads_cnt_std p = 0; p < pairsCount; p++)
                bppRank.push_back(p);
            std::stable_sort(bppRank.begin(), bppRank.end(),
                             [&](const uint_reads_cnt_std &idx1, const uint_reads_cnt_std &idx2) -> bool
                             { return pgPos[idx1] < pgPos[idx2]; });

            pgPos.resize(readsTotalCount);
            int64_t nbpPos = 0;
            int64_t offIdx = -1;
            int64_t delFlagIdx = -1;
            int64_t delIdx = -1;
            int64_t nbpPosIdx = -1;
            int64_t refPrev = 0;
            int64_t prev = 0;
            bool match = false;
            for (uint_reads_cnt_std i = 0; i < pairsCount; i++) {
                uint_reads_cnt_std p = bppRank[i];
                if (offsetInUint16Flag[i] == 1) {
                    int64_t delta = offsetInUint16Value[++offIdx];
                    if (offsetIsBaseFirstFlag[offIdx] == 0)
                        delta = -delta;
                    nbpPos = pgPos[p] + delta;
                } else if (deltaInInt16Flag[++delFlagIdx]){
                    int64_t delta = refPrev + deltaInInt16Value[++delIdx];
                    refPrev = delta;
                    prev = delta;
                    if (deltaIsBaseFirstFlag[delIdx] == 0)
                        delta = -delta;
                    nbpPos = pgPos[p] + delta;
                    match = true;
                } else {
                    nbpPos = notBasePairPos[++nbpPosIdx];
                    int64_t delta = nbpPos - pgPos[p];
                    if (delta < 0)
                        delta = -delta;
                    if (!match || refPrev != prev)
                        refPrev = delta;
                    match = false;
                    prev = delta;
                }
                pgPos[pairsCount + p] = nbpPos;
            }
        }
    }
    template void SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_std>(istream &pgrcIn, vector<uint_pg_len_std> &pgPos, uint_reads_cnt_std readsTotalCount, bool singleFileMode);
    template void SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_max>(istream &pgrcIn, vector<uint_pg_len_max> &pgPos, uint_reads_cnt_std readsTotalCount, bool singleFileMode);

    SeparatedPseudoGenomeOutputBuilder::SeparatedPseudoGenomeOutputBuilder(const string pseudoGenomePrefix,
            bool disableRevComp, bool disableMismatches) : pseudoGenomePrefix(pseudoGenomePrefix),
            disableRevComp(disableRevComp), disableMismatches(disableMismatches) {
        initReadsListDests();
    }

    SeparatedPseudoGenomeOutputBuilder::SeparatedPseudoGenomeOutputBuilder(bool disableRevComp, bool disableMismatches)
            : SeparatedPseudoGenomeOutputBuilder("", disableRevComp, disableMismatches) {
        initReadsListDests();
    }

    void SeparatedPseudoGenomeOutputBuilder::initDest(ostream *&dest, const string &fileSuffix) {
        if (dest == 0) {
            if (onTheFlyMode())
                dest = new ofstream(
                        SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(pseudoGenomePrefix, fileSuffix, true));
            else
                dest = new ostringstream();
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::initReadsListDests() {
        if (SeparatedPseudoGenomePersistence::enableReadPositionRepresentation)
            initDest(rlPosDest, SeparatedPseudoGenomeBase::READSLIST_POSITIONS_FILE_SUFFIX);
        else
            initDest(rlOffDest, SeparatedPseudoGenomeBase::READSLIST_OFFSETS_FILE_SUFFIX);
        initDest(rlOrgIdxDest, SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        if (!disableRevComp)
            initDest(rlRevCompDest, SeparatedPseudoGenomeBase::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        if (!disableMismatches) {
            initDest(rlMisCntDest, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
            initDest(rlMisSymDest, SeparatedPseudoGenomeBase::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
            if (SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation)
                initDest(rlMisRevOffDest,
                         SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
            else
                initDest(rlMisPosDest, SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::destToFile(ostream *dest, const string &fileName) {
        if (!onTheFlyMode() && dest)
            PgSAHelpers::writeStringToFile(fileName, ((ostringstream*) dest)->str());
    }

    void SeparatedPseudoGenomeOutputBuilder::freeDest(ostream* &dest) {
        if (dest) {
            if (onTheFlyMode())
                ((ofstream*) dest)->close();
            delete(dest);
            dest = 0;
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::freeDests() {
        freeDest(pgDest);
        freeDest(pgPropDest);
        freeDest(rlPosDest);
        freeDest(rlOrgIdxDest);
        freeDest(rlRevCompDest);
        freeDest(rlMisCntDest);
        freeDest(rlMisSymDest);
        freeDest(rlMisPosDest);

        freeDest(rlOffDest);
        freeDest(rlMisRevOffDest);
    }

    void SeparatedPseudoGenomeOutputBuilder::prebuildAssert(bool requireOnTheFlyMode) {
        if (onTheFlyMode() != requireOnTheFlyMode) {
            if (requireOnTheFlyMode)
                fprintf(stderr, "This build mode is supported only in on-the-fly mode.\n");
            else
                fprintf(stderr, "This build mode is not supported in on-the-fly mode.\n");
            exit(EXIT_FAILURE);
        }
        if (pgh == 0) {
            fprintf(stderr, "Pseudo genome header not initialized in separated Pg builder.\n");
            exit(EXIT_FAILURE);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::buildProps() {
        pgh->setReadsCount(readsCounter);
        rsProp->readsCount = readsCounter;
        rsProp->allReadsLength = rsProp->constantReadLength ? readsCounter * rsProp->maxReadLength : -1;
        if (pgPropDest) {
            delete (pgPropDest);
            pgPropDest = 0;
        }
        initDest(pgPropDest, SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX);
        pgh->write(*pgPropDest);
        rsProp->write(*pgPropDest);
    }

    void SeparatedPseudoGenomeOutputBuilder::build() {
        prebuildAssert(true);
        buildProps();
        writeReadMode(*pgPropDest, plainTextWriteMode);

        freeDests();
        SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(pseudoGenomePrefix, true);
    }

    void SeparatedPseudoGenomeOutputBuilder::build(const string &pgPrefix) {
        if (pgPrefix.empty())
            return;
        prebuildAssert(false);
        buildProps();
        writeReadMode(*pgPropDest, plainTextWriteMode);
        PgSAHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX,
                                      ((ostringstream*) pgPropDest)->str());

        destToFile(rlPosDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_POSITIONS_FILE_SUFFIX);
        destToFile(rlOffDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_OFFSETS_FILE_SUFFIX);
        destToFile(rlOrgIdxDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        destToFile(rlRevCompDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        destToFile(rlMisCntDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
        destToFile(rlMisSymDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
        destToFile(rlMisRevOffDest,
                   pgPrefix + SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
        destToFile(rlMisPosDest, pgPrefix + SeparatedPseudoGenomeBase::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);
    }

    void SeparatedPseudoGenomeOutputBuilder::compressDest(ostream* dest, ostream &pgrcOut, uint8_t coder_type,
                                                          uint8_t coder_level, int coder_param,
                                                          double estimated_compression,
                                                          SymbolsPackingFacility* symPacker) {
        if (onTheFlyMode() || !dest) {
            fprintf(stderr, "Error during compression: an input stream missing.\n");
            exit(EXIT_FAILURE);
        }
        string tmp = ((ostringstream*) dest)->str();
        if (symPacker) {
            PgSAHelpers::writeValue<uint64_t>(pgrcOut, tmp.length(), false);
            tmp = symPacker->packSequence(tmp.data(), tmp.length());
        }
        writeCompressed(pgrcOut, tmp, coder_type, coder_level, coder_param, estimated_compression);
    }

    void SeparatedPseudoGenomeOutputBuilder::compressRlMisRevOffDest(ostream &pgrcOut, uint8_t coder_level,
            bool transposeMode) {
        uint8_t mismatches_dests_count = coder_level == PGRC_CODER_LEVEL_FAST?1:(UINT8_MAX-1);
        if (mismatches_dests_count == 1) {
            PgSAHelpers::writeValue<uint8_t>(pgrcOut, 1);
            compressDest(rlMisRevOffDest, pgrcOut, PPMD7_CODER, coder_level, 3);
            return;
        }
        vector<uint8_t> misCnt2DestIdx; // NOT-TESTED => {0, 1, 2, 3, 4, 5, 6, 7, 7, 9, 9, 9 };
        misCnt2DestIdx.insert(misCnt2DestIdx.end(), UINT8_MAX, mismatches_dests_count);
        for(uint8_t m = 1; m < mismatches_dests_count; m++)
            misCnt2DestIdx[m] = m;

        ostringstream dests[UINT8_MAX];
        istringstream misRevOffSrc(((ostringstream*) rlMisRevOffDest)->str());
        istringstream misCntSrc(((ostringstream*) rlMisCntDest)->str());

        uint8_t misCnt = 0;
        uint16_t revOff = 0;
        for(uint_reads_cnt_max i = 0; i < readsCounter; i++) {
             PgSAHelpers::readValue<uint8_t>(misCntSrc, misCnt, false);
             for(uint8_t m = 0; m < misCnt; m++) {
                 PgSAHelpers::readReadLengthValue(misRevOffSrc, revOff, false);
                 PgSAHelpers::writeReadLengthValue(dests[misCnt2DestIdx[misCnt]], revOff);
             }
        }


        if (transposeMode) {
            for (uint8_t d = 1; d < mismatches_dests_count; d++) {
                if (misCnt2DestIdx[d] == misCnt2DestIdx[d - 1] || misCnt2DestIdx[d] == misCnt2DestIdx[d + 1])
                    continue;
                string matrix = dests[d].str();
                uint64_t readsCount = matrix.size() / d / (bytePerReadLengthMode ? 1 : 2);
                if (bytePerReadLengthMode)
                    dests[d].str(transpose<uint8_t>(matrix, readsCount, d));
                else
                    dests[d].str(transpose<uint16_t>(matrix, readsCount, d));
            }
        }
        while (mismatches_dests_count > 0 && dests[mismatches_dests_count].tellp() == 0)
            mismatches_dests_count--;
        PgSAHelpers::writeValue<uint8_t>(pgrcOut, mismatches_dests_count);
        for(uint8_t m = 1; m < mismatches_dests_count; m++)
            PgSAHelpers::writeValue<uint8_t>(pgrcOut, misCnt2DestIdx[m]);
        for(uint8_t m = 1; m <= mismatches_dests_count; m++) {
            *logout << (int) m << ": ";
            compressDest(&dests[m], pgrcOut, PPMD7_CODER, coder_level, 2);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::compressedBuild(ostream &pgrcOut, uint8_t coder_level, bool ignoreOffDest) {
        prebuildAssert(false);
        buildProps();
        writeReadMode(*pgPropDest, false);
        const string tmp = ((ostringstream*) pgPropDest)->str();
        pgrcOut.write(tmp.data(), tmp.length());

        if (!ignoreOffDest) {
            *logout << "Reads list offsets... ";
            compressDest(rlOffDest, pgrcOut, PPMD7_CODER, coder_level, 3);
        }
        if (!this->disableRevComp) {
            *logout << "Reverse complements info... ";
            compressDest(rlRevCompDest, pgrcOut, PPMD7_CODER, coder_level, 2, COMPRESSION_ESTIMATION_UINT8_BITMAP);
        }
        if (!this->disableMismatches) {
            *logout << "Mismatches counts... ";
            compressDest(rlMisCntDest, pgrcOut, PPMD7_CODER, coder_level, 2, COMPRESSION_ESTIMATION_MIS_CNT);
            *logout << "Mismatched symbols codes... ";
            compressDest(rlMisSymDest, pgrcOut, PPMD7_CODER, coder_level, 2, COMPRESSION_ESTIMATION_MIS_SYM);
            *logout << "Mismatches offsets (rev-coded)... " << endl;
            compressRlMisRevOffDest(pgrcOut, coder_level);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::updateOriginalIndexesIn(SeparatedPseudoGenome *sPg) {
        string tmp = ((ostringstream *) rlOrgIdxDest)->str();
        uint_reads_cnt_std *orgIdxPtr = (uint_reads_cnt_std *) tmp.data();
        sPg->getReadsList()->orgIdx.assign(orgIdxPtr, orgIdxPtr + readsCounter);
        sPg->getReadsList()->readsCount = readsCounter;
    }

    void SeparatedPseudoGenomeOutputBuilder::writeReadEntry(const DefaultReadsListEntry &rlEntry) {
        lastWrittenPos = rlEntry.pos;
        if (SeparatedPseudoGenomePersistence::enableReadPositionRepresentation)
            PgSAHelpers::writeValue<uint_pg_len_max>(*rlPosDest, rlEntry.pos);
        else
            PgSAHelpers::writeReadLengthValue(*rlOffDest, rlEntry.offset);
        PgSAHelpers::writeValue<uint_reads_cnt_std>(*rlOrgIdxDest, rlEntry.idx);
        if (!disableRevComp)
            PgSAHelpers::writeValue<uint8_t>(*rlRevCompDest, rlEntry.revComp?1:0);
        if (!disableMismatches) {
            PgSAHelpers::writeValue<uint8_t>(*rlMisCntDest, rlEntry.mismatchesCount);
            if (rlEntry.mismatchesCount) {
                for (uint8_t i = 0; i < rlEntry.mismatchesCount; i++)
                    PgSAHelpers::writeValue<uint8_t>(*rlMisSymDest, rlEntry.mismatchCode[i]);
                if (SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation) {
                    uint8_t currentPos = pgh->getMaxReadLength() - 1;
                    for (int16_t i = rlEntry.mismatchesCount - 1; i >= 0; i--) {
                        PgSAHelpers::writeReadLengthValue(*rlMisRevOffDest,
                                                                   currentPos - rlEntry.mismatchOffset[i]);
                        currentPos = rlEntry.mismatchOffset[i] - 1;
                    }
                } else {
                    for (uint8_t i = 0; i < rlEntry.mismatchesCount; i++)
                        PgSAHelpers::writeReadLengthValue(*rlMisPosDest, rlEntry.mismatchOffset[i]);
                }
            }
        }
        readsCounter++;
    }

    void SeparatedPseudoGenomeOutputBuilder::writeExtraReadEntry(const DefaultReadsListEntry &rlEntry) {
        writeReadEntry(rlEntry);
        if (rlIt != 0) {
            ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &itEntry = this->rlIt->peekReadEntry();
            itEntry.offset -= rlEntry.offset;
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::setReadsSourceIterator(DefaultReadsListIteratorInterface *rlIt) {
        rlIt->rewind();
        this->rlIt = rlIt;
    }

    uint_pg_len_max SeparatedPseudoGenomeOutputBuilder::writeReadsFromIterator(uint_pg_len_max stopPos) {
        if (iterationPaused) {
            if (rlIt->peekReadEntry().pos >= stopPos)
                return lastWrittenPos;
            writeReadEntry(rlIt->peekReadEntry());
            iterationPaused = false;
        }
        while (rlIt->moveNext()) {
            if (rlIt->peekReadEntry().pos >= stopPos) {
                iterationPaused = true;
                return lastWrittenPos;
            }
            writeReadEntry(rlIt->peekReadEntry());
        }
        return -1;
    }

    void SeparatedPseudoGenomeOutputBuilder::writePseudoGenome(PseudoGenomeBase *pgb, IndexesMapping* orgIndexesMapping, bool revComplPairFile) {

        initDest(pgDest, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX);
        string pg = pgb->getPseudoGenomeVirtual();
        PgSAHelpers::writeArray(*pgDest, (void*) pg.data(), pg.length());

        ReadsListIteratorExtendedWrapperBase* rlIt =
                TemplateUserGenerator::generateReadsListUser<ReadsListIteratorExtendedWrapper, ReadsListIteratorExtendedWrapperBase>(pgb);
        rlIt->applyIndexesMapping(orgIndexesMapping);
        if (revComplPairFile)
            rlIt->applyRevComplPairFileFlag();
        if (pgb->isReadLengthMin())
            bytePerReadLengthMode = true;
        setReadsSourceIterator(rlIt);
        writeReadsFromIterator();

        if (pgh == 0)
            pgh = new PseudoGenomeHeader(pgb);
        if (rsProp == 0)
            rsProp = new ReadsSetProperties(*(pgb->getReadsSetProperties()));
        if (pgh->getReadsCount() != readsCounter) {
            fprintf(stderr, "Incorrect reads count validation while building separated Pg (%u instead of %u).\n",
                    readsCounter, pgh->getReadsCount());
            exit(EXIT_FAILURE);
        }

        delete(rlIt);
    }

    void SeparatedPseudoGenomeOutputBuilder::copyPseudoGenomeProperties(const string &pseudoGenomePrefix) {
        ifstream pgPropSrc(pseudoGenomePrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX,
                           ios_base::in | ios_base::binary);
        if (pgPropSrc.fail()) {
            fprintf(stderr, "Cannot read pseudogenome properties (%s does not open).\n",
                    pseudoGenomePrefix.c_str());
            exit(EXIT_FAILURE);
        }
        if (pgh != 0)
            delete(pgh);
        if (rsProp != 0)
            delete(rsProp);
        this->pgh = new PseudoGenomeHeader(pgPropSrc);
        this->rsProp = new ReadsSetProperties(pgPropSrc);
        pgPropSrc.close();
    }

    void SeparatedPseudoGenomeOutputBuilder::copyPseudoGenomeProperties(SeparatedPseudoGenome* sPg) {
        this->pgh = new PseudoGenomeHeader(sPg);
        this->rsProp = new ReadsSetProperties(*(sPg->getReadsSetProperties()));
        if (sPg->isReadLengthMin())
            bytePerReadLengthMode = true;
    }

    void SeparatedPseudoGenomeOutputBuilder::appendPseudoGenome(const string &pg) {
        if (pgDest) {
            pgh->setPseudoGenomeLength(pgh->getPseudoGenomeLength() + pg.length());
            PgSAHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        } else {
            pgh->setPseudoGenomeLength(pg.length());
            initDest(pgDest, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX);
            PgSAHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::feedSeparatedPseudoGenome(SeparatedPseudoGenome *sPg, bool skipPgSequence) {
        if (!skipPgSequence) {
            initDest(pgDest, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX);
            const string &pg = sPg->getPgSequence();
            PgSAHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        }
        if (sPg->isReadLengthMin())
            bytePerReadLengthMode = true;
        setReadsSourceIterator(sPg->getReadsList());
        if (pgh == 0)
            pgh = new PseudoGenomeHeader(sPg);
        writeReadsFromIterator();

        if (rsProp == 0)
            rsProp = new ReadsSetProperties(*(sPg->getReadsSetProperties()));
        if (pgh->getReadsCount() != readsCounter) {
            fprintf(stderr, "Incorrect reads count validation while building separated Pg (%u instead of %u).\n",
                    readsCounter, pgh->getReadsCount());
            exit(EXIT_FAILURE);
        }
    }

    SeparatedPseudoGenomeOutputBuilder::~SeparatedPseudoGenomeOutputBuilder() {
        if (pgh)
            delete(pgh);
        if (rsProp)
            delete(rsProp);
        freeDests();
    }
}
