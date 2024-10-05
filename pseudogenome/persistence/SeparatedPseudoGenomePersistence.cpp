#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"
#include "../SeparatedPseudoGenomeBase.h"
#include "../../coders/PropsLibrary.h"
#include <memory>

#ifdef __APPLE__
#define PSTLD_HEADER_ONLY
#include "../../utils/pstld.h"
#define parallel_algorithm pstld
#else
#include <parallel/algorithm>
#define parallel_algorithm __gnu_parallel
#endif

#include <cassert>

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, const string &pseudoGenomePrefix,
            IndexesMapping* orgIndexesMapping, bool revComplPairFile) {
        time_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(pseudoGenomePrefix, !revComplPairFile, true);
        builder.writePseudoGenome(pgb, orgIndexesMapping, revComplPairFile);
        builder.build();
        cout << "Writing (" << pseudoGenomePrefix << ") pseudo genome files in " << time_millis() << " msec." << endl << endl;
    }

    void SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(SeparatedPseudoGenome *sPg,
            const string &pseudoGenomePrefix, bool skipPgSequence) {
        time_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(sPg->getReadsList()->isRevCompEnabled(),
                sPg->getReadsList()->areMismatchesEnabled());
        builder.feedSeparatedPseudoGenome(sPg, skipPgSequence);
        builder.build(pseudoGenomePrefix);
        cout << "Writing (" << pseudoGenomePrefix << ") pseudo genome files in " << time_millis() << " msec." << endl << endl;
    }

    void SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(SeparatedPseudoGenome *sPg,
                                                                      ostream* pgrcOut, uint8_t coder_level,
                                                                      bool ignoreOffDest, bool skipPgSequence) {
        time_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(!sPg->getReadsList()->isRevCompEnabled(),
                                                   !sPg->getReadsList()->areMismatchesEnabled());
        builder.feedSeparatedPseudoGenome(sPg, true);
        builder.compressedBuild(*pgrcOut, coder_level, ignoreOffDest);
        *logout << "Compressed Pg reads list in " << time_millis() << " msec." << endl << endl;
        if (!skipPgSequence) {
            time_checkpoint();
            auto pgCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, LZMA_DATAPERIODCODE_8_t);
            writeCompressed(*pgrcOut, sPg->getPgSequence().data(), sPg->getPgSequence().size(), pgCoderProps.get(),
                    COMPRESSION_ESTIMATION_BASIC_DNA);
            *logout << "Compressed Pg sequence in " << time_millis() << " msec." << endl << endl;
        }
    }

    SeparatedPseudoGenome* SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(const string &pgPrefix,
            bool skipReadsList){
        PseudoGenomeHeader* pgh = nullptr;
        ReadsSetProperties* prop = nullptr;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string pgSequence = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(pgPrefix);
        ExtendedReadsListWithConstantAccessOption* caeRl = skipReadsList?nullptr:
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
        PgHelpers::readArray(pgSrc, &buffer[0], size);
        return buffer;
    }


    void SeparatedPseudoGenomePersistence::writePseudoGenomeSequence(string &pgSequence, string pgPrefix) {
        ofstream pgDest =
                getPseudoGenomeElementDest(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX, true);
        PgHelpers::writeArray(pgDest, (void*) pgSequence.data(), pgSequence.length());
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
        writeReadMode(pair1OffsetsDest, PgHelpers::plainTextWriteMode);
        writeReadMode(pair1IndexesDest, PgHelpers::plainTextWriteMode);
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
        time_checkpoint();
        vector<uint_reads_cnt_std> orgIdxs;
        for(string pgFilePrefix: pgFilePrefixes)
            SeparatedPseudoGenomePersistence::appendIndexesFromPg(pgFilePrefix, orgIdxs);

        SeparatedPseudoGenomePersistence::writePairMapping(pgFilePrefixes[0], orgIdxs);
        *logout << "... dumping pairs completed in " << time_millis() << " msec. " << endl;
    }

    void SeparatedPseudoGenomePersistence::compressReadsOrder(ostream &pgrcOut,
            const vector<uint_reads_cnt_std>& orgIdxs, uint8_t coder_level,
            bool completeOrderInfo, bool ignorePairOrderInformation, bool singleFileMode) {
        time_checkpoint();
        uint_reads_cnt_std readsCount = orgIdxs.size();
        int lzma_reads_dataperiod_param = readsCount <= UINT32_MAX ? LZMA_DATAPERIODCODE_32_t : LZMA_DATAPERIODCODE_64_t;
        vector<uint_reads_cnt_std> rev(readsCount);
        for (uint_reads_cnt_std i = 0; i < readsCount; i++)
            rev[orgIdxs[i]] = i;
        if (completeOrderInfo && singleFileMode) {
            *logout << "Reverse index of original indexes... ";
            auto revIdxCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, lzma_reads_dataperiod_param);
            writeCompressed(pgrcOut, (char *) rev.data(), rev.size() * sizeof(uint_reads_cnt_std),
                    revIdxCoderProps.get());
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
            vector<CompressionJob> cJobs;
            auto rlRelOffFlagCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 11);
            cJobs.emplace_back("Uint8 reads list relative offsets of pair reads (flag)... ",
                               (unsigned char *) offsetInUint8Flag.data(), offsetInUint8Flag.size() * sizeof(uint8_t),
                               rlRelOffFlagCoderProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
            auto rlRelOffValCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 2);
            cJobs.emplace_back("Uint8 reads list relative offsets of pair reads (value)... ",
                               (unsigned char *) offsetInUint8Value.data(), offsetInUint8Value.size() * sizeof(uint8_t),
                            rlRelOffValCoderProps.get());
            auto relOffDelFlagCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 10);
            cJobs.emplace_back("Relative offsets deltas of pair reads (flag)... ",
                            (unsigned char *) deltaInInt8Flag.data(), deltaInInt8Flag.size() * sizeof(uint8_t),
                            relOffDelFlagCoderProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
            auto relOffDelValCoderProps = getRelativeOffsetDeltasOfPairsValueCoderProps(coder_level);
            cJobs.emplace_back("Relative offsets deltas of pair reads (value)... ",
                            (unsigned char*) deltaInInt8Value.data(), deltaInInt8Value.size() * sizeof(uint8_t),
                            relOffDelValCoderProps.get());
//            writeCompressed(pgrcOut, (char *) deltaInInt8Value.data(), deltaInInt8Value.size() * sizeof(int8_t), PPMD7_CODER, coder_level, 2);
            auto rlRelOffCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, lzma_reads_dataperiod_param);
            double estimated_reads_ratio = simpleUintCompressionEstimate(readsCount, readsCount <= UINT32_MAX?UINT32_MAX:UINT64_MAX);
            cJobs.emplace_back("Full reads list relative offsets of pair reads ... ",
                               (unsigned char *) fullOffset.data(), fullOffset.size() * sizeof(uint_reads_cnt_std),
                               rlRelOffCoderProps.get(), estimated_reads_ratio);
            auto basesOrgIdxsCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, lzma_reads_dataperiod_param);
            auto fileFlagsOffCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 2);
            auto fileFlagsNoOffCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 10);
            if (completeOrderInfo) {
                cJobs.emplace_back("Original indexes of pair bases... ",
                                   (unsigned char *) revPairBaseOrgIdx.data(), revPairBaseOrgIdx.size() * sizeof(uint_reads_cnt_std),
                                basesOrgIdxsCoderProps.get(), estimated_reads_ratio);
            } else if (!ignorePairOrderInformation) {
                cJobs.emplace_back("File flags of pair bases (for offsets)... ",
                                   (unsigned char *) offsetPairBaseFileFlag.data(), offsetPairBaseFileFlag.size() * sizeof(uint8_t),
                                fileFlagsOffCoderProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
                cJobs.emplace_back("File flags of pair bases (for non-offsets)... ",
                                   (unsigned char *) nonOffsetPairBaseFileFlag.data(), nonOffsetPairBaseFileFlag.size() * sizeof(uint8_t),
                                fileFlagsNoOffCoderProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
            }
            CompressionJob::writeCompressedCollectiveParallel(pgrcOut, cJobs);
        }
        *logout << "... compressing order information completed in " << time_millis() << " msec. " << endl;
        *logout << endl;
    }

    void SeparatedPseudoGenomePersistence::decompressReadsOrder(istream &pgrcIn,
                                                                vector<uint_reads_cnt_std> &rlIdxOrder,
                                                                bool completeOrderInfo, bool ignorePairOrderInformation,
                                                                bool singleFileMode) {
        if (singleFileMode) {
            if (!completeOrderInfo)
                return;
            string rlIdxOrderStr;
            readCompressed(pgrcIn, rlIdxOrderStr);
            moveStringToVector(rlIdxOrderStr, rlIdxOrder);
        } else {
            vector<uint8_t> offsetInUint8Flag;
            vector<uint8_t> offsetInUint8Value;
            vector<uint8_t> deltaInInt8Flag;
            vector<int8_t> deltaInInt8Value;
            vector<uint_reads_cnt_std> fullOffset;

            string offsetInUint8FlagStr, offsetInUint8ValueStr, deltaInInt8FlagStr, deltaInInt8ValueStr;
            string fullOffsetStr, revPairBaseOrgIdxStr, offsetPairBaseFileFlagStr, nonOffsetPairBaseFileFlagStr;
            vector<string*> destStrings;
            destStrings.push_back(&offsetInUint8FlagStr);
            destStrings.push_back(&offsetInUint8ValueStr);
            destStrings.push_back(&deltaInInt8FlagStr);
            destStrings.push_back(&deltaInInt8ValueStr);
            destStrings.push_back(&fullOffsetStr);

            if (completeOrderInfo) {
                destStrings.push_back(&revPairBaseOrgIdxStr);
            } else if (!ignorePairOrderInformation) {
                destStrings.push_back(&offsetPairBaseFileFlagStr);
                destStrings.push_back(&nonOffsetPairBaseFileFlagStr);
            }
            readCompressedCollectiveParallel(pgrcIn, destStrings);

            moveStringToVector(offsetInUint8FlagStr, offsetInUint8Flag);
            moveStringToVector(offsetInUint8ValueStr, offsetInUint8Value);
            moveStringToVector(deltaInInt8FlagStr, deltaInInt8Flag);
            moveStringToVector(deltaInInt8ValueStr, deltaInInt8Value);
            moveStringToVector(fullOffsetStr, fullOffset);
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
                moveStringToVector(revPairBaseOrgIdxStr, revPairBaseOrgIdx);
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
                moveStringToVector(offsetPairBaseFileFlagStr, offsetPairBaseFileFlag);
                moveStringToVector(nonOffsetPairBaseFileFlagStr, nonOffsetPairBaseFileFlag);
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
        time_checkpoint();
        uint_reads_cnt_std readsTotalCount = orgIdx2PgPos.size();
        int lzma_pos_dataperiod_param = sizeof(uint_pg_len) == 4 ? LZMA_DATAPERIODCODE_32_t : LZMA_DATAPERIODCODE_64_t;
        double estimated_pos_ratio = simpleUintCompressionEstimate(joinedPgLength, sizeof(uint_pg_len) == 4?UINT32_MAX:UINT64_MAX);
        if (singleFileMode) {
            uint_pg_len_max* const maxPgPosPtr = orgIdx2PgPos.data();
            if (sizeof(uint_pg_len) < sizeof(uint_pg_len_max)) {
                uint_pg_len* PgPosPtr = (uint_pg_len*) maxPgPosPtr;
                for (uint_reads_cnt_std i = 0; i < readsTotalCount; i++)
                    *(PgPosPtr++) = (uint_pg_len) orgIdx2PgPos[i];
            }
            auto readsPgPosProps = getReadsPositionsCoderProps(coder_level, lzma_pos_dataperiod_param);
            writeCompressed(pgrcOut, (char*) maxPgPosPtr, readsTotalCount * sizeof(uint_pg_len),
                            readsPgPosProps.get(), estimated_pos_ratio);
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
            parallel_algorithm::stable_sort(bppRank.begin(), bppRank.end(),
                    [&](const uint_reads_cnt_std &idx1, const uint_reads_cnt_std &idx2) -> bool
                        { return basePairPos[idx1] < basePairPos[idx2]; });
            *logout << "... reordering bases checkpoint: " << time_millis() << " msec. " << endl;
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
            assert(deltaPairEncodingEnabled);
            vector<CompressionJob> cJobs;
            auto basePosProps = getReadsPositionsCoderProps(coder_level, lzma_pos_dataperiod_param);
            cJobs.emplace_back("Base pair position... ", (unsigned char *) basePairPos.data(), basePairPos.size() * sizeof(uint_pg_len),
                            basePosProps.get(), estimated_pos_ratio);
            auto pairRelOffFlagProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 11);
            cJobs.emplace_back("Uint16 relative offset of pair positions (flag)... ",
                               (unsigned char *) offsetInUint16Flag.data(), offsetInUint16Flag.size() * sizeof(uint8_t),
                            pairRelOffFlagProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
            auto ppmd2CoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 2);
            auto fse12CoderProps = getDefaultFSECoderProps(12);
            auto pairRelOffPosFlagProps = getSelectorCoderProps( { ppmd2CoderProps.get(), fse12CoderProps.get() } );
            cJobs.emplace_back("Is uint16 relative offset of pair positions positive (flag)... ",
                               (unsigned char *) offsetIsBaseFirstFlag.data(), offsetIsBaseFirstFlag.size() * sizeof(uint8_t),
                            pairRelOffPosFlagProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
            auto pairRelOffValProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 2);
            auto range2CoderProps = getDefaultRangeCoderProps(256, 2);
            auto selPairRelOffValCoderProps = getSelectorCoderProps({ range2CoderProps.get(), pairRelOffValProps.get() });
            cJobs.emplace_back("Uint16 relative offset of pair positions (value)... ",
                               (unsigned char*) offsetInUint16Value.data(), offsetInUint16Value.size() * sizeof(uint16_t),
                            selPairRelOffValCoderProps.get());
            auto pairRelOffDeltaFlagProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 10);
            auto pairRelOffDeltaPosFlagProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 5);
            auto pairRelOffDeltaValProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 3);
            if (deltaPairEncodingEnabled) {
                cJobs.emplace_back("Relative offset deltas of pair positions (flag)... ",
                                   (unsigned char *) deltaInInt16Flag.data(), deltaInInt16Flag.size() * sizeof(uint8_t),
                                pairRelOffDeltaFlagProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
                cJobs.emplace_back("Is relative offset (for deltas stream) of pair positions positive (flag)... ",
                                   (unsigned char *) deltaIsBaseFirstFlag.data(),
                                deltaIsBaseFirstFlag.size() * sizeof(uint8_t),
                                pairRelOffDeltaPosFlagProps.get(), COMPRESSION_ESTIMATION_UINT8_BITMAP);
                cJobs.emplace_back("Relative offset deltas of pair positions (value)... ",
                                   (unsigned char *) deltaInInt16Value.data(), deltaInInt16Value.size() * sizeof(int16_t),
                                pairRelOffDeltaValProps.get());
            }
            auto nonBasePosProps = getReadsPositionsCoderProps(coder_level, lzma_pos_dataperiod_param);
            cJobs.emplace_back("Not-base pair position... ", (unsigned char *) notBasePairPos.data(), notBasePairPos.size() * sizeof(uint_pg_len),
                            nonBasePosProps.get(), estimated_pos_ratio);
            CompressionJob::writeCompressedCollectiveParallel(pgrcOut, cJobs);
        }
        *logout << "... compressing reads positions completed in " << time_millis() << " msec. " << endl;
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
                                                                          PgRCParams* params) {
        uint_reads_cnt_std readsTotalCount = params->readsTotalCount;
        bool singleFileMode = params->singleReadsMode;
        string pgPosStr;
        if (singleFileMode) {
            readCompressed(pgrcIn, pgPosStr);
            moveStringToVector(pgPosStr, pgPos);
        } else {
            bool deltaPairEncodingEnabled = params->isVersionAtLeast(1, 3) ? true : (bool) pgrcIn.get();
            vector<uint8_t> offsetInUint16Flag;
            vector<uint8_t> offsetIsBaseFirstFlag;
            vector<uint16_t> offsetInUint16Value;
            vector<uint8_t> deltaInInt16Flag;
            vector<uint8_t> deltaIsBaseFirstFlag;
            vector<int16_t> deltaInInt16Value;
            vector<uint_pg_len> notBasePairPos;
            pgPos.reserve(readsTotalCount);

            string offsetInUint16FlagStr, offsetIsBaseFirstFlagStr, offsetInUint16ValueStr;
            string deltaInInt16FlagStr, deltaIsBaseFirstFlagStr, deltaInInt16ValueStr, notBasePairPosStr;
            vector<string*> destStrings;
            destStrings.push_back(&pgPosStr);
            destStrings.push_back(&offsetInUint16FlagStr);
            destStrings.push_back(&offsetIsBaseFirstFlagStr);
            destStrings.push_back(&offsetInUint16ValueStr);
            if (deltaPairEncodingEnabled) {
                destStrings.push_back(&deltaInInt16FlagStr);
                destStrings.push_back(&deltaIsBaseFirstFlagStr);
                destStrings.push_back(&deltaInInt16ValueStr);
            }
            destStrings.push_back(&notBasePairPosStr);
            readCompressedCollectiveParallel(pgrcIn, destStrings);
            moveStringToVector(pgPosStr, pgPos);
            moveStringToVector(offsetInUint16FlagStr, offsetInUint16Flag);
            moveStringToVector(offsetIsBaseFirstFlagStr, offsetIsBaseFirstFlag);
            moveStringToVector(offsetInUint16ValueStr, offsetInUint16Value);
            if (deltaPairEncodingEnabled) {
                moveStringToVector(deltaInInt16FlagStr, deltaInInt16Flag);
                moveStringToVector(deltaIsBaseFirstFlagStr, deltaIsBaseFirstFlag);
                moveStringToVector(deltaInInt16ValueStr, deltaInInt16Value);
            }
            moveStringToVector(notBasePairPosStr, notBasePairPos);

            const uint_reads_cnt_std pairsCount = readsTotalCount / 2;
            vector<uint_reads_cnt_std> bppRank;
            bppRank.reserve(pairsCount);
            for (uint_reads_cnt_std p = 0; p < pairsCount; p++)
                bppRank.push_back(p);
            parallel_algorithm::stable_sort(bppRank.begin(), bppRank.end(),
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
    template void SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_std>(istream &pgrcIn, vector<uint_pg_len_std> &pgPos, PgRCParams* params);
    template void SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_max>(istream &pgrcIn, vector<uint_pg_len_max> &pgPos, PgRCParams* params);

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
        if (dest == nullptr) {
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
            PgHelpers::writeStringToFile(fileName, ((ostringstream*) dest)->str());
    }

    void SeparatedPseudoGenomeOutputBuilder::freeDest(ostream* &dest) {
        if (dest) {
            if (onTheFlyMode())
                ((ofstream*) dest)->close();
            delete(dest);
            dest = nullptr;
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
        if (pgh == nullptr) {
            fprintf(stderr, "Pseudo genome header not initialized in separated Pg builder.\n");
            exit(EXIT_FAILURE);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::buildProps() {
        pgh->setReadsCount(readsCounter);
        rsProp->readsCount = readsCounter;
        rsProp->allReadsLength = rsProp->constantReadLength ? (size_t) readsCounter * rsProp->maxReadLength : -1;
        if (pgPropDest) {
            delete (pgPropDest);
            pgPropDest = nullptr;
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
        PgHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_PROPERTIES_SUFFIX,
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

    string SeparatedPseudoGenomeOutputBuilder::toStringAndSeparateZeros(ostream* dest, string& zeroFlags) {
        string tmp = ((ostringstream*) dest)->str();
        zeroFlags.reserve(tmp.size());
        size_t j = 0;
        for (size_t i = 0; i < tmp.size(); i++) {
            bool isZero = tmp[i] == 0;
            zeroFlags.push_back(isZero);
            if (!isZero)
                tmp[j++] = tmp[i];
        }
        tmp.resize(j);
        return tmp;
    }

    string SeparatedPseudoGenomeOutputBuilder::toString(ostream* dest) {
        if (onTheFlyMode() || !dest) {
            fprintf(stderr, "Error during compression: an input stream missing.\n");
            exit(EXIT_FAILURE);
        }
        return ((ostringstream*) dest)->str();
    }

    void SeparatedPseudoGenomeOutputBuilder::compressRlMisRevOffDest(ostream &pgrcOut, uint8_t coder_level,
                                                                     bool separateFirstOffsetMode, bool transposeMode) {
        ostringstream mismatchesPropsOut;
        auto fseCoderProps = getDefaultFSECoderProps();
        uint8_t mismatches_dests_count = coder_level == CODER_LEVEL_FAST?1:(UINT8_MAX-1);
        if (mismatches_dests_count == 1) {
            PgHelpers::writeValue<uint8_t>(mismatchesPropsOut, 1);
            string mismatchesPropsString(mismatchesPropsOut.str());
            writeCompressed(pgrcOut, mismatchesPropsString, fseCoderProps.get());
            auto coderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 3);
            string tmp = toString(rlMisRevOffDest);
            writeCompressed(pgrcOut, tmp, coderProps.get());
            return;
        }
        vector<uint8_t> misCnt2DestIdx; // NOT-TESTED => {0, 1, 2, 3, 4, 5, 6, 7, 7, 9, 9, 9 };
        misCnt2DestIdx.insert(misCnt2DestIdx.end(), UINT8_MAX, mismatches_dests_count);
        for(uint8_t m = 1; m < mismatches_dests_count; m++)
            misCnt2DestIdx[m] = m;

        ostringstream destsFirst[UINT8_MAX];
        ostringstream dests[UINT8_MAX];
        istringstream misRevOffSrc(((ostringstream*) rlMisRevOffDest)->str());
        istringstream misCntSrc(((ostringstream*) rlMisCntDest)->str());

        uint8_t misCnt = 0;
        uint16_t revOff = 0;
        for(uint_reads_cnt_max i = 0; i < readsCounter; i++) {
             PgHelpers::readValue<uint8_t>(misCntSrc, misCnt, false);
             for(uint8_t m = 0; m < misCnt; m++) {
                 PgHelpers::readReadLengthValue(misRevOffSrc, revOff, false);
                 PgHelpers::writeReadLengthValue(dests[misCnt2DestIdx[misCnt]], revOff);
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
        PgHelpers::writeValue<uint8_t>(mismatchesPropsOut, mismatches_dests_count);
        for(uint8_t m = 1; m < mismatches_dests_count; m++)
            PgHelpers::writeValue<uint8_t>(mismatchesPropsOut, misCnt2DestIdx[m]);
        string mismatchesPropsString(mismatchesPropsOut.str());

        auto ppmdCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 2);
        vector<unique_ptr<CoderProps>> rangeCoderPropsList(mismatches_dests_count + 1);
        vector<unique_ptr<CoderProps>> firstCoderPropsList(mismatches_dests_count + 1);
        vector<unique_ptr<CoderProps>> seleCoderPropsList(mismatches_dests_count + 1);
        vector<string> misRefOffFirst(mismatches_dests_count + 1);
        vector<string> misRefOff(mismatches_dests_count + 1);
        vector<CompressionJob> cJobs;
        cJobs.emplace_back("mismatches props... ", mismatchesPropsString, fseCoderProps.get());
        for(uint8_t m = 1; m <= mismatches_dests_count; m++) {
            rangeCoderPropsList[m] = getDefaultRangeCoderProps(rsProp->maxReadLength,
                                                               separateFirstOffsetMode ? 1 : m);
            firstCoderPropsList[m] = getSelectorCoderProps(
                    {rangeCoderPropsList[m].get(), fseCoderProps.get(), ppmdCoderProps.get()}, 0.2);
            seleCoderPropsList[m] = getSelectorCoderProps(
                    {rangeCoderPropsList[m].get(), fseCoderProps.get(), ppmdCoderProps.get()}, 0.2);
        }
        for(uint8_t m = 1; m <= mismatches_dests_count; m++) {
            if (separateFirstOffsetMode) {
                misRefOffFirst[m] = toString(&destsFirst[m]);
                cJobs.emplace_back(to_string((int) m) + " (1st): ", misRefOffFirst[m], firstCoderPropsList[m].get());
            }
            if (!separateFirstOffsetMode || m > 1) {
                misRefOff[m] = toString(&dests[m]);
                cJobs.emplace_back(to_string((int) m) + ": ", misRefOff[m], seleCoderPropsList[m].get());
            }
        }
        CompressionJob::writeCompressedCollectiveParallel(pgrcOut, cJobs);
    }

    void SeparatedPseudoGenomeOutputBuilder::compressedBuild(ostream &pgrcOut, uint8_t coder_level, bool ignoreOffDest) {
        prebuildAssert(false);
        buildProps();
        writeReadMode(*pgPropDest, false);
        const string pgPropStr = ((ostringstream*) pgPropDest)->str();
        vector<CompressionJob> cJobs;
        auto fse12CoderProps = getDefaultFSECoderProps(12);
        cJobs.emplace_back("Pseudogenome props... ", pgPropStr, fse12CoderProps.get());
        string rlOffStr, rlRevCompStr, zeroFlags, misCnts, rlMisSym;
        auto ppmd5CoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 5);
        auto rloCoderProps = getSelectorCoderProps( {fse12CoderProps.get(), ppmd5CoderProps.get() });
        if (!ignoreOffDest) {
            rlOffStr = toString(rlOffDest);
            cJobs.emplace_back("Reads list offsets... ", rlOffStr, rloCoderProps.get());
        }
        auto rciCoderProps = ignoreOffDest ? getDefaultFSECoderProps( 12) :
                getDefaultCoderProps(PPMD7_CODER, coder_level, 3);
        if (!this->disableRevComp) {
            rlRevCompStr = toString(rlRevCompDest);
            cJobs.emplace_back("Reverse complements info... ", rlRevCompStr, rciCoderProps.get(),
                               COMPRESSION_ESTIMATION_UINT8_BITMAP);
        }
        auto rangeCoderProps = getDefaultRangeCoderProps();
        auto ppmd13CoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 13);
        auto zerosFlagCoderProps = getSelectorCoderProps( { fse12CoderProps.get(),
                                                            ignoreOffDest ? rangeCoderProps.get() : ppmd13CoderProps.get() } );;
        auto mismatchesCountsCoderProps = getSelectorCoderProps( { fse12CoderProps.get(), rangeCoderProps.get() } );
        auto mismatchesSymbolsCoderProps = getDefaultCoderProps(PPMD7_CODER, coder_level, 8);
        if (!this->disableMismatches) {
            misCnts = toStringAndSeparateZeros(rlMisCntDest, zeroFlags);
            cJobs.emplace_back("Mismatches counts (zero flags)... ", zeroFlags, zerosFlagCoderProps.get(),
                               COMPRESSION_ESTIMATION_MIS_CNT);
            cJobs.emplace_back("Mismatches counts (non-zero values)... ", misCnts,
                               mismatchesCountsCoderProps.get(),
                               COMPRESSION_ESTIMATION_MIS_CNT);
            rlMisSym = toString(rlMisSymDest);
            string basesOrdered = reorderingSymbolsExclusiveMismatchEncoding(rlMisSym);
            pgrcOut.write(basesOrdered.data(), basesOrdered.length());

            cJobs.emplace_back("Mismatched symbols codes... ", rlMisSym, mismatchesSymbolsCoderProps.get(),
                               COMPRESSION_ESTIMATION_MIS_SYM);
        }
        CompressionJob::writeCompressedCollectiveParallel(pgrcOut, cJobs);
        if (!this->disableMismatches) {
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
            PgHelpers::writeValue<uint_pg_len_max>(*rlPosDest, rlEntry.pos);
        else
            PgHelpers::writeReadLengthValue(*rlOffDest, rlEntry.offset);
        PgHelpers::writeValue<uint_reads_cnt_std>(*rlOrgIdxDest, rlEntry.idx);
        if (!disableRevComp)
            PgHelpers::writeValue<uint8_t>(*rlRevCompDest, rlEntry.revComp?1:0);
        if (!disableMismatches) {
            PgHelpers::writeValue<uint8_t>(*rlMisCntDest, rlEntry.mismatchesCount);
            if (rlEntry.mismatchesCount) {
                for (uint8_t i = 0; i < rlEntry.mismatchesCount; i++)
                    PgHelpers::writeValue<uint8_t>(*rlMisSymDest, rlEntry.mismatchCode[i]);
                if (SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation) {
                    uint8_t currentPos = pgh->getMaxReadLength() - 1;
                    for (int16_t i = rlEntry.mismatchesCount - 1; i >= 0; i--) {
                        PgHelpers::writeReadLengthValue(*rlMisRevOffDest,
                                                                   currentPos - rlEntry.mismatchOffset[i]);
                        currentPos = rlEntry.mismatchOffset[i] - 1;
                    }
                } else {
                    for (uint8_t i = 0; i < rlEntry.mismatchesCount; i++)
                        PgHelpers::writeReadLengthValue(*rlMisPosDest, rlEntry.mismatchOffset[i]);
                }
            }
        }
        readsCounter++;
    }

    void SeparatedPseudoGenomeOutputBuilder::writeExtraReadEntry(const DefaultReadsListEntry &rlEntry) {
        writeReadEntry(rlEntry);
        if (rlIt != nullptr) {
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
        PgHelpers::writeArray(*pgDest, (void*) pg.data(), pg.length());

        ReadsListIteratorExtendedWrapperBase* rlIt =
                TemplateUserGenerator::generateReadsListUser<ReadsListIteratorExtendedWrapper, ReadsListIteratorExtendedWrapperBase>(pgb);
        rlIt->applyIndexesMapping(orgIndexesMapping);
        if (revComplPairFile)
            rlIt->applyRevComplPairFileFlag();
        if (pgb->isReadLengthMin())
            bytePerReadLengthMode = true;
        setReadsSourceIterator(rlIt);
        writeReadsFromIterator();

        if (pgh == nullptr)
            pgh = new PseudoGenomeHeader(pgb);
        if (rsProp == nullptr)
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
        if (pgh != nullptr)
            delete(pgh);
        if (rsProp != nullptr)
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
            PgHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        } else {
            pgh->setPseudoGenomeLength(pg.length());
            initDest(pgDest, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX);
            PgHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::feedSeparatedPseudoGenome(SeparatedPseudoGenome *sPg, bool skipPgSequence) {
        if (!skipPgSequence) {
            initDest(pgDest, SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX);
            const string &pg = sPg->getPgSequence();
            PgHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        }
        if (sPg->isReadLengthMin())
            bytePerReadLengthMode = true;
        setReadsSourceIterator(sPg->getReadsList());
        if (pgh == nullptr)
            pgh = new PseudoGenomeHeader(sPg);
        writeReadsFromIterator();

        if (rsProp == nullptr)
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

    string SeparatedPseudoGenomeOutputBuilder::reorderingSymbolsExclusiveMismatchEncoding(string &rlMisSym) {
        size_t counts[5] {};
        for(char c : rlMisSym) {
            uint8_t ctxCode = c;
            uint8_t mismatchValue = cxtCode2MismatchValue(ctxCode);
            counts[mismatchValue]++;
        }
        vector<uint8_t> order = { 0, 1, 2, 3, 4 };
        std::sort(order.begin(), order.end(), [counts](const uint8_t &o1, const uint8_t &o2) -> bool
                                { return counts[o1] > counts[o2]; });
        uint8_t orderRevIdx[5];
        string basesOrder;
        for(int i = 0; i < 5; i++) {
            basesOrder.push_back(value2symbol(order[i]));
            orderRevIdx[order[i]] = i;
        }
        for (char &c : rlMisSym) {
            uint8_t ctxCode = c;
            uint8_t actualValue = orderRevIdx[cxtCode2ActualValue(ctxCode)];
            uint8_t mismatchValue = orderRevIdx[cxtCode2MismatchValue(ctxCode)];
            c = (char) (mismatchValue - (mismatchValue > actualValue?1:0));
        }
        return basesOrder;
    }
}
