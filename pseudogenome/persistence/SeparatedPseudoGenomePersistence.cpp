#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"

#include <cstdio>

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
            const string &pseudoGenomePrefix, ostream* pgrcOut, bool skipPgSequence) {
        clock_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(sPg->getReadsList()->isRevCompEnabled(),
                sPg->getReadsList()->areMismatchesEnabled());
        builder.writeSeparatedPseudoGenome(sPg, skipPgSequence);
        if (!pseudoGenomePrefix.empty()) builder.build(pseudoGenomePrefix);
        if (pgrcOut) builder.compressedBuild(*pgrcOut);
        cout << "Writing (" << pseudoGenomePrefix << ") pseudo genome files in " << clock_millis() << " msec." << endl << endl;
    }

    SeparatedPseudoGenome* SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(const string &pgPrefix,
            bool skipReadsList){
        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string pgSequence = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(pgPrefix);
        ConstantAccessExtendedReadsList* caeRl = skipReadsList?0:
                ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(pgPrefix, pgh->getPseudoGenomeLength());
        SeparatedPseudoGenome* spg = new SeparatedPseudoGenome(std::move(pgSequence), caeRl, prop);
        delete(pgh);
        delete(prop);
        return spg;
    }


    std::ifstream SeparatedPseudoGenomePersistence::getPseudoGenomeSrc(const string &pseudoGenomePrefix) {
        const string pgFile = pseudoGenomePrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX;
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
                getPseudoGenomeElementDest(pgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX, true);
        PgSAHelpers::writeArray(pgDest, (void*) pgSequence.data(), pgSequence.length());
        pgDest.close();
        acceptTemporaryPseudoGenomeElement(pgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX, true);
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
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, PSEUDOGENOME_FILE_SUFFIX, false);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, PSEUDOGENOME_PROPERTIES_SUFFIX, false);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_POSITIONS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_REVERSECOMPL_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHES_COUNT_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_OFFSETS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX, clearReadsListDescriptionIfOverride);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX, clearReadsListDescriptionIfOverride);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX, false);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX, false);
    }

    const string SeparatedPseudoGenomePersistence::PG_FILES_EXTENSION = ".pg";

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX = "" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX = "_prop" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_POSITIONS_FILE_SUFFIX = "_rl_pos" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = "_rl_idx" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX = "_rl_rc" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX = "_rl_mis_cnt" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX = "_rl_mis_sym" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX = "_rl_mis_pos" + PG_FILES_EXTENSION;


    const string SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX = "_rl_off" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX = "_rl_mis_roff" + PG_FILES_EXTENSION;

    const string SeparatedPseudoGenomePersistence::READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX = "_rl_pr_idx" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX = "_rl_pr_off" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX = "_rl_pr_sf" + PG_FILES_EXTENSION;

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX = "_map_off" + PG_FILES_EXTENSION;
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX = "_map_len" + PG_FILES_EXTENSION;

    const string SeparatedPseudoGenomePersistence::TEMPORARY_FILE_SUFFIX = ".temp";

    bool SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = false;
    bool SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = true;

    void SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(const string &pseudoGenomePrefix,
                                                                     PseudoGenomeHeader *&pgh,
                                                                     ReadsSetProperties *&rsProp,
                                                                     bool &plainTextReadMode) {
        ifstream pgPropSrc(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX,
                           ios_base::in | ios_base::binary);
        if (pgPropSrc.fail()) {
            fprintf(stderr, "Cannot read pseudogenome properties (%s does not open).\n",
                    pseudoGenomePrefix.c_str());
            exit(EXIT_FAILURE);
        }
        pgh = new PseudoGenomeHeader(pgPropSrc);
        rsProp = new ReadsSetProperties(pgPropSrc);
        plainTextReadMode = confirmTextReadMode(pgPropSrc);
        pgPropSrc.close();
    }

    void SeparatedPseudoGenomePersistence::appendIndexesFromPg(string pgFilePrefix, vector<uint_reads_cnt_std> &idxs) {
        bool plainTextReadMode;
        PseudoGenomeHeader* pgh;
        ReadsSetProperties* rsProp;
        getPseudoGenomeProperties(pgFilePrefix, pgh, rsProp, plainTextReadMode);
        ifstream orgIdxsSrc = getPseudoGenomeElementSrc(pgFilePrefix, READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
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
        ofstream pair1OffsetsDest = getPseudoGenomeElementDest(pgFilePrefix, READSLIST_PAIR_FIRST_OFFSETS_FILE_SUFFIX, true);
        ofstream pair1SrcFlagDest = getPseudoGenomeElementDest(pgFilePrefix, READSLIST_PAIR_FIRST_SOURCE_FLAG_FILE_SUFFIX, true);
        ofstream pair1IndexesDest = getPseudoGenomeElementDest(pgFilePrefix, READSLIST_PAIR_FIRST_INDEXES_FILE_SUFFIX, true);
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
        cout << "... dumping pairs completed in " << clock_millis() << " msec. " << endl;
    }

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
            initDest(rlPosDest, SeparatedPseudoGenomePersistence::READSLIST_POSITIONS_FILE_SUFFIX);
        else
            initDest(rlOffDest, SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX);
        initDest(rlOrgIdxDest, SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        if (!disableRevComp)
            initDest(rlRevCompDest, SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        if (!disableMismatches) {
            initDest(rlMisCntDest, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
            initDest(rlMisSymDest, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
            if (SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation)
                initDest(rlMisRevOffDest,
                         SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
            else
                initDest(rlMisPosDest, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::destToFile(ostream *dest, const string &fileName) {
        if (!onTheFlyMode() && dest)
            PgSAHelpers::writeStringToFile(fileName, ((ostringstream*) dest)->str());
    }

    void SeparatedPseudoGenomeOutputBuilder::compressDest(ostream* dest, ostream &pgrcOut, uint8_t coder_type,
            uint8_t coder_level, int coder_param) {
        if (onTheFlyMode() || !dest) {
            fprintf(stderr, "Error during compression: an input stream missing.\n");
            exit(EXIT_FAILURE);
        }
        const string tmp = ((ostringstream*) dest)->str();
        writeCompressed(pgrcOut, tmp, coder_type, coder_level, coder_param);
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
        initDest(pgPropDest, SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX);
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
        prebuildAssert(false);
        buildProps();
        writeReadMode(*pgPropDest, plainTextWriteMode);
        PgSAHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX,
                                      ((ostringstream*) pgPropDest)->str());

        destToFile(rlPosDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_POSITIONS_FILE_SUFFIX);
        destToFile(rlOffDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX);
        destToFile(rlOrgIdxDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        destToFile(rlRevCompDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX);
        destToFile(rlMisCntDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
        destToFile(rlMisSymDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
        destToFile(rlMisRevOffDest,
                   pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
        destToFile(rlMisPosDest, pgPrefix + SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);
    }

    void SeparatedPseudoGenomeOutputBuilder::compressedBuild(ostream &pgrcOut) {
        prebuildAssert(false);
        buildProps();
        writeReadMode(*pgPropDest, false);
        const string tmp = ((ostringstream*) pgPropDest)->str();
        pgrcOut.write(tmp.data(), tmp.length());

        int lzma_coder_param = this->rsProp->maxReadLength <= UINT8_MAX?PGRC_DATAPERIODCODE_8_t:PGRC_DATAPERIODCODE_16_t;
        cout << "Reads list offsets... ";
//        compressDest(rlOffDest, pgrcOut, LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM, lzma_coder_param);
        compressDest(rlOffDest, pgrcOut, PPMD7_CODER, PGRC_CODER_LEVEL_MAXIMUM, 3);
        cout << "Reverse complements info... ";
//        compressDest(rlRevCompDest, pgrcOut, LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM, PGRC_DATAPERIODCODE_8_t);
        compressDest(rlRevCompDest, pgrcOut, PPMD7_CODER, PGRC_CODER_LEVEL_MAXIMUM, 3);
        cout << "Mismatches counts... ";
//        compressDest(rlMisCntDest, pgrcOut, LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM, lzma_coder_param);
        compressDest(rlMisCntDest, pgrcOut, PPMD7_CODER, PGRC_CODER_LEVEL_MAXIMUM, 3);
        cout << "Mismatched symbols codes... ";
//        compressDest(rlMisSymDest, pgrcOut, LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM, lzma_coder_param);
        compressDest(rlMisSymDest, pgrcOut, PPMD7_CODER, PGRC_CODER_LEVEL_MAXIMUM, 3);
        cout << "Mismatches offsets (rev-coded)... ";
//        compressDest(rlMisRevOffDest, pgrcOut, LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM, lzma_coder_param);
        compressDest(rlMisRevOffDest, pgrcOut, PPMD7_CODER, PGRC_CODER_LEVEL_MAXIMUM, 3);
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

        initDest(pgDest, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX);
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
            fprintf(stderr, "Incorrect reads count validation while building separated Pg (%llu instead of %llu).\n",
                    readsCounter, pgh->getReadsCount());
            exit(EXIT_FAILURE);
        }

        delete(rlIt);
    }

    void SeparatedPseudoGenomeOutputBuilder::copyPseudoGenomeProperties(const string &pseudoGenomePrefix) {
        ifstream pgPropSrc(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX,
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
            initDest(pgDest, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX);
            PgSAHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::writeSeparatedPseudoGenome(SeparatedPseudoGenome *sPg, bool skipPgSequence) {
        if (!skipPgSequence) {
            initDest(pgDest, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX);
            const string &pg = sPg->getPgSequence();
            PgSAHelpers::writeArray(*pgDest, (void *) pg.data(), pg.length());
        }
        if (sPg->isReadLengthMin())
            bytePerReadLengthMode = true;
        setReadsSourceIterator(sPg->getReadsList());
        writeReadsFromIterator();

        if (pgh == 0)
            pgh = new PseudoGenomeHeader(sPg);
        if (rsProp == 0)
            rsProp = new ReadsSetProperties(*(sPg->getReadsSetProperties()));
        if (pgh->getReadsCount() != readsCounter) {
            fprintf(stderr, "Incorrect reads count validation while building separated Pg (%llu instead of %llu).\n",
                    readsCounter, pgh->getReadsCount());
            exit(EXIT_FAILURE);
        }
    }

    SeparatedPseudoGenomeOutputBuilder::~SeparatedPseudoGenomeOutputBuilder() {
        if (pgh)
            delete(pgh);
        if (rsProp)
            delete(rsProp);

    }
}
