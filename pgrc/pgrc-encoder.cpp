#include <sys/stat.h>
#include "pgrc-encoder.h"

#include "../matching/ReadsMatchers.h"
#include "../matching/SimplePgMatcher.h"
#include "../pseudogenome/TemplateUserGenerator.h"
#include "../readsset/persistance/ReadsSetPersistence.h"
#include "../pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "../pseudogenome/generator/ParallelGreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "../pseudogenome/persistence/PseudoGenomePersistence.h"
#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    uint_read_len_max probeReadsLength(const string &srcFastqFile);

    time_t getTimeInSec(chrono::steady_clock::time_point end_t, chrono::steady_clock::time_point begin_t) {
        chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(end_t - begin_t);
        return (double)time_span.count() / 1000000.0;
    }

    string getTimeInSec(chrono::steady_clock::time_point end_t, chrono::steady_clock::time_point begin_t, int decimalPlaces) {
        chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(end_t - begin_t);
        return toString((double)time_span.count() / 1000000000.0, decimalPlaces);
    }

    void PgRCEncoder::prepareChainData() {
        params->initCompressionParameters();
        params->qualityDivision = params->error_limit_in_promils < 1000;
        params->generatorDivision = params->gen_quality_coef > 0;
        params->readLength = probeReadsLength(params->srcFastqFile);
        if (params->pairFastqFile.empty() && !params->preserveOrderMode)
            params->singleReadsMode = true;

        if (std::ifstream(params->pgRCFileName))
            fprintf(stderr, "Warning: file %s already exists\n", params->pgRCFileName.data());
        pgrcOut = fstream(params->pgRCFileName + TEMPORARY_FILE_SUFFIX, ios::out | ios::binary | ios::trunc);
        pgrcOut.write(PGRC_HEADER, strlen(PGRC_HEADER));
        pgrcOut.put(PGRC_VERSION_MODE);
        pgrcOut.put(PGRC_VERSION_MAJOR);
        pgrcOut.put(PGRC_VERSION_MINOR);
        pgrcOut.put(PGRC_VERSION_REVISION);
        pgrcOut.put(params->compressionLevel);
        const char pgrc_mode = params->singleReadsMode ? PGRC_SE_MODE :
                               (params->preserveOrderMode?(params->pairFastqFile.empty()?PGRC_ORD_SE_MODE:PGRC_ORD_PE_MODE):
                     (params->ignorePairOrderInformation?PGRC_MIN_PE_MODE:PGRC_PE_MODE));
        pgrcOut.put(pgrc_mode);
        pgrcOut.put(params->separateNReads);
        params->revComplPairFile = false;
        if (pgrc_mode == PGRC_PE_MODE || pgrc_mode == PGRC_ORD_PE_MODE) {
            if (!params->disableRevComplPairFileMode)
                params->revComplPairFile = true;
            pgrcOut.put(params->revComplPairFile);
        }

        string tmpDirectoryName = params->pgRCFileName;
        if (params->qualityDivision)
            tmpDirectoryName = tmpDirectoryName + "_q" + toString(params->error_limit_in_promils) + (params->singleReadsMode?"s":"");
        tmpDirectoryName = tmpDirectoryName + (params->nReadsLQ?"_n":"") + (params->separateNReads?"_N":"") + "_g"
                + params->gen_quality_str;
        bool enablePreReadsMatching = params->preReadsExactMatchingChars > 0;
        tmpDirectoryName = tmpDirectoryName +
                (enablePreReadsMatching?("_l" + (((char) tolower(params->preMatchingMode))
                    + ((toupper(params->preMatchingMode) == params->preMatchingMode)?string("s"):string(""))
                    + toString(params->preReadsExactMatchingChars))):"")
                + "_m" + (((char) tolower(params->matchingMode))
                + (toupper(params->matchingMode) == params->matchingMode?string("s"):string(""))
                + toString(params->readsExactMatchingChars))
                + "_M" + toString(params->minCharsPerMismatch)
                + "_p" + toString(params->targetPgMatchLength);
        if (params->extraFilesForValidation && params->skipStages == 0) {
            mode_t mode = 0777;
#ifdef __MINGW32__
            int nError = mkdir(tmpDirectoryName.data());
#else
            int nError = mkdir(tmpDirectoryName.data(), mode);
#endif
            if (nError != 0) {
                srand(time(nullptr));
                tmpDirectoryName = tmpDirectoryName + "_" + toString(rand() % 100000);
#ifdef __MINGW32__
                int nError = mkdir(tmpDirectoryName.data());
#else
                nError = mkdir(tmpDirectoryName.data(), mode);
#endif
                if (nError != 0) {
                    fprintf(stderr, "Error creating folder %s\n", tmpDirectoryName.data());
                    exit(EXIT_FAILURE);
                }
            }
        }
        pgrcOut.write(tmpDirectoryName.data(), tmpDirectoryName.length());
        pgrcOut << endl;

        params->lqDivisionFile = tmpDirectoryName + "/" + BAD_INFIX + DIVISION_EXTENSION;
        params->nDivisionFile = tmpDirectoryName + "/" + N_INFIX + DIVISION_EXTENSION;
        params->pgHqPrefix = tmpDirectoryName + "/" + GOOD_INFIX;

        params->pgFilesPrefixesWithM = tmpDirectoryName + "/";
        params->pgMappedHqPrefix = params->pgFilesPrefixesWithM + GOOD_INFIX;
        params->pgMappedLqPrefix = params->pgFilesPrefixesWithM + BAD_INFIX;
        params->pgSeqFinalHqPrefix = params->pgFilesPrefixesWithM + GOOD_INFIX;
        params->pgSeqFinalLqPrefix = params->pgFilesPrefixesWithM + BAD_INFIX;
        params->pgNPrefix = params->pgFilesPrefixesWithM + N_INFIX;
        params->mappedLqDivisionFile = params->pgFilesPrefixesWithM + BAD_INFIX + DIVISION_EXTENSION;
    }

    void PgRCEncoder::executePgRCChain() {
        start_t = chrono::steady_clock::now();
        prepareChainData();
        uint8_t stageCount = 0;
        if (params->skipStages < ++stageCount && params->qualityDivision) {
            runQualityBasedDivision();
            if (params->disableInMemoryMode || params->endAtStage == stageCount) {
                persistReadsQualityDivision();
                data.disposeChainData();
            }
        }
        div_t = chrono::steady_clock::now();
        if (params->skipStages < ++stageCount && params->endAtStage >= stageCount) {
            prepareForPgGeneratorBaseReadsDivision();
            if (params->generatorDivision)
                runPgGeneratorBasedReadsDivision();
            if (params->disableInMemoryMode || params->endAtStage == stageCount) {
                persistReadsQualityDivision();
                data.disposeChainData();
            }

        }
        pgDiv_t = chrono::steady_clock::now();
        if (params->skipStages < ++stageCount && params->endAtStage >= stageCount) {
            prepareForHqPgGeneration();
            runHQPgGeneration();
            if (params->disableInMemoryMode || params->endAtStage == stageCount) {
                persistHQPg();
                persistReadsQualityDivision();
                data.disposeChainData();
            }
        }
        good_t = chrono::steady_clock::now();
        if (params->skipStages < ++stageCount && params->endAtStage >= stageCount) {
            prepareForMappingLQReadsOnHQPg();
            runMappingLQReadsOnHQPg();
            if (params->disableInMemoryMode || params->endAtStage == stageCount) {
                persistHQPgSequence();
                persistMappedReadsQualityDivision();
                data.disposeChainData();
            } else {
                //// Already done during runMappingLQReadsOnHQPg()
//                compressMappedHQPgReadsList();
                if (!params->singleReadsMode) {
                    if (params->preserveOrderMode)
                        data.orgIdx2PgPos = std::move(data.hqPg->getReadsList()->pos);
                    else
                        data.rlIdxOrder = std::move(data.hqPg->getReadsList()->orgIdx);
                }
                data.hqPg->disposeReadsList();
            }
        }
        match_t = chrono::steady_clock::now();
        if (params->skipStages < ++stageCount && params->endAtStage >= stageCount) {
            prepareForLQPgAndNPgGeneration();
            runLQPgGeneration();
            if (params->disableInMemoryMode || params->endAtStage == stageCount) {
                persistLQPg();
                if (data.hqPg)
                    persistHQPgSequence();
            } else {
                compressLQPgReadsList();
            }
            if (!params->singleReadsMode) {
                if (params->preserveOrderMode) {
                    ExtendedReadsListWithConstantAccessOption *const pgRl = data.lqPg->getReadsList();
                    uint_pg_len_max pos = data.hqPg->getPseudoGenomeLength();
                    for (uint_reads_cnt_std i = 0; i < pgRl->readsCount; i++) {
                        pos += pgRl->off[i];
                        data.orgIdx2PgPos[pgRl->orgIdx[i]] = pos;
                    }
                } else
                    data.rlIdxOrder.insert(data.rlIdxOrder.end(), data.lqPg->getReadsList()->orgIdx.begin(),
                                           data.lqPg->getReadsList()->orgIdx.end());
            }
            data.lqPg->disposeReadsList();
            if (params->separateNReads) {
                runNPgGeneration();
                if (params->disableInMemoryMode || params->endAtStage == stageCount) {
                    persistNPg();
                } else {
                    compressNPgReadsList();
                }
                if (!params->singleReadsMode) {
                    if (params->preserveOrderMode) {
                        ExtendedReadsListWithConstantAccessOption *const pgRl = data.nPg->getReadsList();
                        uint_pg_len_max pos = data.hqPg->getPseudoGenomeLength() + data.lqPg->getPseudoGenomeLength();
                        for (uint_reads_cnt_std i = 0; i < pgRl->readsCount; i++) {
                            pos += pgRl->off[i];
                            data.orgIdx2PgPos[pgRl->orgIdx[i]] = pos;
                        }
                    } else
                        data.rlIdxOrder.insert(data.rlIdxOrder.end(), data.nPg->getReadsList()->orgIdx.begin(),
                                               data.nPg->getReadsList()->orgIdx.end());
                }
                data.nPg->disposeReadsList();
            }
            delete(data.divReadsSets);
            data.divReadsSets = nullptr;
        }
        bad_t = chrono::steady_clock::now();
        if (!params->singleReadsMode && params->skipStages < ++stageCount && params->endAtStage >= stageCount) {
            if (params->preserveOrderMode) {
                const uint_pg_len_max joinedPgLength =
                        data.hqPg->getPseudoGenomeLength() + data.lqPg->getPseudoGenomeLength() +
                        data.nPg->getPseudoGenomeLength();
                bool isJoinedPgLengthStd = joinedPgLength <= UINT32_MAX;
                if (isJoinedPgLengthStd)
                    SeparatedPseudoGenomePersistence::compressReadsPgPositions<uint_pg_len_std>(pgrcOut,
                            data.orgIdx2PgPos, joinedPgLength, params->compressionLevel, params->pairFastqFile.empty());
                else
                    SeparatedPseudoGenomePersistence::compressReadsPgPositions<uint_pg_len_max>(pgrcOut,
                            data.orgIdx2PgPos, joinedPgLength, params->compressionLevel, params->pairFastqFile.empty());
                data.orgIdx2PgPos.clear();
            } else {
                SeparatedPseudoGenomePersistence::compressReadsOrder(pgrcOut, data.rlIdxOrder,
                        params->compressionLevel, params->preserveOrderMode, params->ignorePairOrderInformation,
                        params->pairFastqFile.empty());
                data.rlIdxOrder.clear();
            }
            if (params->extraFilesForValidation) {
                if (params->separateNReads)
                    SeparatedPseudoGenomePersistence::dumpPgPairs({params->pgMappedHqPrefix, params->pgMappedLqPrefix,
                                                                   params->pgNPrefix});
                else
                    SeparatedPseudoGenomePersistence::dumpPgPairs({params->pgMappedHqPrefix, params->pgMappedLqPrefix});
            }
        }
        order_t = chrono::steady_clock::now();
        if (params->skipStages < ++stageCount && params->endAtStage >= stageCount) {
            prepareForPgMatching();
            string emptySequence;
            SimplePgMatcher::matchPgsInPg(data.hqPg->getPgSequence(), data.lqPg->getPgSequence(),
                    params->separateNReads?data.nPg->getPgSequence():emptySequence, params->separateNReads,
                    pgrcOut, params->compressionLevel, params->extraFilesForValidation?params->pgSeqFinalHqPrefix:"",
                    params->extraFilesForValidation?params->pgSeqFinalLqPrefix:"",
                    params->extraFilesForValidation?params->pgNPrefix:"",
                    params->targetPgMatchLength);
        }
        finalizeCompression();
        data.disposeChainData();
#ifdef DEVELOPER_BUILD
        generateReport();
#endif
    }

    void PgRCEncoder::runQualityBasedDivision() {
        ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                params->srcFastqFile, params->pairFastqFile, params->revComplPairFile);
        data.divReadsSets =
                DividedPCLReadsSets::getQualityDivisionBasedReadsSets(allReadsIterator, params->readLength,
                        params->error_limit_in_promils / 1000.0, params->simplified_suffix_mode,
                        params->separateNReads, params->nReadsLQ);
        delete (allReadsIterator);
    }

    void PgRCEncoder::persistReadsQualityDivision() {
        data.divReadsSets->getLqReadsIndexesMapping()->saveMapping(params->lqDivisionFile);
        if (params->separateNReads)
            data.divReadsSets->getNReadsIndexesMapping()->saveMapping(params->nDivisionFile);
    }

    void PgRCEncoder::prepareForPgGeneratorBaseReadsDivision() {
        if (!data.divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    params->srcFastqFile, params->pairFastqFile, params->revComplPairFile);
            if (params->qualityDivision) {
                data.divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                        allReadsIterator, params->readLength, params->lqDivisionFile, params->nReadsLQ,
                        params->separateNReads ? params->nDivisionFile : "");
            } else {
                data.divReadsSets = DividedPCLReadsSets::getSimpleDividedPCLReadsSets(allReadsIterator, params->readLength,
                        params->separateNReads, params->nReadsLQ);
            }
            delete (allReadsIterator);
        }
    }

    void PgRCEncoder::runPgGeneratorBasedReadsDivision() {
        cout << "HQ ";
        data.divReadsSets->getHqReadsSet()->printout();
        const vector<bool>& isReadHqInHqReadsSet =
                (numberOfThreads > 1 && data.divReadsSets->getHqReadsSet()->readsCount() > PARALLEL_PG_GENERATION_READS_COUNT_THRESHOLD)?
                ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getHQReads(
                        data.divReadsSets->getHqReadsSet(), params->gen_quality_coef):
                GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getHQReads(
                        data.divReadsSets->getHqReadsSet(), params->gen_quality_coef);
        data.divReadsSets->moveLqReadsFromHqReadsSetsToLqReadsSets(isReadHqInHqReadsSet);
    }

    void PgRCEncoder::prepareForHqPgGeneration() {
        if (!data.divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    params->srcFastqFile, params->pairFastqFile, params->revComplPairFile);
            data.divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, params->readLength, params->lqDivisionFile, params->nReadsLQ,
                    params->separateNReads ? params->nDivisionFile : "");
            delete (allReadsIterator);
        }
    }

    void PgRCEncoder::runHQPgGeneration() {
        cout << "HQ ";
        data.divReadsSets->getHqReadsSet()->printout();
        if (numberOfThreads > 1 &&
            data.divReadsSets->getHqReadsSet()->readsCount() > PARALLEL_PG_GENERATION_READS_COUNT_THRESHOLD)
            data.hqPg = ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                    data.divReadsSets->getHqReadsSet());
        else
            data.hqPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                    data.divReadsSets->getHqReadsSet());
        data.divReadsSets->disposeHqReadsSet();
        IndexesMapping* hq2IndexesMapping = data.divReadsSets->generateHqReadsIndexesMapping();
        data.hqPg->applyIndexesMapping(hq2IndexesMapping);
        delete(hq2IndexesMapping);
    }

    void PgRCEncoder::persistHQPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(data.hqPg, params->pgHqPrefix);
    }

    void PgRCEncoder::prepareForMappingLQReadsOnHQPg() {
        if (!data.divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    params->srcFastqFile, params->pairFastqFile, params->revComplPairFile);
            data.divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, params->readLength, params->lqDivisionFile, params->nReadsLQ,
                    params->separateNReads ? params->nDivisionFile : "", true);
            delete (allReadsIterator);
        }
        if (!data.hqPg)
            data.hqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(params->pgHqPrefix);
    }

    void PgRCEncoder::runMappingLQReadsOnHQPg() {
        cout << "LQ ";
        data.divReadsSets->getLqReadsSet()->printout();
        if (params->separateNReads) {
            cout << "N ";
            data.divReadsSets->getNReadsSet()->printout();
        }
        ConstantLengthReadsSetInterface* readsSet =
                params->separateNReads?((ConstantLengthReadsSetInterface*) new SumOfConstantLengthReadsSets(
                    data.divReadsSets->getLqReadsSet(), data.divReadsSets->getNReadsSet())):
                        data.divReadsSets->getLqReadsSet();
        IndexesMapping* mapping = params->separateNReads?((IndexesMapping*) new SumOfMappings(
                    data.divReadsSets->getLqReadsIndexesMapping(), data.divReadsSets->getNReadsIndexesMapping())):
                        data.divReadsSets->getLqReadsIndexesMapping();
        if (!params->forceConstantParamsMode)
            params->readsExactMatchingChars += matchingCharsCorrection(data.hqPg->getPseudoGenomeLength());
        bool dumpInfoFlag = false;
        const vector<bool>& isReadMappedIntoHqPg = mapReadsIntoPg(
                data.hqPg, true, params->preserveOrderMode, readsSet,
                !params->pairFastqFile.empty() && !params->singleReadsMode, params->revComplPairFile,
                DefaultReadsMatcher::DISABLED_PREFIX_MODE,
                params->preReadsExactMatchingChars, params->readsExactMatchingChars,
                params->minCharsPerMismatch, params->preMatchingMode, params->matchingMode,
                dumpInfoFlag, pgrcOut, params->compressionLevel,
                params->extraFilesForValidation?params->pgMappedHqPrefix:"", mapping);
        uint_reads_cnt_max nBegIdx = data.divReadsSets->getLqReadsSet()->readsCount();
        data.divReadsSets->removeReadsFromLqReadsSet(isReadMappedIntoHqPg);
        if (params->separateNReads) {
            delete(mapping);
            delete(readsSet);
            data.divReadsSets->removeReadsFromNReadsSet(isReadMappedIntoHqPg, nBegIdx);
        }
    }

    void PgRCEncoder::persistMappedReadsQualityDivision() {
        data.divReadsSets->getLqReadsIndexesMapping()->saveMapping(params->mappedLqDivisionFile);
        if (params->separateNReads)
            data.divReadsSets->getNReadsIndexesMapping()->saveMapping(params->nDivisionFile);
    }

    void PgRCEncoder::compressMappedHQPgReadsList() {
        cout << "Error: unimplemented standalone compressMEMMappedPgSequences!" << endl;
        exit(EXIT_FAILURE);
    }

    void PgRCEncoder::persistHQPgSequence() {
        SeparatedPseudoGenomePersistence::writePseudoGenomeSequence(data.hqPg->getPgSequence(), params->pgMappedHqPrefix);
    }

    void PgRCEncoder::prepareForLQPgAndNPgGeneration() {
        if (!data.divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    params->srcFastqFile, params->pairFastqFile, params->revComplPairFile);
            data.divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, params->readLength, params->mappedLqDivisionFile, params->nReadsLQ,
                    params->separateNReads ? params->nDivisionFile : "", true);
            delete (allReadsIterator);
        }
    }

    void PgRCEncoder::runLQPgGeneration() {
        cout << "LQ ";
        data.divReadsSets->getLqReadsSet()->printout();
        if (numberOfThreads > 1 &&
            data.divReadsSets->getLqReadsSet()->readsCount() > PARALLEL_PG_GENERATION_READS_COUNT_THRESHOLD)
            data.lqPg = ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                    data.divReadsSets->getLqReadsSet());
        else
            data.lqPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                    data.divReadsSets->getLqReadsSet());
        data.lqPg->applyIndexesMapping(data.divReadsSets->getLqReadsIndexesMapping());
        data.divReadsSets->disposeLqReadsSet();
    }

    void PgRCEncoder::persistLQPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(data.lqPg, params->pgMappedLqPrefix);
    }

    void PgRCEncoder::compressLQPgReadsList() {
        SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(data.lqPg, &(pgrcOut),
                params->compressionLevel, params->preserveOrderMode);
        if (params->extraFilesForValidation)
            SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(data.lqPg, params->pgMappedLqPrefix, true);
    }

    void PgRCEncoder::runNPgGeneration() {
        cout << "N ";
        data.divReadsSets->getNReadsSet()->printout();
        if (numberOfThreads > 1 &&
            data.divReadsSets->getNReadsSet()->readsCount() > PARALLEL_PG_GENERATION_READS_COUNT_THRESHOLD)
            data.nPg = ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                    data.divReadsSets->getNReadsSet());
        else
            data.nPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                    data.divReadsSets->getNReadsSet());
        data.nPg->applyIndexesMapping(data.divReadsSets->getNReadsIndexesMapping());
        data.divReadsSets->disposeNReadsSet();
    }

    void PgRCEncoder::persistNPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(data.nPg, params->pgNPrefix);
    }

    void PgRCEncoder::compressNPgReadsList() {
        SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(data.nPg, &(pgrcOut),
                params->compressionLevel, params->preserveOrderMode);
        if (params->extraFilesForValidation)
            SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(data.nPg, params->pgNPrefix, true);
    }

    void PgRCEncoder::prepareForPgMatching() {
        if (!data.hqPg)
            data.hqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(params->pgHqPrefix, true);
        if (!data.lqPg)
            data.lqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(params->pgMappedLqPrefix, true);
        if (params->separateNReads && !data.nPg)
            data.nPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(params->pgNPrefix, true);
    }


    void PgRCEncoder::compressMEMMappedPgSequences() {
        cout << "Error: unimplemented standalone compressMEMMappedPgSequences!" << endl;
        exit(EXIT_FAILURE);
    }

    void PgRCEncoder::generateReport() {
        string outputfile = "pgrc_res.txt";
        bool hasHeader = (bool) std::ifstream(outputfile);
        fstream fout(outputfile, ios::out | ios::binary | ios::app);
        if (!hasHeader)
            fout << "srcFastq\tpairFastq\trcPairFile\tpgPrefix\tq[%o]\tg[%]\tm\tM\tp\tt\tsize[B]\ttotal[s]\tdiv[s]\tPgDiv[s]\tgood[s]\treadsMatch[s]\tbad&N[s]\torder[s]\tpgSeq-s[s]" << endl;

        fout << params->srcFastqFile << "\t" << params->pairFastqFile << "\t"
            << (params->revComplPairFile?"yes":"no") << "\t"
            << params->pgRCFileName << "\t" << toString(params->error_limit_in_promils) << (params->singleReadsMode?"s":"") << "\t"
            << params->gen_quality_str << "\t";
        if (params->preReadsExactMatchingChars > 0)
            fout << (char) tolower(params->preMatchingMode)
                << ((toupper(params->preMatchingMode) == params->preMatchingMode)?string("s"):string(""))
                << (int) params->preReadsExactMatchingChars;
        fout << (char) tolower(params->matchingMode)
            << ((toupper(params->matchingMode) == params->matchingMode)?string("s"):string(""))
            << (int) params->readsExactMatchingChars << "\t" << (int) params->minCharsPerMismatch << "\t";
        fout << params->targetPgMatchLength << "\t" << numberOfThreads << "\t";
        fout << pgRCSize << "\t";
        fout << getTimeInSec(chrono::steady_clock::now(), start_t, 2) << "\t";
        fout << getTimeInSec(div_t, start_t, 2) << "\t";
        fout << getTimeInSec(pgDiv_t, div_t, 2) << "\t";
        fout << getTimeInSec(good_t, pgDiv_t, 2) << "\t";
        fout << getTimeInSec(match_t, good_t, 2) << "\t";
        fout << getTimeInSec(bad_t, match_t, 2) << "\t";
        fout << getTimeInSec(order_t, bad_t, 2) << "\t";
        fout << getTimeInSec(chrono::steady_clock::now(), order_t, 2) << endl;
    }

    void PgRCEncoder::finalizeCompression() {
        pgRCSize = pgrcOut.tellp();
        cout << endl << "Created PgRC of size " << pgRCSize << " bytes in "
             << toString((double) time_millis(start_t ) / 1000, 2) << " s." << endl;
        pgrcOut.close();
        string pgRCTempFileName = params->pgRCFileName + TEMPORARY_FILE_SUFFIX;
        if (std::ifstream(params->pgRCFileName))
            remove(params->pgRCFileName.c_str());
        if (rename(pgRCTempFileName.c_str(), params->pgRCFileName.c_str()) != 0) {
            fprintf(stderr, "Error preparing output file: %s\n", params->pgRCFileName.c_str());
            exit(EXIT_FAILURE);
        };
    }

    uint_read_len_max probeReadsLength(const string &srcFastqFile) {
        ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt = ReadsSetPersistence::createManagedReadsIterator(
                srcFastqFile);
        readsIt->moveNext();
        uint_read_len_max readsLength = readsIt->getReadLength();
        delete(readsIt);
        return readsLength;
    }


}