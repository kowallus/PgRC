#include <sys/stat.h>
#include "PgRCManager.h"

#include "matching/ReadsMatchers.h"
#include "matching/SimplePgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    static const char *const BAD_INFIX = "bad";
    static const char *const GOOD_INFIX = "good";
    static const char *const N_INFIX = "N";
    static const char *const DIVISION_EXTENSION = ".div";

    static const char *const TEMPORARY_FILE_SUFFIX = ".temp";
    static const char *const PGRC_HEADER = "PgRC";

    uint_read_len_max probeReadsLength(const string &srcFastqFile);

    time_t getTimeInSec(chrono::steady_clock::time_point end_t, chrono::steady_clock::time_point begin_t) {
        chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(end_t - begin_t);
        return (double)time_span.count() / 1000000.0;
    }

    string getTimeInSec(chrono::steady_clock::time_point end_t, chrono::steady_clock::time_point begin_t, int decimalPlaces) {
        chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(end_t - begin_t);
        return toString((double)time_span.count() / 1000000.0, decimalPlaces); }

    void PgRCManager::initCompressionParameters() {
        setPreMatchingMode('c');
        switch (compressionLevel) {
            case PGRC_CODER_LEVEL_FAST:
                setQualityBasedDivisionErrorLimitInPromils(1);
                setPgGeneratorBasedDivisionOverlapThreshold_str("0");
                setMatchingMode('C');
                setPreReadsExactMatchingChars(0);
                setReadSeedLength(54);
                setMinCharsPerMismatch(2);
                setMinimalPgReverseComplementedRepeatLength(25);
                break;
            case PGRC_CODER_LEVEL_NORMAL:
                setQualityBasedDivisionErrorLimitInPromils(1000);
                setPgGeneratorBasedDivisionOverlapThreshold_str("65");
                setPreReadsExactMatchingChars(0);
                setMatchingMode('c');
                setReadSeedLength(38);
                setMinCharsPerMismatch(6);
                setMinimalPgReverseComplementedRepeatLength(50);
                break;
            case PGRC_CODER_LEVEL_MAX:
                setQualityBasedDivisionErrorLimitInPromils(1000);
                setPgGeneratorBasedDivisionOverlapThreshold_str("65");
                setPreReadsExactMatchingChars(64);
                setMatchingMode('c');
                setReadSeedLength(30);
                setMinCharsPerMismatch(6);
                setMinimalPgReverseComplementedRepeatLength(50);
                break;
            default:
                fprintf(stderr, "Error: unknown compression level: %d.", compressionLevel);
                exit(EXIT_FAILURE);
        }
    }

    void PgRCManager::prepareChainData() {
        initCompressionParameters();
        qualityDivision = error_limit_in_promils < 1000;
        generatorDivision = gen_quality_coef > 0;
        readLength = probeReadsLength(srcFastqFile);
        if (pairFastqFile.empty() && !preserveOrderMode)
            singleReadsMode = true;

        if (std::ifstream(pgRCFileName))
            fprintf(stderr, "Warning: file %s already exists\n", pgRCFileName.data());
        pgrcOut = fstream(pgRCFileName + TEMPORARY_FILE_SUFFIX, ios::out | ios::binary | ios::trunc);
        pgrcOut.write(PGRC_HEADER, strlen(PGRC_HEADER));
        const char pgrc_mode = singleReadsMode ? PGRC_SE_MODE :
                               (preserveOrderMode?(pairFastqFile.empty()?PGRC_ORD_SE_MODE:PGRC_ORD_PE_MODE):
                     (ignorePairOrderInformation?PGRC_MIN_PE_MODE:PGRC_PE_MODE));
        pgrcOut.put(pgrc_mode);
        pgrcOut.put(separateNReads);
        revComplPairFile = false;
        if (pgrc_mode == PGRC_PE_MODE || pgrc_mode == PGRC_ORD_PE_MODE) {
            if (!disableRevComplPairFileMode)
                revComplPairFile = true;
            pgrcOut.put(revComplPairFile);
        }

        string tmpDirectoryName = pgRCFileName;
        if (qualityDivision)
            tmpDirectoryName = tmpDirectoryName + "_q" + toString(error_limit_in_promils);
        tmpDirectoryName = tmpDirectoryName + (nReadsLQ?"_n":"") + (separateNReads?"_N":"") + "_g" + gen_quality_str;
        bool enablePreReadsMatching = preReadsExactMatchingChars > 0;
        tmpDirectoryName = tmpDirectoryName +
                (enablePreReadsMatching?("_l" + (((char) tolower(preMatchingMode))
                                                 + ((toupper(preMatchingMode) == preMatchingMode)?string("s"):string(""))
                                                 + toString(preReadsExactMatchingChars))):"")
                + "_m" + (((char) tolower(matchingMode))
                + (toupper(matchingMode) == matchingMode?string("s"):string(""))
                + toString(readsExactMatchingChars))
                + "_M" + toString(minCharsPerMismatch)
                + "_p" + toString(targetPgMatchLength);
        if (extraFilesForValidation && skipStages == 0) {
            mode_t mode = 0777;
            int nError = mkdir(tmpDirectoryName.data(), mode);
            if (nError != 0) {
                srand(time(NULL));
                tmpDirectoryName = tmpDirectoryName + "_" + toString(rand() % 100000);
                nError = mkdir(tmpDirectoryName.data(), mode);
                if (nError != 0) {
                    fprintf(stderr, "Error creating folder %s\n", tmpDirectoryName.data());
                    exit(EXIT_FAILURE);
                }
            }
        }
        pgrcOut.write(tmpDirectoryName.data(), tmpDirectoryName.length());
        pgrcOut << endl;

        lqDivisionFile = tmpDirectoryName + "/" + BAD_INFIX + DIVISION_EXTENSION;
        nDivisionFile = tmpDirectoryName + "/" + N_INFIX + DIVISION_EXTENSION;
        pgHqPrefix = tmpDirectoryName + "/" + GOOD_INFIX;

        pgFilesPrefixesWithM = tmpDirectoryName + "/";
        pgMappedHqPrefix = pgFilesPrefixesWithM + GOOD_INFIX;
        pgMappedLqPrefix = pgFilesPrefixesWithM + BAD_INFIX;
        pgSeqFinalHqPrefix = pgFilesPrefixesWithM + GOOD_INFIX;
        pgSeqFinalLqPrefix = pgFilesPrefixesWithM + BAD_INFIX;
        pgNPrefix = pgFilesPrefixesWithM + N_INFIX;
        mappedLqDivisionFile = pgFilesPrefixesWithM + BAD_INFIX + DIVISION_EXTENSION;
    }

    void PgRCManager::executePgRCChain() {
        start_t = chrono::steady_clock::now();
        prepareChainData();
        stageCount = 0;
        if (skipStages < ++stageCount && qualityDivision) {
            runQualityBasedDivision();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistReadsQualityDivision();
                disposeChainData();
            }
        }
        div_t = chrono::steady_clock::now();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForPgGeneratorBaseReadsDivision();
            if (generatorDivision)
                runPgGeneratorBasedReadsDivision();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistReadsQualityDivision();
                disposeChainData();
            }

        }
        pgDiv_t = chrono::steady_clock::now();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForHqPgGeneration();
            runHQPgGeneration();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistHQPg();
                persistReadsQualityDivision();
                disposeChainData();
            }
        }
        good_t = chrono::steady_clock::now();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForMappingLQReadsOnHQPg();
            runMappingLQReadsOnHQPg();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistHQPgSequence();
                persistMappedReadsQualityDivision();
                disposeChainData();
            } else {
                //// Already done during runMappingLQReadsOnHQPg()
//                compressMappedHQPgReadsList();
                if (!singleReadsMode) {
                    if (preserveOrderMode)
                        orgIdx2PgPos = std::move(hqPg->getReadsList()->pos);
                    else
                        rlIdxOrder = std::move(hqPg->getReadsList()->orgIdx);
                }
                hqPg->disposeReadsList();
            }
        }
        match_t = chrono::steady_clock::now();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForLQPgAndNPgGeneration();
            runLQPgGeneration();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistLQPg();
                if (hqPg)
                    persistHQPgSequence();
            } else {
                compressLQPgReadsList();
            }
            if (!singleReadsMode) {
                if (preserveOrderMode) {
                    ExtendedReadsListWithConstantAccessOption *const pgRl = lqPg->getReadsList();
                    uint_pg_len_max pos = hqPg->getPseudoGenomeLength();
                    for (uint_reads_cnt_std i = 0; i < pgRl->readsCount; i++) {
                        pos += pgRl->off[i];
                        orgIdx2PgPos[pgRl->orgIdx[i]] = pos;
                    }
                } else
                    rlIdxOrder.insert(rlIdxOrder.end(), lqPg->getReadsList()->orgIdx.begin(), lqPg->getReadsList()->orgIdx.end());
            }
            lqPg->disposeReadsList();
            if (separateNReads) {
                runNPgGeneration();
                if (disableInMemoryMode || endAtStage == stageCount) {
                    persistNPg();
                } else {
                    compressNPgReadsList();
                }
                if (!singleReadsMode) {
                    if (preserveOrderMode) {
                        ExtendedReadsListWithConstantAccessOption *const pgRl = nPg->getReadsList();
                        uint_pg_len_max pos = hqPg->getPseudoGenomeLength() + lqPg->getPseudoGenomeLength();
                        for (uint_reads_cnt_std i = 0; i < pgRl->readsCount; i++) {
                            pos += pgRl->off[i];
                            orgIdx2PgPos[pgRl->orgIdx[i]] = pos;
                        }
                    } else
                        rlIdxOrder.insert(rlIdxOrder.end(), nPg->getReadsList()->orgIdx.begin(),
                                   nPg->getReadsList()->orgIdx.end());
                }
                nPg->disposeReadsList();
            }
            delete(divReadsSets);
            divReadsSets = 0;
        }
        bad_t = chrono::steady_clock::now();
        if (!singleReadsMode && skipStages < ++stageCount && endAtStage >= stageCount) {
            if (preserveOrderMode) {
                const uint_pg_len_max joinedPgLength =
                        hqPg->getPseudoGenomeLength() + lqPg->getPseudoGenomeLength() + nPg->getPseudoGenomeLength();
                bool isJoinedPgLengthStd = joinedPgLength <= UINT32_MAX;
                if (isJoinedPgLengthStd)
                    SeparatedPseudoGenomePersistence::compressReadsPgPositions<uint_pg_len_std>(pgrcOut, orgIdx2PgPos,
                            joinedPgLength, compressionLevel, pairFastqFile.empty());
                else
                    SeparatedPseudoGenomePersistence::compressReadsPgPositions<uint_pg_len_max>(pgrcOut, orgIdx2PgPos,
                            joinedPgLength, compressionLevel, pairFastqFile.empty());
                orgIdx2PgPos.clear();
            } else {
                SeparatedPseudoGenomePersistence::compressReadsOrder(pgrcOut, rlIdxOrder, compressionLevel, preserveOrderMode,
                                                                     ignorePairOrderInformation, pairFastqFile.empty());
                rlIdxOrder.clear();
            }
            if (extraFilesForValidation) {
                if (separateNReads)
                    SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedHqPrefix, pgMappedLqPrefix, pgNPrefix});
                else
                    SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedHqPrefix, pgMappedLqPrefix});
            }
        }
        order_t = chrono::steady_clock::now();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForPgMatching();
            string emptySequence;
            SimplePgMatcher::matchPgsInPg(hqPg->getPgSequence(), lqPg->getPgSequence(),
                    separateNReads?nPg->getPgSequence():emptySequence, separateNReads, pgrcOut, compressionLevel,
                                          extraFilesForValidation?pgSeqFinalHqPrefix:"",
                                          extraFilesForValidation?pgSeqFinalLqPrefix:"",
                                          extraFilesForValidation?pgNPrefix:"", targetPgMatchLength);
        }
        finalizeCompression();
        disposeChainData();
#ifdef DEVELOPER_BUILD
        generateReport();
#endif
    }

    void PgRCManager::runQualityBasedDivision() {
        ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                srcFastqFile, pairFastqFile, revComplPairFile);
        divReadsSets =
                DividedPCLReadsSets::getQualityDivisionBasedReadsSets(allReadsIterator, readLength, error_limit_in_promils / 1000.0,
                        separateNReads, nReadsLQ);
        delete (allReadsIterator);
    }

    void PgRCManager::persistReadsQualityDivision() {
        divReadsSets->getLqReadsIndexesMapping()->saveMapping(lqDivisionFile);
        if (separateNReads)
            divReadsSets->getNReadsIndexesMapping()->saveMapping(nDivisionFile);
    }

    void PgRCManager::prepareForPgGeneratorBaseReadsDivision() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    srcFastqFile, pairFastqFile, revComplPairFile);
            if (qualityDivision) {
                divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                        allReadsIterator, readLength, lqDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "");
            } else {
                divReadsSets = DividedPCLReadsSets::getSimpleDividedPCLReadsSets(allReadsIterator, readLength,
                                                                                     separateNReads, nReadsLQ);
            }
            delete (allReadsIterator);
        }
    }

    void PgRCManager::runPgGeneratorBasedReadsDivision() {
        cout << "HQ ";
        divReadsSets->getHqReadsSet()->printout();
        const vector<bool>& isReadHqInHqReadsSet = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getHQReads(
                divReadsSets->getHqReadsSet(), gen_quality_coef);
        divReadsSets->moveLqReadsFromHqReadsSetsToLqReadsSets(isReadHqInHqReadsSet);
    }

    void PgRCManager::prepareForHqPgGeneration() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    srcFastqFile, pairFastqFile, revComplPairFile);
            divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, readLength, lqDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "");
            delete (allReadsIterator);
        }
    }

    void PgRCManager::runHQPgGeneration() {
        cout << "HQ ";
        divReadsSets->getHqReadsSet()->printout();
        hqPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                divReadsSets->getHqReadsSet());
        divReadsSets->disposeHqReadsSet();
        IndexesMapping* hq2IndexesMapping = divReadsSets->generateHqReadsIndexesMapping();
        hqPg->applyIndexesMapping(hq2IndexesMapping);
        delete(hq2IndexesMapping);
    }

    void PgRCManager::persistHQPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(hqPg, pgHqPrefix);
    }

    void PgRCManager::prepareForMappingLQReadsOnHQPg() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    srcFastqFile, pairFastqFile, revComplPairFile);
            divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, readLength, lqDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "", true);
            delete (allReadsIterator);
        }
        if (!hqPg)
            hqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgHqPrefix);
    }

    void PgRCManager::runMappingLQReadsOnHQPg() {
        cout << "LQ ";
        divReadsSets->getLqReadsSet()->printout();
        if (separateNReads) {
            cout << "N ";
            divReadsSets->getNReadsSet()->printout();
        }
        ConstantLengthReadsSetInterface* readsSet =
                separateNReads?((ConstantLengthReadsSetInterface*) new SumOfConstantLengthReadsSets(
                        divReadsSets->getLqReadsSet(), divReadsSets->getNReadsSet())): divReadsSets->getLqReadsSet();
        IndexesMapping* mapping = separateNReads?((IndexesMapping*) new SumOfMappings(
                divReadsSets->getLqReadsIndexesMapping(), divReadsSets->getNReadsIndexesMapping())):divReadsSets->getLqReadsIndexesMapping();
        if (!forceConstantParamsMode)
            readsExactMatchingChars += matchingCharsCorrection(hqPg->getPseudoGenomeLength());
        const vector<bool>& isReadMappedIntoHqPg = mapReadsIntoPg(
                hqPg, true, preserveOrderMode, readsSet, !pairFastqFile.empty() && !singleReadsMode, revComplPairFile,
                DefaultReadsMatcher::DISABLED_PREFIX_MODE,
                preReadsExactMatchingChars, readsExactMatchingChars,
                minCharsPerMismatch, preMatchingMode, matchingMode,
                false, pgrcOut, compressionLevel, extraFilesForValidation?pgMappedHqPrefix:"", mapping);
        uint_reads_cnt_max nBegIdx = divReadsSets->getLqReadsSet()->readsCount();
        divReadsSets->removeReadsFromLqReadsSet(isReadMappedIntoHqPg);
        if (separateNReads) {
            delete(mapping);
            delete(readsSet);
            divReadsSets->removeReadsFromNReadsSet(isReadMappedIntoHqPg, nBegIdx);
        }
    }

    void PgRCManager::persistMappedReadsQualityDivision() {
        divReadsSets->getLqReadsIndexesMapping()->saveMapping(mappedLqDivisionFile);
        if (separateNReads)
            divReadsSets->getNReadsIndexesMapping()->saveMapping(nDivisionFile);
    }

    void PgRCManager::compressMappedHQPgReadsList() {
        cout << "Error: unimplemented standalone compressMEMMappedPgSequences!" << endl;
        exit(EXIT_FAILURE);
    }

    void PgRCManager::persistHQPgSequence() {
        SeparatedPseudoGenomePersistence::writePseudoGenomeSequence(hqPg->getPgSequence(), pgMappedHqPrefix);
    }

    void PgRCManager::prepareForLQPgAndNPgGeneration() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                    srcFastqFile, pairFastqFile, revComplPairFile);
            divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, readLength, mappedLqDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "", true);
            delete (allReadsIterator);
        }
    }

    void PgRCManager::runLQPgGeneration() {
        cout << "LQ ";
        divReadsSets->getLqReadsSet()->printout();
        lqPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                divReadsSets->getLqReadsSet());
        lqPg->applyIndexesMapping(divReadsSets->getLqReadsIndexesMapping());
        divReadsSets->disposeLqReadsSet();
    }

    void PgRCManager::persistLQPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(lqPg, pgMappedLqPrefix);
    }

    void PgRCManager::compressLQPgReadsList() {
        SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(lqPg, &pgrcOut, compressionLevel, preserveOrderMode);
        if (extraFilesForValidation)
            SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(lqPg, pgMappedLqPrefix, true);
    }

    void PgRCManager::runNPgGeneration() {
        cout << "N ";
        divReadsSets->getNReadsSet()->printout();
        nPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                divReadsSets->getNReadsSet());
        nPg->applyIndexesMapping(divReadsSets->getNReadsIndexesMapping());
        divReadsSets->disposeNReadsSet();
    }

    void PgRCManager::persistNPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(nPg, pgNPrefix);
    }

    void PgRCManager::compressNPgReadsList() {
        SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(nPg, &pgrcOut, compressionLevel, preserveOrderMode);
        if (extraFilesForValidation)
            SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(nPg, pgNPrefix, true);
    }

    void PgRCManager::prepareForPgMatching() {
        if (!hqPg)
            hqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgHqPrefix, true);
        if (!lqPg)
            lqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgMappedLqPrefix, true);
        if (separateNReads && !nPg)
            nPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgNPrefix, true);
    }


    void PgRCManager::compressMEMMappedPgSequences() {
        cout << "Error: unimplemented standalone compressMEMMappedPgSequences!" << endl;
        exit(EXIT_FAILURE);
    }

    void PgRCManager::generateReport() {
        string outputfile = "pgrc_res.txt";
        bool hasHeader = (bool) std::ifstream(outputfile);
        fstream fout(outputfile, ios::out | ios::binary | ios::app);
        if (!hasHeader)
            fout << "srcFastq\tpairFastq\trcPairFile\tpgPrefix\tq[%o]\tg[%]\tm\tM\tp\tsize[B]\ttotal[s]\tdiv[s]\tPgDiv[s]\tgood[s]\treadsMatch[s]\tbad&N[s]\torder[s]\tpgSeq-s[s]" << endl;

        fout << srcFastqFile << "\t" << pairFastqFile << "\t" << (revComplPairFile?"yes":"no") << "\t"
             << pgRCFileName << "\t" << toString(error_limit_in_promils) << "\t" << gen_quality_str << "\t";

        if (preReadsExactMatchingChars > 0)
            fout << (char) tolower(preMatchingMode) << ((toupper(preMatchingMode) == preMatchingMode)?string("s"):string("")) << (int) preReadsExactMatchingChars;
        fout << (char) tolower(matchingMode) << ((toupper(matchingMode) == matchingMode)?string("s"):string("")) << (int) readsExactMatchingChars << "\t" << (int) minCharsPerMismatch << "\t" << targetPgMatchLength << "\t";
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

    void PgRCManager::finalizeCompression() {
        pgRCSize = pgrcOut.tellp();
        cout << endl << "Created PgRC of size " << pgRCSize << " bytes in "
             << toString((double) time_millis(start_t ) / 1000, 2) << " s." << endl;
        pgrcOut.close();
        string pgRCTempFileName = pgRCFileName + TEMPORARY_FILE_SUFFIX;
        if (std::ifstream(pgRCFileName))
            remove(pgRCFileName.c_str());
        if (rename(pgRCTempFileName.c_str(), pgRCFileName.c_str()) != 0) {
            fprintf(stderr, "Error preparing output file: %s\n", pgRCFileName.c_str());
            exit(EXIT_FAILURE);
        };
    }

    void PgRCManager::disposeChainData() {
        if (divReadsSets) {
            delete (divReadsSets);
            divReadsSets = 0;
        }
        if (hqPg) {
            delete (hqPg);
            hqPg = 0;
        }
        if (lqPg) {
            delete (lqPg);
            lqPg = 0;
        }
        if (nPg) {
            delete (nPg);
            nPg = 0;
        }
    }

    void PgRCManager::decompressPgRC() {
        start_t = chrono::steady_clock::now();
        string tmpDirectoryPath = pgRCFileName + "/";
        ifstream pgrcIn(pgRCFileName);
        char pgrc_mode;
        if (pgrcIn) {
            for(int i = 0; i < strlen(PGRC_HEADER); i++) {
                if (PGRC_HEADER[i] != pgrcIn.get()) {
                    fprintf(stderr, "Error processing header.\n");
                    exit(EXIT_FAILURE);
                }
            }
            pgrc_mode = pgrcIn.get();
            if (pgrc_mode < PGRC_SE_MODE || pgrc_mode > PGRC_MIN_PE_MODE) {
                fprintf(stderr, "Unsupported decompression mode: %d.\n", pgrc_mode);
                exit(EXIT_FAILURE);
            }
            separateNReads = (bool) pgrcIn.get();
            if (pgrc_mode == PGRC_PE_MODE || pgrc_mode == PGRC_ORD_PE_MODE)
                revComplPairFile = (bool) pgrcIn.get();

            pgrcIn >> tmpDirectoryPath;
            tmpDirectoryPath = tmpDirectoryPath + "/";
            pgrcIn.get();
        }

        pgMappedHqPrefix = tmpDirectoryPath + GOOD_INFIX;
        pgMappedLqPrefix = tmpDirectoryPath + BAD_INFIX;
        pgSeqFinalHqPrefix = tmpDirectoryPath + GOOD_INFIX;
        pgSeqFinalLqPrefix = tmpDirectoryPath + BAD_INFIX;
        pgNPrefix = tmpDirectoryPath + N_INFIX;

        preserveOrderMode = pgrc_mode == PGRC_ORD_SE_MODE || pgrc_mode == PGRC_ORD_PE_MODE;
        ignorePairOrderInformation = pgrc_mode == PGRC_MIN_PE_MODE;
        singleReadsMode = pgrc_mode == PGRC_SE_MODE || pgrc_mode == PGRC_ORD_SE_MODE;
        if (pgrcIn)
            loadAllPgs(pgrcIn);
        else
            loadAllPgs();
        cout << "... loaded Pgs (checkpoint: " << time_millis(start_t) << " msec.)" << endl;

        uint_reads_cnt_max nonNPgReadsCount = hqPg->getReadsSetProperties()->readsCount
                                              + lqPg->getReadsSetProperties()->readsCount;
        uint_reads_cnt_max nPgReadsCount = nPg?nPg->getReadsSetProperties()->readsCount:0;
        uint_reads_cnt_max readsTotalCount = nonNPgReadsCount + nPgReadsCount;
        if (revComplPairFile)
            if (isJoinedPgLengthStd)
                applyRevComplPairFileToPgs<uint_pg_len_std>(orgIdx2StdPgPos);
            else
                applyRevComplPairFileToPgs<uint_pg_len_max>(orgIdx2PgPos);

        if (srcFastqFile.empty()) {
            if (singleReadsMode && !preserveOrderMode) {
                if (ENABLE_PARALLEL_DECOMPRESSION && dnaStreamSize() > CHUNK_SIZE_IN_BYTES)
                    writeAllReadsInSEModeParallel(pgRCFileName);
                else
                    writeAllReadsInSEMode(pgRCFileName);
            }
            else if (!preserveOrderMode)
                writeAllReadsInPEMode(pgRCFileName);
            else if (isJoinedPgLengthStd)
                writeAllReadsInORDMode<uint_pg_len_std>(pgRCFileName, orgIdx2StdPgPos);
            else
                writeAllReadsInORDMode<uint_pg_len_max>(pgRCFileName, orgIdx2PgPos);
            cout << "Decompressed ";
        } else {
            preparePgsForValidation();
            validateAllPgs();
            validatePgsOrder();
            cout << "Validated ";
        }

        cout << readsTotalCount << " reads in " << time_millis(start_t) << " msec." << endl;

        disposeChainData();
    }

    void PgRCManager::pushOutToQueue(string &out) {
        std::lock_guard<std::mutex> _(this->mut);
        out_queue.push(std::move(out));
        data_cond.notify_one();
        out.resize(0);
        out.reserve(CHUNK_SIZE_IN_BYTES);
    }

    void PgRCManager::finishWritingParallel() {
        std::lock_guard<std::mutex> _(this->mut);
        out_queue.push(std::move(""));
        data_cond.notify_one();
    }

    void PgRCManager::writeAllReadsInSEModeParallel(const string &outPrefix) {
        cout << "... parallel mode" << endl;
        hqPg->getReadsList()->enableConstantAccess(true);
        lqPg->getReadsList()->enableConstantAccess(true);
        if (nPg) nPg->getReadsList()->enableConstantAccess(true);
        std::thread writing(&PgRCManager::writeFromQueue, this, outPrefix);
        string res, read;
        read.resize(readLength);
        uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
        res.reserve(CHUNK_SIZE_IN_BYTES);
        for(uint_reads_cnt_max i = 0; i < hqReadsCount; i++) {
            if (res.size() > res_size_guard) {
                pushOutToQueue(res);
            }
            hqPg->getRead_Unsafe(i, (char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        for(uint_reads_cnt_max i = 0; i < lqReadsCount; i++) {
            if (res.size() > res_size_guard) {
                pushOutToQueue(res);
            }
            lqPg->getRead_RawSequence(i, (char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        for(uint_reads_cnt_max i = 0; i < nPgReadsCount; i++) {
            if (res.size() > res_size_guard) {
                pushOutToQueue(res);
            }
            nPg->getRead_RawSequence(i, (char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        pushOutToQueue(res);
        cout << "... finished loading queue (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        finishWritingParallel();
        writing.join();
    }

    void PgRCManager::writeFromQueue(const string &outPrefix) {
        fstream fout(outPrefix + "_out", ios_base::out | ios_base::binary | std::ios::trunc);
        string out;
        do {
            std::unique_lock<std::mutex> lk(mut);
            data_cond.wait(
                    lk,[this]{return !out_queue.empty();});
            out = out_queue.front();
            out_queue.pop();
            lk.unlock();
            fout << out;
        } while (!out.empty());
        fout.close();
    }

    void PgRCManager::writeAllReadsInSEMode(const string &outPrefix) const {
        fstream fout(outPrefix + "_out", ios_base::out | ios_base::binary | std::ios::trunc);
        string res, read;
        read.resize(readLength);
        uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
        uint64_t totalSize = dnaStreamSize();
        res.reserve(totalSize < res_size_guard?totalSize:res_size_guard + (readLength + 1));
        for(uint_reads_cnt_max i = 0; i < hqReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            hqPg->getNextRead_Unsafe((char *) read.data());
            res.append(read);
            res.push_back('\n');
        }
        for(uint_reads_cnt_max i = 0; i < lqReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            lqPg->getNextRead_RawSequence((char*) read.data());
            res.append(read);
            res.push_back('\n');
        }
        for(uint_reads_cnt_max i = 0; i < nPgReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            nPg->getNextRead_RawSequence((char*) read.data());
            res.append(read);
            res.push_back('\n');
        }
        fout << res;
        fout.close();
    }

    void PgRCManager::writeAllReadsInPEMode(const string &outPrefix) const {
        hqPg->getReadsList()->enableConstantAccess(true);
        lqPg->getReadsList()->enableConstantAccess(true);
        if (nPg) nPg->getReadsList()->enableConstantAccess(true);
        const uint8_t PE_PARTS_COUNT = 2;

        for(uint8_t p = 0; p < PE_PARTS_COUNT; p++) {
            fstream fout(outPrefix + "_out" + "_" + toString(p + 1),
                         ios_base::out | ios_base::binary | std::ios::trunc);
            string res, read;
            read.resize(readLength);
            char* readPtr = (char *) read.data();
            uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
            uint64_t totalSize = dnaStreamSize();
            res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (readLength + 1));
            uint_reads_cnt_max i = 0;
            for (i = p; i < readsTotalCount; i += PE_PARTS_COUNT) {
                if (res.size() > res_size_guard) {
                    fout << res;
                    res.resize(0);
                }
                uint_reads_cnt_std idx = rlIdxOrder[i];
                if (idx < hqReadsCount)
                    hqPg->getRead(idx, readPtr);
                else {
                    if (idx < nonNPgReadsCount)
                        lqPg->getRead_RawSequence(idx - hqReadsCount, readPtr);
                    else
                        nPg->getRead_RawSequence(idx - nonNPgReadsCount, readPtr);
                    if (p)
                        PgSAHelpers::reverseComplementInPlace(readPtr, readLength);
                }
                res.append(read);
                res.push_back('\n');
            }
            fout << res;
            fout.close();
        }
    }

    template <typename uint_pg_len>
    void PgRCManager::writeAllReadsInORDMode(const string &outPrefix, vector<uint_pg_len> &orgIdx2PgPos) const {
        uint8_t parts = singleReadsMode?1:2;

        for(uint8_t p = 0; p < parts; p++) {
            fstream fout(outPrefix + "_out" + (singleReadsMode?"":("_" + toString(p + 1))),
                    ios_base::out | ios_base::binary | std::ios::trunc);
            string res, read;
            read.resize(readLength);
            char* readPtr = (char *) read.data();
            uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
            uint64_t totalSize = dnaStreamSize();
            res.reserve(totalSize < res_size_guard ? totalSize : res_size_guard + (readLength + 1));
            uint_reads_cnt_max i = 0;
            uint_reads_cnt_std endI = (readsTotalCount / parts) * (p + 1);
            for (i = (readsTotalCount / parts) * p; i < endI; i++) {
                if (res.size() > res_size_guard) {
                    fout << res;
                    res.resize(0);
                }
                uint_pg_len pos = orgIdx2PgPos[i];
                if (pos < hqPgLen)
                    hqPg->getNextRead_Unsafe(readPtr, pos);
                else {
                    if (pos < nonNPgLen)
                        lqPg->getRawSequenceOfReadLength(readPtr, pos - hqPgLen);
                    else
                        nPg->getRawSequenceOfReadLength(readPtr, pos - nonNPgLen);
                    if (p)
                        PgSAHelpers::reverseComplementInPlace(readPtr, readLength);
                }
                res.append(read);
                res.push_back('\n');
            }
            fout << res;
            fout.close();
        }
    }
    template void PgRCManager::writeAllReadsInORDMode<uint_pg_len_std>(const string &outPrefix, vector<uint_pg_len_std> &orgIdx2PgPos) const;
    template void PgRCManager::writeAllReadsInORDMode<uint_pg_len_max>(const string &outPrefix, vector<uint_pg_len_max> &orgIdx2PgPos) const;

    uint_reads_cnt_max PgRCManager::dnaStreamSize() const {
        return (hqPg->getReadsSetProperties()->readsCount + lqPg->getReadsSetProperties()->readsCount) *
               (hqPg->getReadsSetProperties()->maxReadLength + 1);
    }

    void PgRCManager::preparePgsForValidation() const {
        if (preserveOrderMode) {
            lqPg->getReadsList()->orgIdx.clear();
            nPg->getReadsList()->orgIdx.clear();
            uint8_t parts = singleReadsMode?1:2;
            for (uint_reads_cnt_std i = 0; i < readsTotalCount; i++) {
                uint_pg_len_max pos = isJoinedPgLengthStd ? orgIdx2StdPgPos[i] : orgIdx2PgPos[i];
                uint_reads_cnt_std orgIdx = i < readsTotalCount / parts ? i * parts : (i - readsTotalCount / parts) * parts + 1;
                if (pos < hqPgLen)
                    hqPg->getReadsList()->pos.push_back(pos);
                else if (pos < nonNPgLen) {
                    lqPg->getReadsList()->pos.push_back(pos - hqPgLen);
                    lqPg->getReadsList()->orgIdx.push_back(orgIdx);
                } else {
                    nPg->getReadsList()->pos.push_back(pos - nonNPgLen);
                    nPg->getReadsList()->orgIdx.push_back(orgIdx);
                }
            }
            hqPg->getReadsList()->enableConstantAccess(true);
        } else {
            hqPg->getReadsList()->enableConstantAccess(true);
            lqPg->getReadsList()->enableConstantAccess(true);
            if (nPg) nPg->getReadsList()->enableConstantAccess(true);
        }
    }

    void PgRCManager::validateAllPgs() {
        vector<uint_reads_cnt_max> orgIdx2rlIdx = getAllPgsOrgIdxs2RlIdx();
        ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                srcFastqFile, pairFastqFile);

        vector<bool> validated(readsTotalCount, false);
        uint_reads_cnt_max notValidatedCount = 0;
        uint_reads_cnt_max errorsCount = 0;
        for(uint_reads_cnt_max i = 0; i < readsTotalCount; i++) {
            if (!allReadsIterator->moveNext()) {
                cout << "The number of compressed reads is too big (" << (readsTotalCount - i) << " reads more)." << endl;
                break;
            }
            string read;
            uint_reads_cnt_max idx = orgIdx2rlIdx[i];
            //uint_reads_cnt_max idx = (i % 2) * (readsTotalCount / 2) + (i / 2);
            if (validated[idx]) {
                notValidatedCount++;
                continue;
            }
            validated[idx] = true;
            if (idx < hqReadsCount)
                read = hqPg->getRead(idx);
            else {
                if (idx < nonNPgReadsCount)
                    read = lqPg->getRead(idx - hqReadsCount);
                else
                    read = nPg->getRead(idx - nonNPgReadsCount);
                if (!singleReadsMode && !ignorePairOrderInformation && (i % 2 == 1))
                    PgSAHelpers::reverseComplementInPlace(read);
            }
            if (read != allReadsIterator->getRead())
                errorsCount++;
        }
        uint_reads_cnt_max missingCount = 0;
        while (allReadsIterator->moveNext())
            missingCount++;
        if (missingCount)
            cout << "The number of compressed reads is too small (" << (missingCount) << " reads missing)." << endl;
        if (notValidatedCount)
            cout << notValidatedCount << " of compressed reads could not been properly validated." << endl;
        if (errorsCount)
            cout << "Found " << errorsCount << " errors in compressed reads." << endl;
        if (!missingCount && !notValidatedCount && !errorsCount)
            cout << "Validation successful!" << endl;

        delete (allReadsIterator);

    }

    const vector<uint_reads_cnt_max> PgRCManager::getAllPgsOrgIdxs2RlIdx() const {
        vector<uint_reads_cnt_max> orgIdx2rlIdx;
        orgIdx2rlIdx.resize(readsTotalCount);
        for(uint_reads_cnt_max i = 0; i < hqPg->getReadsSetProperties()->readsCount; i++)
            orgIdx2rlIdx[hqPg->getReadsList()->orgIdx[i]] = i;
        for(uint_reads_cnt_max i = 0; i < lqPg->getReadsSetProperties()->readsCount; i++)
            orgIdx2rlIdx[lqPg->getReadsList()->orgIdx[i]] = hqPg->getReadsSetProperties()->readsCount + i;
        for(uint_reads_cnt_max i = 0; i < nPgReadsCount; i++)
            orgIdx2rlIdx[nPg->getReadsList()->orgIdx[i]] = nonNPgReadsCount + i;
        return orgIdx2rlIdx;
    }

    const uint_reads_cnt_max PgRCManager::getAllPgsOrgIdx(uint_reads_cnt_max idx) const {
        if (idx < hqReadsCount)
            return hqPg->getReadsList()->orgIdx[idx];
        else if (idx < nonNPgReadsCount)
            return lqPg->getReadsList()->orgIdx[idx - hqReadsCount];
        else
            return nPg->getReadsList()->orgIdx[idx - nonNPgReadsCount];
    }

    void PgRCManager::validatePgsOrder() {
        if (!preserveOrderMode && singleReadsMode)
            return;

        vector<bool> validated(readsTotalCount, false);
        uint_reads_cnt_max notValidatedCount = 0;
        uint_reads_cnt_max errorsCount = 0;
        if (preserveOrderMode) {
            rlIdxOrder = getAllPgsOrgIdxs2RlIdx();
            for(uint_reads_cnt_std i = 0; i < readsTotalCount; i++) {
                uint_reads_cnt_std rlIdx = rlIdxOrder[i];
                if (validated[rlIdx]) {
                    notValidatedCount++;
                    continue;
                }
                if (i != getAllPgsOrgIdx(rlIdx))
                    errorsCount++;
            }
        } else {
            for(uint_reads_cnt_std p = 0; p < readsTotalCount / 2; p++) {
                uint_reads_cnt_std rlIdx = rlIdxOrder[p * 2];
                uint_reads_cnt_std rlPairIdx = rlIdxOrder[p * 2 + 1];
                if (validated[rlIdx]) notValidatedCount++;
                if (validated[rlPairIdx]) notValidatedCount++;
                if (!validated[rlIdx] && !validated[rlPairIdx]) {
                    validated[rlIdx] = true;
                    validated[rlPairIdx] = true;
                    uint_reads_cnt_std orgIdx = getAllPgsOrgIdx(rlIdx);
                    uint_reads_cnt_std orgPairIdx = getAllPgsOrgIdx(rlPairIdx);
                    uint_reads_cnt_std smallerIdx = orgIdx < orgPairIdx ? orgIdx : orgPairIdx;
                    uint_reads_cnt_std largerIdx = orgIdx >= orgPairIdx ? orgIdx : orgPairIdx;
                    if (largerIdx - smallerIdx != 1 || smallerIdx % 2)
                        errorsCount++;
                    else if (!ignorePairOrderInformation && smallerIdx != orgIdx)
                        errorsCount++;
                }
            }
        }
        if (notValidatedCount)
            cout << "Order of " << notValidatedCount << " compressed reads could not been properly validated." << endl;
        if (errorsCount)
            cout << "Found " << errorsCount << " errors in compressed reads order." << endl;
        if (!notValidatedCount && !errorsCount)
            cout << "Order validation successful!" << endl;
    }

    template<typename uint_pg_len>
    void PgRCManager::applyRevComplPairFileToPgs(vector<uint_pg_len> &orgIdx2PgPos) {
        if (preserveOrderMode) {
            uint_reads_cnt_std hqRlIdx = 0;
            const uint_reads_cnt_max pairsCount = readsTotalCount / 2;
            for (uint_reads_cnt_max i = 0; i < pairsCount; i++) {
                    uint_pg_len pgPos = orgIdx2PgPos[i];
                    if (pgPos < hqPgLen)
                        hqRlIdx++;
            }
            for (uint_reads_cnt_max i = pairsCount; i < readsTotalCount; i++) {
                uint_pg_len pgPos = orgIdx2PgPos[i];
                if (pgPos < hqPgLen) {
                    hqPg->getReadsList()->revComp[hqRlIdx] = !hqPg->getReadsList()->revComp[hqRlIdx];
                    hqRlIdx++;
                }
            }
        } else {
            for (uint_reads_cnt_max i = 1; i < readsTotalCount; i += 2) {
                uint_reads_cnt_std idx = rlIdxOrder[i];
                if (idx < hqReadsCount)
                    hqPg->getReadsList()->revComp[idx] = !hqPg->getReadsList()->revComp[idx];
            }
        }
    }
    template void PgRCManager::applyRevComplPairFileToPgs<uint_pg_len_std>(vector<uint_pg_len_std> &orgIdx2PgPos);
    template void PgRCManager::applyRevComplPairFileToPgs<uint_pg_len_max>(vector<uint_pg_len_max> &orgIdx2PgPos);

    void PgRCManager::loadAllPgs(istream& pgrcIn) {
        chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
        PseudoGenomeHeader hqPgh(pgrcIn);
        ReadsSetProperties hqRsProp(pgrcIn);
        if (confirmTextReadMode(pgrcIn)) {
            cout << "Reads list text mode unsupported during decompression." << endl;
            exit(EXIT_FAILURE);
        }
        ExtendedReadsListWithConstantAccessOption* hqCaeRl =
                ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgrcIn,
                        &hqPgh, &hqRsProp, srcFastqFile.empty()?"":pgSeqFinalHqPrefix, preserveOrderMode);
        PseudoGenomeHeader lqPgh(pgrcIn);
        ReadsSetProperties lqRsProp(pgrcIn);
        if (confirmTextReadMode(pgrcIn)) {
            cout << "Reads list text mode unsupported during decompression." << endl;
            exit(EXIT_FAILURE);
        }
        ExtendedReadsListWithConstantAccessOption* lqCaeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgrcIn,
                    &lqPgh, &lqRsProp, srcFastqFile.empty()?"":pgSeqFinalLqPrefix, preserveOrderMode, true, true);
        ExtendedReadsListWithConstantAccessOption* nCaeRl = 0;
        ReadsSetProperties nRsProp;
        PseudoGenomeHeader nPgh;
        if (separateNReads) {
            nPgh = PseudoGenomeHeader(pgrcIn);
            nRsProp = ReadsSetProperties(pgrcIn);
            if (confirmTextReadMode(pgrcIn)) {
                cout << "Reads list text mode unsupported during decompression." << endl;
                exit(EXIT_FAILURE);
            }
            nCaeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgrcIn,
                    &nPgh, &nRsProp, srcFastqFile.empty()?"":pgNPrefix, preserveOrderMode, true, true);
        }
        readLength = hqRsProp.maxReadLength;
        hqReadsCount = hqRsProp.readsCount;
        lqReadsCount = lqRsProp.readsCount;
        nonNPgReadsCount = hqReadsCount + lqReadsCount;
        nPgReadsCount = separateNReads?nRsProp.readsCount:0;
        readsTotalCount = nonNPgReadsCount + nPgReadsCount;

        hqPgLen = hqPgh.getPseudoGenomeLength();
        nonNPgLen = hqPgLen + lqPgh.getPseudoGenomeLength();
        if (preserveOrderMode) {
            isJoinedPgLengthStd = nonNPgLen + nPgh.getPseudoGenomeLength() <= UINT32_MAX;
            if (isJoinedPgLengthStd)
                SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_std>(pgrcIn, orgIdx2StdPgPos,
                        readsTotalCount, singleReadsMode);
            else
                SeparatedPseudoGenomePersistence::decompressReadsPgPositions<uint_pg_len_max>(pgrcIn, orgIdx2PgPos,
                        readsTotalCount, singleReadsMode);
        } else {
            SeparatedPseudoGenomePersistence::decompressReadsOrder(pgrcIn, rlIdxOrder,
                                                                   preserveOrderMode, ignorePairOrderInformation, singleReadsMode);
        }
        cout << "... loaded Pgs Reads Lists (checkpoint: " << time_millis(start_t) << " msec.)" << endl;
        string hqPgSeq, lqPgSeq, nPgSeq;
        SimplePgMatcher::restoreMatchedPgs(pgrcIn, hqPgLen, hqPgSeq, lqPgSeq, nPgSeq);
        hqPg = new SeparatedPseudoGenome(move(hqPgSeq), hqCaeRl, &hqRsProp);
        lqPg = new SeparatedPseudoGenome(move(lqPgSeq), lqCaeRl, &lqRsProp);
        nPg = new SeparatedPseudoGenome(move(nPgSeq), nCaeRl, &nRsProp);
    }

    void PgRCManager::loadAllPgs() {
        string hqPgSeq = SimplePgMatcher::restoreAutoMatchedPg(pgSeqFinalHqPrefix, true);
        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgSeqFinalHqPrefix, pgh, prop, plainTextReadMode);
        ExtendedReadsListWithConstantAccessOption* caeRl =
                ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgSeqFinalHqPrefix, pgh->getPseudoGenomeLength());
        hqPg = new SeparatedPseudoGenome(move(hqPgSeq), caeRl, prop);
        delete(pgh);
        delete(prop);

        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgSeqFinalLqPrefix, pgh, prop, plainTextReadMode);
        string lqPgSeq = SimplePgMatcher::restoreMatchedPg(hqPg->getPgSequence(), pgSeqFinalLqPrefix, true, plainTextReadMode);
        caeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgSeqFinalLqPrefix, pgh->getPseudoGenomeLength());
        lqPg = new SeparatedPseudoGenome(move(lqPgSeq), caeRl, prop);
        delete(pgh);
        delete(prop);

        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgNPrefix, pgh, prop, plainTextReadMode);
        string nPgSeq = SimplePgMatcher::restoreMatchedPg(hqPg->getPgSequence(), pgNPrefix, true, plainTextReadMode);
        caeRl = ExtendedReadsListWithConstantAccessOption::loadConstantAccessExtendedReadsList(pgNPrefix, pgh->getPseudoGenomeLength());
        nPg = new SeparatedPseudoGenome(move(nPgSeq), caeRl, prop);
        delete(pgh);
        delete(prop);

        readLength = hqPg->getReadsSetProperties()->maxReadLength;
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