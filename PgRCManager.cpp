#include <sys/stat.h>
#include "PgRCManager.h"

#include "matching/ReadsMatchers.h"
#include "matching/SimplePgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"

namespace PgTools {

    static const char *const BAD_INFIX = "bad";
    static const char *const GOOD_INFIX = "good";
    static const char *const N_INFIX = "N";
    static const char *const DIVISION_EXTENSION = ".div";

    static const char *const TEMPORARY_FILE_SUFFIX = ".temp";
    static const char *const PGRC_HEADER = "PgRC";

    uint_read_len_max probeReadsLength(const string &srcFastqFile);
    clock_t getTimeInSec(clock_t end_t, clock_t begin_t) { return ((end_t - begin_t) / CLOCKS_PER_SEC); }

    void PgRCManager::prepareChainData() {
        qualityDivision = error_limit_in_promils < 1000;
        readLength = probeReadsLength(srcFastqFile);
        if (pairFastqFile.empty() && !preserveOrderMode)
            singleReadsMode = true;

        if (std::ifstream(pgRCFileName))
            fprintf(stderr, "Warning: file %s already exists\n", pgRCFileName.data());
        pgrcOut = fstream(pgRCFileName + TEMPORARY_FILE_SUFFIX, ios::out | ios::binary | ios::trunc);
        pgrcOut.write(PGRC_HEADER, strlen(PGRC_HEADER));
        pgrcOut.put(singleReadsMode?PGRC_SE_MODE:
            (preserveOrderMode?(pairFastqFile.empty()?PGRC_ORD_SE_MODE:PGRC_ORD_PE_MODE):PGRC_PE_MODE));
        pgrcOut.put(separateNReads);

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
        start_t = clock();
        prepareChainData();
        stageCount = 0;
        if (skipStages < ++stageCount && qualityDivision) {
            runQualityBasedDivision();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistReadsQualityDivision();
                disposeChainData();
            }
        }
        div_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForPgGeneratorBaseReadsDivision();
            runPgGeneratorBasedReadsDivision();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistReadsQualityDivision();
                disposeChainData();
            }
        }
        pgDiv_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForHqPgGeneration();
            runHQPgGeneration();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistHQPg();
                persistReadsQualityDivision();
                disposeChainData();
            }
        }
        good_t = clock();
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
                if (!singleReadsMode)
                    orgIdxs = std::move(hqPg->getReadsList()->orgIdx);
                hqPg->disposeReadsList();
            }
        }
        match_t = clock();
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
            if (!singleReadsMode)
                orgIdxs.insert(orgIdxs.end(), lqPg->getReadsList()->orgIdx.begin(), lqPg->getReadsList()->orgIdx.end());
            lqPg->disposeReadsList();
            if (separateNReads) {
                runNPgGeneration();
                persistNPg();
                if (!singleReadsMode)
                    orgIdxs.insert(orgIdxs.end(), nPg->getReadsList()->orgIdx.begin(), nPg->getReadsList()->orgIdx.end());
                nPg->disposeReadsList();
            }
            delete(divReadsSets);
            divReadsSets = 0;
        }
        bad_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForPgMatching();
            string emptySequence;
            //DefaultPgMatcher::matchPgInPgFile(pgMappedHqPrefix, pgMappedHqPrefix, readsLength, pgHqPrefix, true, false);
            SimplePgMatcher::matchPgsInPg(hqPg->getPgSequence(), lqPg->getPgSequence(),
                    separateNReads?nPg->getPgSequence():emptySequence, pgrcOut, compressionLevel,
                                          extraFilesForValidation?pgSeqFinalHqPrefix:"",
                                          extraFilesForValidation?pgSeqFinalLqPrefix:"",
                                          extraFilesForValidation?pgNPrefix:"", targetPgMatchLength);
//            testCompressSequences();
        }
        gooder_t = clock();
        if (!singleReadsMode && skipStages < ++stageCount && endAtStage >= stageCount) {
            SeparatedPseudoGenomePersistence::compressReadsOrder(pgrcOut, orgIdxs, compressionLevel, preserveOrderMode,
                    pairFastqFile.empty());
            if (extraFilesForValidation) {
                if (separateNReads)
                    SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedHqPrefix, pgMappedLqPrefix, pgNPrefix});
                else
                    SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedHqPrefix, pgMappedLqPrefix});
            }
        }
        finalizeCompression();
        disposeChainData();
        reportTimes();
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
        divReadsSets->getHqReadsSet()->printout();
        hqPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                divReadsSets->getHqReadsSet());
        divReadsSets->disposeHqReadsSet();
        IndexesMapping* hq2IndexesMapping = divReadsSets->generateHqReadsIndexesMapping();
        hqPg->applyIndexesMapping(hq2IndexesMapping);
        delete(hq2IndexesMapping);
        if (revComplPairFile)
            hqPg->applyRevComplPairFile();
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
        divReadsSets->getLqReadsSet()->printout();
        if (separateNReads)
            divReadsSets->getNReadsSet()->printout();
        ConstantLengthReadsSetInterface* readsSet =
                separateNReads?((ConstantLengthReadsSetInterface*) new SumOfConstantLengthReadsSets(
                        divReadsSets->getLqReadsSet(), divReadsSets->getNReadsSet())): divReadsSets->getLqReadsSet();
        IndexesMapping* mapping = separateNReads?((IndexesMapping*) new SumOfMappings(
                divReadsSets->getLqReadsIndexesMapping(), divReadsSets->getNReadsIndexesMapping())):divReadsSets->getLqReadsIndexesMapping();
        const vector<bool>& isReadMappedIntoHqPg = mapReadsIntoPg(
                hqPg, true, readsSet, DefaultReadsMatcher::DISABLED_PREFIX_MODE,
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
        divReadsSets->getLqReadsSet()->printout();
        lqPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                divReadsSets->getLqReadsSet());
        lqPg->applyIndexesMapping(divReadsSets->getLqReadsIndexesMapping());
        divReadsSets->disposeLqReadsSet();
        if (revComplPairFile)
            lqPg->applyRevComplPairFile();
    }

    void PgRCManager::persistLQPg() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(lqPg, pgMappedLqPrefix);
    }

    void PgRCManager::compressLQPgReadsList() {
        SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(lqPg, &pgrcOut, compressionLevel);
        if (extraFilesForValidation)
            SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(lqPg, pgMappedLqPrefix, true);
    }

    void PgRCManager::runNPgGeneration() {
        divReadsSets->getNReadsSet()->printout();
        nPg = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
                divReadsSets->getNReadsSet());
        nPg->applyIndexesMapping(divReadsSets->getNReadsIndexesMapping());
        divReadsSets->disposeNReadsSet();
        if (revComplPairFile)
            nPg->applyRevComplPairFile();
    }

    void PgRCManager::persistNPg() {
        SeparatedPseudoGenomePersistence::compressSeparatedPseudoGenomeReadsList(nPg, &pgrcOut, compressionLevel);
        if (extraFilesForValidation)
            SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(nPg, pgNPrefix);
    }

    void PgRCManager::prepareForPgMatching() {
        if (!hqPg)
            hqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgHqPrefix, true);
        if (!lqPg)
            lqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgMappedLqPrefix, true);
    }


    void PgRCManager::compressMEMMappedPgSequences() {
        cout << "Error: unimplemented standalone compressMEMMappedPgSequences!" << endl;
        exit(EXIT_FAILURE);
    }

    void PgRCManager::reportTimes() {
        string outputfile = "pgrc_res.txt";
        bool hasHeader = (bool) std::ifstream(outputfile);
        fstream fout(outputfile, ios::out | ios::binary | ios::app);
        if (!hasHeader)
            fout << "srcFastq\tpairFastq\trcPairFile\tpgPrefix\tq[%o]\tg[%o]\tm\tM\tp\ttotal[s]\tdiv[s]\tPgDiv[s]\tgood[s]\treadsMatch[s]\tbad[s]\tpgMatch[s]\tpost[s]" << endl;

        fout << srcFastqFile << "\t" << pairFastqFile << "\t" << (revComplPairFile?"yes":"no") << "\t"
             << pgRCFileName << "\t" << toString(error_limit_in_promils) << "\t" << gen_quality_str << "\t";

        if (preReadsExactMatchingChars > 0)
            fout << (char) tolower(preMatchingMode) << ((toupper(preMatchingMode) == preMatchingMode)?string("s"):string("")) << (int) preReadsExactMatchingChars;
        fout << (char) tolower(matchingMode) << ((toupper(matchingMode) == matchingMode)?string("s"):string("")) << (int) readsExactMatchingChars << "\t" << (int) minCharsPerMismatch << "\t" << targetPgMatchLength << "\t";
        fout << getTimeInSec(clock(), start_t) << "\t";
        fout << getTimeInSec(div_t, start_t) << "\t";
        fout << getTimeInSec(pgDiv_t, div_t) << "\t";
        fout << getTimeInSec(good_t, pgDiv_t) << "\t";
        fout << getTimeInSec(match_t, good_t) << "\t";
        fout << getTimeInSec(bad_t, match_t) << "\t";
        fout << getTimeInSec(gooder_t, bad_t) << "\t";
        fout << getTimeInSec(clock(), gooder_t) << endl;
    }

    void PgRCManager::finalizeCompression() {
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
        start_t = clock();
        string tmpDirectoryPath = pgRCFileName + "/";
        ifstream pgrcIn(pgRCFileName);
        if (pgrcIn) {
            for(int i = 0; i < strlen(PGRC_HEADER); i++) {
                if (PGRC_HEADER[i] != pgrcIn.get()) {
                    fprintf(stderr, "Error processing header.\n");
                    exit(EXIT_FAILURE);
                }
            }
            char pgrc_mode = pgrcIn.get();
            if (pgrc_mode != PGRC_SE_MODE) {
                fprintf(stderr, "Unsupported decompression mode: %d.\n", pgrc_mode);
                exit(EXIT_FAILURE);
            }
            separateNReads = (bool) pgrcIn.get();

            pgrcIn >> tmpDirectoryPath;
            tmpDirectoryPath = tmpDirectoryPath + "/";
            pgrcIn.get();
        }

        pgMappedHqPrefix = tmpDirectoryPath + GOOD_INFIX;
        pgMappedLqPrefix = tmpDirectoryPath + BAD_INFIX;
        pgSeqFinalHqPrefix = tmpDirectoryPath + GOOD_INFIX;
        pgSeqFinalLqPrefix = tmpDirectoryPath + BAD_INFIX;
        pgNPrefix = tmpDirectoryPath + N_INFIX;

        if (pgrcIn)
            loadAllPgs(pgrcIn);
        else
            loadAllPgs();
        cout << "... loaded Pgs (checkpoint: " << clock_millis(start_t) << " msec.)" << endl;

        if (srcFastqFile.empty()) {
            if (ENABLE_PARALLEL_DECOMPRESSION && dnaStreamSize() < CHUNK_SIZE_IN_BYTES)
                writeAllReadsInSEModeParallel(tmpDirectoryPath);
            else
                writeAllReadsInSEMode(tmpDirectoryPath);
            cout << "Decompressed ";
        } else {
            validateAllPgs();
            cout << "Validated ";
        }

        cout << (hqPg->getReadsSetProperties()->readsCount + lqPg->getReadsSetProperties()->readsCount +
                (nPg?nPg->getReadsSetProperties()->readsCount:0)) <<
                " reads in " << clock_millis(start_t) << " msec." << endl;

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

    void PgRCManager::writeAllReadsInSEModeParallel(const string &tmpDirectoryPath) {
        cout << "... parallel mode" << endl;
        std::thread writing(&PgRCManager::writeFromQueue, this, tmpDirectoryPath);
        string res;
        uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
        res.reserve(CHUNK_SIZE_IN_BYTES);
        for(uint_reads_cnt_max i = 0; i < hqPg->getReadsSetProperties()->readsCount; i++) {
            if (res.size() > res_size_guard) {
                pushOutToQueue(res);
            }
            res.append(hqPg->getRead(i));
            res.push_back('\n');
        }
        for(uint_reads_cnt_max i = 0; i < lqPg->getReadsSetProperties()->readsCount; i++) {
            if (res.size() > res_size_guard) {
                pushOutToQueue(res);
            }
            res.append(lqPg->getRead(i));
            res.push_back('\n');
        }
        uint_reads_cnt_max nPgReadsCount = nPg?nPg->getReadsSetProperties()->readsCount:0;
        for(uint_reads_cnt_max i = 0; i < nPgReadsCount; i++) {
            if (res.size() > res_size_guard) {
                pushOutToQueue(res);
            }
            res.append(nPg->getRead(i));
            res.push_back('\n');
        }
        pushOutToQueue(res);
        cout << "... finished loading queue (checkpoint: " << clock_millis(start_t) << " msec.)" << endl;
        finishWritingParallel();
        writing.join();
    }

    void PgRCManager::writeFromQueue(const string &tmpDirectoryPath) {
        fstream fout(tmpDirectoryPath + "out", ios_base::out | ios_base::binary | std::ios::trunc);
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

    void PgRCManager::writeAllReadsInSEMode(const string &tmpDirectoryPath) const {
        fstream fout(tmpDirectoryPath + "out", ios_base::out | ios_base::binary | std::ios::trunc);
        string res;
        uint64_t res_size_guard = CHUNK_SIZE_IN_BYTES;
        uint64_t totalSize = dnaStreamSize();
        res.reserve(totalSize < res_size_guard?totalSize:res_size_guard + (hqPg->getReadsSetProperties()->maxReadLength + 1));
        for(uint_reads_cnt_max i = 0; i < hqPg->getReadsSetProperties()->readsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            res.append(hqPg->getRead(i));
            res.push_back('\n');
        }
        for(uint_reads_cnt_max i = 0; i < lqPg->getReadsSetProperties()->readsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            res.append(lqPg->getRead(i));
            res.push_back('\n');
        }
        uint_reads_cnt_max nPgReadsCount = nPg?nPg->getReadsSetProperties()->readsCount:0;
        for(uint_reads_cnt_max i = 0; i < nPgReadsCount; i++) {
            if (res.size() > res_size_guard) {
                fout << res;
                res.resize(0);
            }
            res.append(nPg->getRead(i));
            res.push_back('\n');
        }
        fout << res;
        fout.close();
    }

    uint_reads_cnt_max PgRCManager::dnaStreamSize() const {
        return (hqPg->getReadsSetProperties()->readsCount + lqPg->getReadsSetProperties()->readsCount) *
               (hqPg->getReadsSetProperties()->maxReadLength + 1);
    }

    void PgRCManager::validateAllPgs() {
        uint_reads_cnt_max nonNPgReadsCount = hqPg->getReadsSetProperties()->readsCount
                                              + lqPg->getReadsSetProperties()->readsCount;
        uint_reads_cnt_max nPgReadsCount = nPg?nPg->getReadsSetProperties()->readsCount:0;
        uint_reads_cnt_max readsTotalCount = nonNPgReadsCount + nPgReadsCount;
        vector<uint_reads_cnt_max> orgIdx2rlIdx;
        orgIdx2rlIdx.resize(readsTotalCount);
        for(uint_reads_cnt_max i = 0; i < hqPg->getReadsSetProperties()->readsCount; i++)
            orgIdx2rlIdx[hqPg->getReadsList()->orgIdx[i]] = i;
        for(uint_reads_cnt_max i = 0; i < lqPg->getReadsSetProperties()->readsCount; i++)
            orgIdx2rlIdx[lqPg->getReadsList()->orgIdx[i]] = hqPg->getReadsSetProperties()->readsCount + i;
        for(uint_reads_cnt_max i = 0; i < nPgReadsCount; i++)
            orgIdx2rlIdx[nPg->getReadsList()->orgIdx[i]] = nonNPgReadsCount + i;

        ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
                srcFastqFile, pairFastqFile, revComplPairFile);

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
            if (validated[idx]) {
                notValidatedCount++;
                continue;
            }
            validated[idx] = true;
            if (idx < hqPg->getReadsSetProperties()->readsCount)
                read = hqPg->getRead(idx);
            else if (idx < nonNPgReadsCount)
                read = lqPg->getRead(idx - hqPg->getReadsSetProperties()->readsCount);
            else
                read = nPg->getRead(idx - nonNPgReadsCount);
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

    void PgRCManager::loadAllPgs(istream& pgrcIn) {
        PseudoGenomeHeader hqPgh(pgrcIn);
        ReadsSetProperties hqRsProp(pgrcIn);
        if (confirmTextReadMode(pgrcIn)) {
            cout << "Reads list text mode unsupported during decompression." << endl;
            exit(EXIT_FAILURE);
        }
        ConstantAccessExtendedReadsList* hqCaeRl =
                ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(pgrcIn,
                        &hqPgh, &hqRsProp, srcFastqFile.empty()?"":pgSeqFinalHqPrefix);
        PseudoGenomeHeader lqPgh(pgrcIn);
        ReadsSetProperties lqRsProp(pgrcIn);
        if (confirmTextReadMode(pgrcIn)) {
            cout << "Reads list text mode unsupported during decompression." << endl;
            exit(EXIT_FAILURE);
        }
        ConstantAccessExtendedReadsList* lqCaeRl =
                ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(pgrcIn,
                        &lqPgh, &lqRsProp, srcFastqFile.empty()?"":pgSeqFinalLqPrefix);
        ConstantAccessExtendedReadsList* nCaeRl = 0;
        ReadsSetProperties nRsProp;
        if (separateNReads) {
            PseudoGenomeHeader nPgh(pgrcIn);
            nRsProp = ReadsSetProperties(pgrcIn);
            if (confirmTextReadMode(pgrcIn)) {
                cout << "Reads list text mode unsupported during decompression." << endl;
                exit(EXIT_FAILURE);
            }
            nCaeRl = ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(pgrcIn,
                            &nPgh, &nRsProp, srcFastqFile.empty()?"":pgNPrefix);
        }
        string hqPgSeq, lqPgSeq, nPgSeq;
        SimplePgMatcher::restoreMatchedPgs(pgrcIn, hqPgSeq, lqPgSeq, nPgSeq);
        hqPg = new SeparatedPseudoGenome(move(hqPgSeq), hqCaeRl, &hqRsProp);
        lqPg = new SeparatedPseudoGenome(move(lqPgSeq), lqCaeRl, &lqRsProp);
        nPg = new SeparatedPseudoGenome(move(nPgSeq), nCaeRl, &nRsProp);
    }

    void PgRCManager::loadAllPgs() {
        string hqPgSeq = SimplePgMatcher::restoreAutoMatchedPg(pgSeqFinalHqPrefix, true);

        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgSeqFinalHqPrefix, pgh, prop, plainTextReadMode);
        ConstantAccessExtendedReadsList* caeRl =
                ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(pgSeqFinalHqPrefix, pgh->getPseudoGenomeLength());
        hqPg = new SeparatedPseudoGenome(move(hqPgSeq), caeRl, prop);
        delete(pgh);
        delete(prop);

        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgSeqFinalLqPrefix, pgh, prop, plainTextReadMode);

        string lqPgSeq = SimplePgMatcher::restoreMatchedPg(hqPg->getPgSequence(), pgSeqFinalLqPrefix, true, plainTextReadMode);
        caeRl = ConstantAccessExtendedReadsList::loadConstantAccessExtendedReadsList(pgSeqFinalLqPrefix, pgh->getPseudoGenomeLength());
        lqPg = new SeparatedPseudoGenome(move(lqPgSeq), caeRl, prop);
        delete(pgh);
        delete(prop);
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