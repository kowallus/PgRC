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

    uint_read_len_max probeReadsLength(const string &srcFastqFile);
    clock_t getTimeInSec(clock_t end_t, clock_t begin_t) { return ((end_t - begin_t) / CLOCKS_PER_SEC); }

    void PgRCManager::prepareChainData() {
        qualityDivision = error_limit_in_promils < 1000;
        readLength = probeReadsLength(srcFastqFile);

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


        mode_t mode = 0777;
        int nError = mkdir(tmpDirectoryName.data(), mode);
        if (nError != 0) {
            srand (time(NULL));
            tmpDirectoryName = tmpDirectoryName + "_" + toString(rand()%100000);
            nError = mkdir(tmpDirectoryName.data(), mode);
            if (nError != 0) {
                fprintf(stderr, "Error creating folder %s\n", tmpDirectoryName.data());
                exit(EXIT_FAILURE);
            }
        }

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
//                persistMappedHQPgReadsList();
                hqPg->disposeReadsList();
            }
        }
        match_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForLQPgAndNPgGeneration();
            if (separateNReads) {
                runNPgGeneration();
                persistNPg();
                delete(nPg);
                nPg = 0;
            }
            runLQPgGeneration();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistLQPg();
                if (hqPg)
                    persistHQPgSequence();
            } else {
                persistLQPgReadsList();
            }
            delete(divReadsSets);
            divReadsSets = 0;
            lqPg->disposeReadsList();
        }
        bad_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForPgMatching();
            //DefaultPgMatcher::matchPgInPgFile(pgMappedHqPrefix, pgMappedHqPrefix, readsLength, pgHqPrefix, true, false);
            SimplePgMatcher::matchPgInPgFiles(hqPg->getPgSequence(), lqPg->getPgSequence(),
                                              pgSeqFinalHqPrefix, pgSeqFinalLqPrefix, targetPgMatchLength);
        }
        gooder_t = clock();
        if (pairFastqFile != "" && skipStages < ++stageCount && endAtStage >= stageCount) {
            if (separateNReads)
                SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedHqPrefix, pgMappedLqPrefix, pgNPrefix});
            else
                SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedHqPrefix, pgMappedLqPrefix});
        }
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
        const vector<bool>& isLqReadMappedIntoHqPg = mapReadsIntoPg(
                hqPg, true, divReadsSets->getLqReadsSet(), DefaultReadsMatcher::DISABLED_PREFIX_MODE,
                preReadsExactMatchingChars, readsExactMatchingChars,
                minCharsPerMismatch, preMatchingMode, matchingMode,
                false, pgMappedHqPrefix, divReadsSets->getLqReadsIndexesMapping());
        divReadsSets->removeReadsFromLqReadsSet(isLqReadMappedIntoHqPg);
    }

    void PgRCManager::persistMappedReadsQualityDivision() {
        divReadsSets->getLqReadsIndexesMapping()->saveMapping(mappedLqDivisionFile);
        if (separateNReads)
            divReadsSets->getNReadsIndexesMapping()->saveMapping(nDivisionFile);
    }

    void PgRCManager::persistMappedHQPgReadsList() {
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(hqPg, pgMappedHqPrefix, true);
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

    void PgRCManager::persistLQPgReadsList() {
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
        SeparatedPseudoGenomePersistence::writeSeparatedPseudoGenome(nPg, pgNPrefix);
    }

    void PgRCManager::prepareForPgMatching() {
        if (!hqPg)
            hqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgHqPrefix, true);
        if (!lqPg)
            lqPg = SeparatedPseudoGenomePersistence::loadSeparatedPseudoGenome(pgMappedLqPrefix, true);
    }


    void PgRCManager::persistMEMMappedPgSequences() {

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

    uint_read_len_max probeReadsLength(const string &srcFastqFile) {
        ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt = ReadsSetPersistence::createManagedReadsIterator(
                srcFastqFile);
        readsIt->moveNext();
        uint_read_len_max readsLength = readsIt->getReadLength();
        delete(readsIt);
        return readsLength;
    }


}