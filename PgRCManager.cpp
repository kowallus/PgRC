#include "PgRCManager.h"

#include "matching/ReadsMatchers.h"
#include "matching/SimplePgMatcher.h"
#include "pseudogenome/TemplateUserGenerator.h"
#include "readsset/persistance/ReadsSetPersistence.h"
#include "pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "pseudogenome/persistence/PseudoGenomePersistence.h"
#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    static const char *const BAD_INFIX = "_bad";
    static const char *const GOOD_INFIX = "_good";
    static const char *const N_INFIX = "_N";
    static const char *const DIVISION_EXTENSION = ".div";

    uint_read_len_max probeReadsLength(const string &srcFastqFile);
    clock_t getTimeInSec(clock_t end_t, clock_t begin_t) { return ((end_t - begin_t) / CLOCKS_PER_SEC); }


    void PgRCManager::prepareChainData() {
        qualityDivision = error_limit_in_promils < 1000;
        readLength = probeReadsLength(srcFastqFile);
        targetMismatches = readLength / targetCharsPerMismatch;
        maxMismatches = readLength / maxCharsPerMismatch;

        if (qualityDivision)
            pgFilesPrefixes = pgFilesPrefixes + "_q" + toString(error_limit_in_promils);
        pgFilesPrefixes = pgFilesPrefixes + (nReadsLQ?"_n":"") + (separateNReads?"_N":"") + "_g" + gen_quality_str;
        lqDivisionFile = pgFilesPrefixes + BAD_INFIX + DIVISION_EXTENSION;
        nDivisionFile = pgFilesPrefixes + N_INFIX + DIVISION_EXTENSION;
        pgGoodPrefix = pgFilesPrefixes + GOOD_INFIX;
        pgFilesPrefixesWithM = pgFilesPrefixes + "_m" + toString(targetCharsPerMismatch)
                                      + "_M" + mismatchesMode + toString(maxCharsPerMismatch) + "_p" + toString(minimalPgMatchLength);
        pgMappedGoodPrefix = pgFilesPrefixesWithM + GOOD_INFIX;
        pgMappedBadPrefix = pgFilesPrefixesWithM + BAD_INFIX;
        pgNPrefix = pgFilesPrefixesWithM + N_INFIX;
        mappedBadDivisionFile = pgFilesPrefixesWithM + DIVISION_EXTENSION;
        if (skipIntermediateOutput) {
            pgGoodPrefix = pgMappedGoodPrefix;
            lqDivisionFile = mappedBadDivisionFile;
        }
    }

    void PgRCManager::executePgRCChain() {
        start_t = clock();
        prepareChainData();
        stageCount = 0;
        if (skipStages < ++stageCount && qualityDivision) {
            runQualityBasedDivision();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistQualityBasedDivision();
                disposeChainData();
            }
        }
        div_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForPgGeneratorBaseReadsDivision();
            runPgGeneratorBasedReadsDivision();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistPgGeneratorBasedReadsDivision();
                disposeChainData();
            }
        }
        pgDiv_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForHqPgGeneration();
            runHQPgGeneration();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistHQPg();
                disposeChainData();
            }
        }
        good_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForMappingLQReadsOnHQPg();
            runMappingLQReadsOnHQPg();
            if (disableInMemoryMode || endAtStage == stageCount) {
                persistHQPg();
                disposeChainData();
            } else {
                saveHQPgReadsList();
                extractHQPgSequence();
                freeHQPg();
            }
        }
        match_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            prepareForLQPgAndNPgGeneration();
            if (separateNReads) {
                runNPgGeneration();
                saveNPg();
            }
            runLQPgGeneration();
            if (disableInMemoryMode || endAtStage == stageCount) {
                saveLQPg();
                if (skipStages < stageCount - 1)
                    saveHQPgSequence();
            } else {
                saveLQPgReadsList();
                extractLQPgSequence();
                freeLQPg();
            }
            delete(divReadsSets);
            divReadsSets = 0;
        }
        bad_t = clock();
        if (skipStages < ++stageCount && endAtStage >= stageCount) {
            //DefaultPgMatcher::matchPgInPgFile(pgMappedGoodPrefix, pgMappedGoodPrefix, readsLength, pgGoodPrefix, true, false);
            SimplePgMatcher::matchPgInPgFiles(pgMappedGoodPrefix, pgMappedBadPrefix, minimalPgMatchLength, true);
            saveMEMMappedPgSequences();
        }
        gooder_t = clock();
        if (pairFastqFile != "" && skipStages < ++stageCount && endAtStage >= stageCount) {
            if (separateNReads)
                SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedGoodPrefix, pgMappedBadPrefix, pgNPrefix});
            else
                SeparatedPseudoGenomePersistence::dumpPgPairs({pgMappedGoodPrefix, pgMappedBadPrefix});
        }
        disposeChainData();
        reportTimes();
    }

    void PgRCManager::runQualityBasedDivision() {
        ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(
                srcFastqFile, pairFastqFile);
        divReadsSets =
                DividedPCLReadsSets::getQualityDivisionBasedReadsSets(allReadsIterator, readLength, error_limit_in_promils / 1000.0,
                        separateNReads, nReadsLQ);
        delete (allReadsIterator);
    }

    void PgRCManager::persistQualityBasedDivision() {
        divReadsSets->getLqReadsIndexesMapping()->saveMapping(lqDivisionFile);
        if (separateNReads)
            divReadsSets->getNReadsIndexesMapping()->saveMapping(nDivisionFile);
    }

    void PgRCManager::prepareForPgGeneratorBaseReadsDivision() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(
                    srcFastqFile, pairFastqFile);
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

    void PgRCManager::persistPgGeneratorBasedReadsDivision() {
        divReadsSets->getLqReadsIndexesMapping()->saveMapping(lqDivisionFile);
    }

    void PgRCManager::prepareForHqPgGeneration() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(
                    srcFastqFile, pairFastqFile);
            divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, readLength, lqDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "");
            delete (allReadsIterator);
        }
    }

    void PgRCManager::runHQPgGeneration() {
        divReadsSets->getHqReadsSet()->printout();
        PseudoGenomeBase *goodPgb = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(
                divReadsSets->getHqReadsSet());
        IndexesMapping* good2IndexesMapping = divReadsSets->generateHqReadsIndexesMapping();
        SeparatedPseudoGenomePersistence::writePseudoGenome(goodPgb, pgGoodPrefix, good2IndexesMapping,
                                                            revComplPairFile);
        delete(good2IndexesMapping);
    }

    void PgRCManager::persistHQPg() {

    }

    void PgRCManager::prepareForMappingLQReadsOnHQPg() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(
                    srcFastqFile, pairFastqFile);
            divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, readLength, lqDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "");
            delete (allReadsIterator);
        }
    }

    void PgRCManager::runMappingLQReadsOnHQPg() {
        divReadsSets->getLqReadsSet()->printout();
        const vector<bool>& isLqReadMappedIntoHqPg = mapReadsIntoPg(
                pgGoodPrefix, true, divReadsSets->getLqReadsSet(), DefaultReadsMatcher::DISABLED_PREFIX_MODE,
                targetMismatches, maxMismatches, mismatchesMode, 0,
                false, pgMappedGoodPrefix, divReadsSets->getLqReadsIndexesMapping());
        divReadsSets->removeReadsFromLqReadsSet(isLqReadMappedIntoHqPg);
        divReadsSets->getLqReadsIndexesMapping()->saveMapping(mappedBadDivisionFile);
    }

    void PgRCManager::saveHQPgReadsList() {

    }

    void PgRCManager::saveHQPgSequence() {

    }

    void PgRCManager::extractHQPgSequence() {

    }

    void PgRCManager::freeHQPg() {

    }

    void PgRCManager::prepareForLQPgAndNPgGeneration() {
        if (!divReadsSets) {
            ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(
                    srcFastqFile, pairFastqFile);
            divReadsSets = DividedPCLReadsSets::loadDivisionReadsSets(
                    allReadsIterator, readLength, mappedBadDivisionFile, nReadsLQ, separateNReads ? nDivisionFile : "", true);
            delete (allReadsIterator);
        }
    }

    void PgRCManager::runLQPgGeneration() {
        divReadsSets->getLqReadsSet()->printout();
        PseudoGenomeBase *goodPgb = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(
                divReadsSets->getLqReadsSet());
        SeparatedPseudoGenomePersistence::writePseudoGenome(goodPgb, pgMappedBadPrefix,
                divReadsSets->getLqReadsIndexesMapping(), revComplPairFile);
    }

    void PgRCManager::saveLQPg() {

    }

    void PgRCManager::saveLQPgReadsList() {

    }

    void PgRCManager::extractLQPgSequence() {

    }

    void PgRCManager::freeLQPg() {

    }

    void PgRCManager::runNPgGeneration() {
        divReadsSets->getNReadsSet()->printout();
        PseudoGenomeBase *goodPgb = GreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(
                divReadsSets->getNReadsSet());
        SeparatedPseudoGenomePersistence::writePseudoGenome(goodPgb, pgNPrefix,
                                                            divReadsSets->getNReadsIndexesMapping(), revComplPairFile);
    }

    void PgRCManager::saveNPg() {

    }

    void PgRCManager::saveMEMMappedPgSequences() {

    }

    void PgRCManager::reportTimes() {
        string outputfile = "pgrc_res.txt";
        bool hasHeader = (bool) std::ifstream(outputfile);
        fstream fout(outputfile, ios::out | ios::binary | ios::app);
        if (!hasHeader)
            fout << "srcFastq\tpairFastq\trcPairFile\tpgPrefix\tq[%o]\tg[%o]\tm\tM\tp\ttotal[s]\tdiv[s]\tPgDiv[s]\tgood[s]\treadsMatch[s]\tbad[s]\tpgMatch[s]\tpost[s]" << endl;

        fout << srcFastqFile << "\t" << pairFastqFile << "\t" << (revComplPairFile?"yes":"no") << "\t"
             << pgFilesPrefixes << "\t" << toString(error_limit_in_promils) << "\t" << gen_quality_str << "\t"
             << (int) targetMismatches << "\t" << mismatchesMode << (int) maxMismatches << "\t" << minimalPgMatchLength << "\t";
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