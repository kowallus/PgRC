#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"

#include <cstdio>

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, const string &pseudoGenomePrefix,
            string divisionFile, bool divisionComplement, bool revComplPairFile) {
        clock_checkpoint();
        SeparatedPseudoGenomeOutputBuilder builder(pseudoGenomePrefix, !revComplPairFile, true);
        builder.writePseudoGenome(pgb, divisionFile, divisionComplement, revComplPairFile);
        builder.build();
        cout << "Writing (" << pseudoGenomePrefix << ") pseudo genome files in " << clock_millis() << " msec." << endl << endl;
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

    string SeparatedPseudoGenomePersistence::getPseudoGenome(const string &pseudoGenomePrefix) {
        std::ifstream pgSrc = getPseudoGenomeSrc(pseudoGenomePrefix);
        pgSrc.seekg(0, std::ios::end);
        size_t size = pgSrc.tellg();
        std::string buffer(size, ' ');
        pgSrc.seekg(0);
        pgSrc.read(&buffer[0], size);
        return buffer;
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

        std::ofstream destPgEl(pgElFile, std::ios::out | std::ios::binary);
        return destPgEl;
    }


    bool SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElement(const string &pseudoGenomePrefix,
                                                                              const string &fileSuffix) {
        string pgElFile = pseudoGenomePrefix + fileSuffix;
        string pgElTempFile = pgElFile + TEMPORARY_FILE_SUFFIX;
        if (std::ifstream(pgElTempFile)) {
            if (std::ifstream(pgElFile))
                remove(pgElFile.c_str());
            return rename(pgElTempFile.c_str(), pgElFile.c_str()) == 0;
        }
        return false;
    }

    void SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(const string &pseudoGenomePrefix) {
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, PSEUDOGENOME_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, PSEUDOGENOME_PROPERTIES_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_POSITIONS_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_REVERSECOMPL_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHES_COUNT_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);

        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_OFFSETS_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
    }

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX = ".pg";
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX = ".pg.prop";
    const string SeparatedPseudoGenomePersistence::READSLIST_POSITIONS_FILE_SUFFIX = ".pg.rl.pos";
    const string SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = ".pg.rl.idx";
    const string SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX = ".pg.rl.rc";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_COUNT_FILE_SUFFIX = ".pg.rl.mis.cnt";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHED_SYMBOLS_FILE_SUFFIX = ".pg.rl.mis.sym";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX = ".pg.rl.mis.pos";


    const string SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX = ".pg.rl.off";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX = ".pg.rl.mis.roff";

    const string SeparatedPseudoGenomePersistence::TEMPORARY_FILE_SUFFIX = ".temp";

    bool SeparatedPseudoGenomePersistence::enableReadPositionRepresentation = false;
    bool SeparatedPseudoGenomePersistence::enableRevOffsetMismatchesRepresentation = false;

    SeparatedPseudoGenomeOutputBuilder::SeparatedPseudoGenomeOutputBuilder(const string &pseudoGenomePrefix,
            bool disableRevComp, bool disableMismatches) : pseudoGenomePrefix(pseudoGenomePrefix),
            disableRevComp(disableRevComp), disableMismatches(disableMismatches) {
        initReadsListDests();
    }

    void SeparatedPseudoGenomeOutputBuilder::initDest(ofstream *&dest, const string &fileSuffix) {
        if (dest == 0)
            dest = new ofstream(SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(pseudoGenomePrefix, fileSuffix, true));
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
                initDest(rlMisRevOffDest, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_REVOFFSETS_FILE_SUFFIX);
            else
                initDest(rlMisPosDest, SeparatedPseudoGenomePersistence::READSLIST_MISMATCHES_POSITIONS_FILE_SUFFIX);
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::freeDest(ofstream* &dest) {
        if (dest) {
            dest->close();
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

    void SeparatedPseudoGenomeOutputBuilder::build() {
        if (pgh == 0) {
            fprintf(stderr, "Pseudo genome header not initialized in separated Pg builder.\n");
            exit(EXIT_FAILURE);

        }
        pgh->setReadsCount(readsCounter);
        initDest(pgPropDest, SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX);
        pgh->write(*pgPropDest);
        writeReadMode(*pgPropDest, plainTextWriteMode);

        freeDests();
        SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(pseudoGenomePrefix);
    }

    void SeparatedPseudoGenomeOutputBuilder::writeReadEntry(const DefaultReadsListEntry &rlEntry) {
        if (SeparatedPseudoGenomePersistence::enableReadPositionRepresentation)
            PgSAHelpers::writeValue<uint_pg_len_max>(*rlPosDest, rlEntry.pos);
        else
            PgSAHelpers::writeValue<uint_read_len_max>(*rlOffDest, rlEntry.offset);
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
                        PgSAHelpers::writeValue<uint_read_len_max>(*rlMisRevOffDest,
                                                                   currentPos - rlEntry.mismatchOffset[i]);
                        currentPos = rlEntry.mismatchOffset[i] - 1;
                    }
                } else {
                    for (uint8_t i = 0; i < rlEntry.mismatchesCount; i++)
                        PgSAHelpers::writeValue<uint_read_len_max>(*rlMisPosDest, rlEntry.mismatchOffset[i]);
                }
            }
        }
        readsCounter++;
    }

    void SeparatedPseudoGenomeOutputBuilder::setReadsSourceIterator(DefaultReadsListIteratorInterface *rlIt) {
        this->rlIt = rlIt;
    }

    void SeparatedPseudoGenomeOutputBuilder::writeReadsFromIterator(uint_pg_len_max stopPos) {
        if (iterationPaused) {
            if (rlIt->peekReadEntry().pos >= stopPos)
                return;
            writeReadEntry(rlIt->peekReadEntry());
            iterationPaused = false;
        }
        while (rlIt->moveNext()) {
            if (rlIt->peekReadEntry().pos >= stopPos) {
                iterationPaused = true;
                break;
            }
            writeReadEntry(rlIt->peekReadEntry());
        }
    }

    void SeparatedPseudoGenomeOutputBuilder::writePseudoGenome(PseudoGenomeBase *pgb, string divisionFile,
            bool divisionComplement, bool revComplPairFile) {

        initDest(pgDest, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX);
        (*pgDest) << pgb->getPseudoGenomeVirtual();

        ReadsListIteratorExtendedWrapperBase* rlIt =
                TemplateUserGenerator::generateReadsListUser<ReadsListIteratorExtendedWrapper, ReadsListIteratorExtendedWrapperBase>(pgb);
        if (divisionFile != "")
            rlIt->applyDivision(divisionFile, divisionComplement);
        if (revComplPairFile)
            rlIt->applyRevComplPairFileFlag();

        setReadsSourceIterator(rlIt);
        writeReadsFromIterator();

        pgh = new PseudoGenomeHeader(pgb);
        if (pgh->getReadsCount() != readsCounter) {
            fprintf(stderr, "Incorrect reads count validation while building separated Pg (%llu instead of %llu).\n",
                    readsCounter, pgh->getReadsCount());
            exit(EXIT_FAILURE);
        }

        delete(rlIt);
    }

    void SeparatedPseudoGenomeOutputBuilder::copyPseudoGenomeHeader(const string &pseudoGenomePrefix) {
        ifstream pgPropSrc(pseudoGenomePrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX,
                           ios_base::in | ios_base::binary);
        if (pgPropSrc.fail()) {
            fprintf(stderr, "Cannot read pseudogenome properties (%s does not open).\n",
                    pseudoGenomePrefix.c_str());
            exit(EXIT_FAILURE);
        }
        this->pgh = new PseudoGenomeHeader(pgPropSrc);
        pgPropSrc.close();
    }

}