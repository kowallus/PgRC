#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"

#include <cstdio>

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, const string &pseudoGenomePrefix, string divisionFile, bool divisionComplement) {
        clock_checkpoint();

        SeparatedPseudoGenomeOutputBuilder builder(pseudoGenomePrefix);
        std::ofstream* destPgProp = builder.getPseudoGenomeElementDest(PSEUDOGENOME_PROPERTIES_SUFFIX);
        PgSAIndex::PseudoGenomePersistence::writePseudoGenomeHeader(pgb, *destPgProp);

        std::ofstream* destPg = builder.getPseudoGenomeElementDest(PSEUDOGENOME_FILE_SUFFIX);
        (*destPg) << pgb->getPseudoGenomeVirtual();
        destPg->close();

        SeparatedReadsListWriterBase* srlwb = TemplateUserGenerator::generateReadsListUser<SeparatedReadsListWriter, SeparatedReadsListWriterBase>(pgb);
        srlwb->writeReadsList(builder.getPseudoGenomeElementDest(READSLIST_OFFSETS_FILE_SUFFIX),
                builder.getPseudoGenomeElementDest(READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX), divisionFile, divisionComplement);

        builder.build();
        cout << "Writing pseudo genome files in " << clock_millis() << " msec." << endl;
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
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_OFFSETS_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_REVERSECOMPL_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHESCOUNT_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHEDSYMBOLS_FILE_SUFFIX);
        acceptTemporaryPseudoGenomeElement(pseudoGenomePrefix, READSLIST_MISMATCHESPOS_FILE_SUFFIX);
    }

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX = ".pg";
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX = ".pg.prop";
    const string SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX = ".pg.rl.off";
    const string SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = ".pg.rl.idx";
    const string SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX = ".pg.rl.rc";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHESCOUNT_FILE_SUFFIX = ".pg.rl.mis.cnt";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHEDSYMBOLS_FILE_SUFFIX = ".pg.rl.mis.sym";
    const string SeparatedPseudoGenomePersistence::READSLIST_MISMATCHESPOS_FILE_SUFFIX = ".pg.rl.mis.pos";

    const string SeparatedPseudoGenomePersistence::TEMPORARY_FILE_SUFFIX = ".temp";

    SeparatedPseudoGenomeOutputBuilder::SeparatedPseudoGenomeOutputBuilder(const string &pseudoGenomePrefix)
            : pseudoGenomePrefix(pseudoGenomePrefix) {}

    std::ofstream* SeparatedPseudoGenomeOutputBuilder::getSingletonDest(ofstream *&dest, const string &fileSuffix) {
        if (dest == 0) {
            dest = new ofstream(SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(pseudoGenomePrefix, fileSuffix, true));
        }
        return dest;
    }

    std::ofstream* SeparatedPseudoGenomeOutputBuilder::getPseudoGenomeElementDest(const string &fileSuffix) {
        if (fileSuffix == SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX)
            return getSingletonDest(pgDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX)
            return getSingletonDest(pgPropDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX)
            return getSingletonDest(rlOffDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX)
            return getSingletonDest(rlOrgIdxDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::READSLIST_REVERSECOMPL_FILE_SUFFIX)
            return getSingletonDest(rlRevCompDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::READSLIST_MISMATCHESCOUNT_FILE_SUFFIX)
            return getSingletonDest(rlMisCntDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::READSLIST_MISMATCHEDSYMBOLS_FILE_SUFFIX)
            return getSingletonDest(rlMisSymDest, fileSuffix);
        else if (fileSuffix == SeparatedPseudoGenomePersistence::READSLIST_MISMATCHESPOS_FILE_SUFFIX)
            return getSingletonDest(rlMisPosDest, fileSuffix);

        fprintf(stderr, "Unsupported pseudogenome element file suffix: %s\n",
                fileSuffix.c_str());
        exit(EXIT_FAILURE);

        return 0;
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
        freeDest(rlOffDest);
        freeDest(rlOrgIdxDest);
        freeDest(rlRevCompDest);
        freeDest(rlMisCntDest);
        freeDest(rlMisSymDest);
        freeDest(rlMisPosDest);
    }

    void SeparatedPseudoGenomeOutputBuilder::build() {
        freeDests();
        SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(pseudoGenomePrefix);
    }

}