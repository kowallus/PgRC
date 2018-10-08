#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, const string &pseudoGenomePrefix, string divisionFile, bool divisionComplement) {
        clock_checkpoint();

        std::ofstream destPgProp(pseudoGenomePrefix + PSEUDOGENOME_PROPERTIES_SUFFIX, std::ios::out | std::ios::binary);
        PgSAIndex::PseudoGenomePersistence::writePseudoGenomeHeader(pgb, destPgProp);
        destPgProp.close();
        std::ofstream destPg(pseudoGenomePrefix + PSEUDOGENOME_FILE_SUFFIX, std::ios::out | std::ios::binary);
        destPg << pgb->getPseudoGenomeVirtual();
        destPg.close();

        SeparatedReadsListWriterBase* srlwb = TemplateUserGenerator::generateReadsListUser<SeparatedReadsListWriter, SeparatedReadsListWriterBase>(pgb);
        srlwb->writeReadsList(pseudoGenomePrefix, divisionFile, divisionComplement);
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

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX = ".pg";
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX = ".pgprop" ;
    const string SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX = ".rloff" ;
    const string SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = ".rlidx" ;

}