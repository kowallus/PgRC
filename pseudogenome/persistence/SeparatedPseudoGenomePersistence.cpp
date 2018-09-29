#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, string pseudoGenomePrefix) {
        std::ofstream destPgProp(pseudoGenomePrefix + PSEUDOGENOME_PROPERTIES_SUFFIX, std::ios::out | std::ios::binary);
        PgSAIndex::PseudoGenomePersistence::writePseudoGenomeHeader(pgb, destPgProp);
        destPgProp.close();
        std::ofstream destPg(pseudoGenomePrefix + PSEUDOGENOME_FILE_SUFFIX, std::ios::out | std::ios::binary);
        destPg << pgb->getPseudoGenomeVirtual();
        destPg.close();
        SeparatedReadsListWriterBase* srlwb = TemplateUserGenerator::generateReadsListUser<SeparatedReadsListWriter, SeparatedReadsListWriterBase>(pgb);

        srlwb->writeReadsList(pseudoGenomePrefix);
    }

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX = ".pg";
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX = ".pgprop" ;
    const string SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX = ".rloff" ;
    const string SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = ".rlidx" ;

}