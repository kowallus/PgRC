#include "SeparatedPseudoGenomePersistence.h"
#include "PseudoGenomePersistence.h"
#include "../TemplateUserGenerator.h"

namespace PgTools {

    void SeparatedPseudoGenomePersistence::writePseudoGenome(PseudoGenomeBase *pgb, string pseudoGenomePrefix, string divisionFile, bool divisionComplement) {
        clock_checkpoint();

        std::ofstream destPgProp(pseudoGenomePrefix + PSEUDOGENOME_PROPERTIES_SUFFIX, std::ios::out | std::ios::binary);
        PgSAIndex::PseudoGenomePersistence::writePseudoGenomeHeader(pgb, destPgProp);
        destPgProp.close();
        std::ofstream destPg(pseudoGenomePrefix + PSEUDOGENOME_FILE_SUFFIX, std::ios::out | std::ios::binary);
        destPg << pgb->getPseudoGenomeVirtual();
        destPg.close();
        SeparatedReadsListWriterBase* srlwb = TemplateUserGenerator::generateReadsListUser<SeparatedReadsListWriter, SeparatedReadsListWriterBase>(pgb);
        const uint_reads_cnt_max readsCount = pgb->getReadsSetProperties()->readsCount;

        vector<uint_reads_cnt_max> orgIndexesMapping(readsCount);
        if (divisionFile == "") { 
            for(uint64_t i = 0; i < readsCount; i++)
                orgIndexesMapping[i] = i;
        } else {
            std::ifstream divSource(divisionFile, std::ios::in | std::ios::binary);
            if (divSource.fail()) {
                fprintf(stderr, "cannot open reads division file %s\n", divisionFile.c_str());
                exit(EXIT_FAILURE);
            }
            uint64_t i = 0;
            uint64_t counter = 0;
            uint64_t currentDivIdx = 0;
            readValue(divSource, currentDivIdx);
            while (i < readsCount) {
                if (divisionComplement) {
                    while (counter == currentDivIdx) {
                        readValue(divSource, currentDivIdx);
                        counter++;
                    }
                    orgIndexesMapping[i++] = counter++;
                } else {
                    while (counter != currentDivIdx)
                        counter++;
                    orgIndexesMapping[i++] = counter++;
                    readValue(divSource, currentDivIdx);
                }
            }
        }

        srlwb->writeReadsList(pseudoGenomePrefix, orgIndexesMapping);
        cout << "Writing pseudo genome files in " << clock_millis() << " msec." << endl;
    }

    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX = ".pg";
    const string SeparatedPseudoGenomePersistence::PSEUDOGENOME_PROPERTIES_SUFFIX = ".pgprop" ;
    const string SeparatedPseudoGenomePersistence::READSLIST_OFFSETS_FILE_SUFFIX = ".rloff" ;
    const string SeparatedPseudoGenomePersistence::READSLIST_ORIGINAL_INDEXES_FILE_SUFFIX = ".rlidx" ;

}