#ifndef PSEUDOGENOMEPERSISTENCE_H
#define PSEUDOGENOMEPERSISTENCE_H

#include "../DefaultPseudoGenome.h"
#include "../PackedPseudoGenome.h"
#include "../../utils/helper.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    class PseudoGenomePersistence
    {
        private:
            PseudoGenomePersistence();;

        public:
            virtual ~PseudoGenomePersistence();;

            static void writePseudoGenome(PseudoGenomeBase* pgb, std::ostream& dest);

            static void writePseudoGenome(PseudoGenomeBase* pgb, string pseudoGenomePrefix);

            static bool isValidPseudoGenome(string file);
            
            template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len >
            static PseudoGenomeBase* readPseudoGenomeTemplate(std::istream& src, PseudoGenomeHeader& pgh);

            static PseudoGenomeBase* readPseudoGenome(std::istream& src);

            static PseudoGenomeBase* readPseudoGenome(string pseudoGenomeFile);

            static void writePseudoGenomeHeader(PseudoGenomeBase *pgb, std::ostream &dest);

            static PseudoGenomeBase* checkAndReadPseudoGenome(string pgFile);
    };

}

#endif // PSEUDOGENOMEPERSISTENCE_H
