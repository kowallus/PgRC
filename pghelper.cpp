#include "pghelper.h"

namespace PgTools {

    PseudoGenomeBase* openPg(string pgFile) {
        PseudoGenomeBase* pgb = 0;

        if (PseudoGenomePersistence::isValidPseudoGenome(pgFile))
            pgb = PseudoGenomePersistence::readPseudoGenome(pgFile);

        if (pgb == 0) {
            fprintf(stderr, "Failed loading Pg\n");
            exit(EXIT_FAILURE);
        }

        return pgb;
    }

}