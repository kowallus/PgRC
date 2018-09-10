#ifndef PGTOOLS_PGHELPER_H
#define PGTOOLS_PGHELPER_H

#include <iostream>
#include <fstream>
#include <string>

#include "pseudogenome/persistence/PseudoGenomePersistence.h"

using namespace std;

namespace PgTools {

    PseudoGenomeBase* openPg(string pgFile);

};


#endif //PGTOOLS_PGHELPER_H
