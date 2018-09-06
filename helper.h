#ifndef PGTOOLS_HELPER_H
#define PGTOOLS_HELPER_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

namespace pgTools {

    static const string PSEUDOGENOME_HEADER = "PGEN";

    string getPgFromPgenFile(ifstream &pgSrc);

};


#endif //PGTOOLS_HELPER_H
