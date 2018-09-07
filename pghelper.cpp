#include "pghelper.h"

namespace pgTools {

    static const int PSEUDOGENOME_PRECEDING_LINES = 13;

    string getPgFromPgenFile(ifstream &pgSrc) {
        string line;
        pgSrc >> line;
        if (line == PSEUDOGENOME_HEADER) {
            for (int i = 0; i < PSEUDOGENOME_PRECEDING_LINES; i++)
                pgSrc >> line;
        }
        return line;
    }

}