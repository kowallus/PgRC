#ifndef PGTOOLS_PGMATCHERBASE_H
#define PGTOOLS_PGINPGMATCHERBASE_H

namespace PgTools {

    class PgMatcherBase {
    public:
        virtual ~PgMatcherBase() {}

        virtual void exactMatchPg(string text, ofstream &offsetsDest, uint32_t minMatchLength, bool textFromSamePg) = 0;
    };

}

#endif //PGTOOLS_PGMATCHERBASE_H
