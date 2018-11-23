#ifndef PGTOOLS_DEFAULTPGMATCHER_H
#define PGTOOLS_DEFAULTPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"

namespace PgTools {

    using namespace PgSAIndex;

    class DefaultPgMatcher {
    private:
        const string srcPgPrefix;
        PseudoGenomeHeader* pgh = 0;
        bool plainTextReadMode = false;

    public:
        DefaultPgMatcher(const string& srcPgPrefix);

        void exactMatchPg(const string& text, ofstream &offsetsDest, uint32_t minMatchLength,
                bool textFromSamePg, bool textIsRevComplOfPg);

        virtual ~DefaultPgMatcher();
    };

}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
