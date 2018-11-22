#ifndef PGTOOLS_DEFAULTPGMATCHER_H
#define PGTOOLS_DEFAULTPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"
#include "PgMatcherBase.h"

namespace PgTools {

    using namespace PgSAIndex;

    class DefaultPgMatcher: public PgMatcherBase {
    private:
        const string srcPgPrefix;
        PseudoGenomeHeader* pgh = 0;
        bool plainTextReadMode = false;

    public:
        DefaultPgMatcher(const string& srcPgPrefix);

        void exactMatchPg(const string& text, ofstream &offsetsDest, uint32_t minMatchLength, bool textFromSamePg) override;

        virtual ~DefaultPgMatcher();
    };

}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
