#ifndef PGTOOLS_SIMPLEPGMATCHER_H
#define PGTOOLS_SIMPLEPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"
#include "../pseudogenome/persistence/SeparatedExtendedReadsList.h"
#include "TextMatchers.h"

namespace PgTools {

    using namespace PgSAIndex;

    struct PgMatch;

    class SimplePgMatcher {
    private:

        const string srcPgPrefix;
        uint32_t targetMatchLength;
        uint32_t minMatchLength;

        TextMatcher* matcher;
        string srcPg;

        string targetPgPrefix;
        string destPg;
        vector<TextMatch> textMatches;
        bool revComplMatching;

        void exactMatchPg();

        void correctDestPositionDueToRevComplMatching();
        void resolveMappingCollisionsInTheSameText();

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        SimplePgMatcher(const string& srcPgPrefix, uint32_t targetMatchLength);

        virtual ~SimplePgMatcher();

        void markAndRemoveExactMatches(const string &destPgPrefix, bool revComplMatching);

        static void matchPgInPgFiles(const string &goodPgPrefix, const string &badPgPrefix, uint_pg_len_max targetMatchLength,
                             bool revComplMatching);

    };
}

#endif //PGTOOLS_SIMPLEPGMATCHER_H
