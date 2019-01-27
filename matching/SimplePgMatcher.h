#ifndef PGTOOLS_SIMPLEPGMATCHER_H
#define PGTOOLS_SIMPLEPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"
#include "../pseudogenome/readslist/SeparatedExtendedReadsList.h"
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
        const string& srcPg;

        string targetPgPrefix;
        string destPg;
        vector<TextMatch> textMatches;
        bool revComplMatching;

        void exactMatchPg();

        void correctDestPositionDueToRevComplMatching();
        void resolveMappingCollisionsInTheSameText();

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        SimplePgMatcher(const string& srcPgPrefix, const string& srcPg, uint32_t targetMatchLength);

        virtual ~SimplePgMatcher();

        void markAndRemoveExactMatches(const string &destPgPrefix, const string &queryPg, bool revComplMatching);

        static void matchPgInPgFiles(string& hqPgSequence, string& lqPgSequence,
                const string &hqPgPrefix, const string &lqPgPrefix, uint_pg_len_max targetMatchLength);

    };
}

#endif //PGTOOLS_SIMPLEPGMATCHER_H
