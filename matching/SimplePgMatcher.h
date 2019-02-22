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

        TextMatcher* matcher;
        const string& srcPg;

        string targetPgPrefix;
        uint64_t destPgLength;
        vector<TextMatch> textMatches;
        bool revComplMatching;

        void exactMatchPg(string& destPg, uint32_t minMatchLength);

        void correctDestPositionDueToRevComplMatching();
        void resolveMappingCollisionsInTheSameText();

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        SimplePgMatcher(const string& srcPgPrefix, const string& srcPg, uint32_t targetMatchLength,
                uint32_t minMatchLength = UINT32_MAX);

        virtual ~SimplePgMatcher();

        void markAndRemoveExactMatches(const string &destPgPrefix, string &destPg, bool revComplMatching,
                uint32_t minMatchLength = UINT32_MAX);

        static void matchPgInPgFiles(string& hqPgSequence, string& lqPgSequence,
                const string &hqPgPrefix, const string &lqPgPrefix, uint_pg_len_max targetMatchLength,
                uint32_t minMatchLength = UINT32_MAX);

        static string restoreMatchedPg(string &srcPg, const string &destPgPrefix, bool revComplMatching,
                                       bool plainTextReadMode, bool srcIsDest = false);

        static string restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching);
    };
}

#endif //PGTOOLS_SIMPLEPGMATCHER_H
