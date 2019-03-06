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
        uint32_t targetMatchLength;

        TextMatcher* matcher;
        const string& srcPg;

        uint64_t destPgLength;
        vector<TextMatch> textMatches;
        bool revComplMatching;

        void exactMatchPg(string& destPg, bool destPgIsSrcPg, uint32_t minMatchLength);

        void correctDestPositionDueToRevComplMatching();
        void resolveMappingCollisionsInTheSameText();

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        SimplePgMatcher(const string& srcPg, uint32_t targetMatchLength,
                uint32_t minMatchLength = UINT32_MAX);

        virtual ~SimplePgMatcher();

        void markAndRemoveExactMatches(bool destPgIsSrcPg,
                string &destPg, bool revComplMatching, uint32_t minMatchLength = UINT32_MAX);

        string pgMapped;
        string pgMapOff;
        string pgMapLen;

        static void matchPgsInPg(string &hqPgSequence, string &lqPgSequence, ostream &pgrcOut, uint8_t coder_level,
                                 const string &hqPgPrefix, const string &lqPgPrefix, uint_pg_len_max targetMatchLength,
                                 uint32_t minMatchLength = UINT32_MAX);

        static void restoreMatchedPgs(istream &pgrcIn, string &hqPgSequence, string &lqPgSequence);

        static string restoreMatchedPg(string &srcPg, size_t srcLen, const string& destPg,
                istream &pgMapOffSrc, istream &pgMapLenSrc,
                bool revComplMatching, bool plainTextReadMode, bool srcIsDest = false);

        static string restoreMatchedPg(string &srcPg, const string& destPgPrefix,
                                       bool revComplMatching, bool plainTextReadMode);

        static string restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching);

        void writeMatchingResult(const string &pgPrefix);
    };
}

#endif //PGTOOLS_SIMPLEPGMATCHER_H
