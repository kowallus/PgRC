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

        TextMatcher* matcher = 0;
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
                string &destPg, string &resPgMapOff, string& resPgMapLen,
                bool revComplMatching, uint32_t minMatchLength = UINT32_MAX);

        static void matchPgsInPg(string &hqPgSequence, string &lqPgSequence, string &nPgSequence,
                                    bool separateNReads, ostream &pgrcOut, uint8_t coder_level,
                                    const string &hqPgPrefix, const string &lqPgPrefix, const string &nPgPrefix,
                                    uint_pg_len_max targetMatchLength, uint32_t minMatchLength = UINT32_MAX);

        static void restoreMatchedPgs(istream &pgrcIn, uint_pg_len_max orgHqPgLen,
                string &hqPgSequence, string &lqPgSequence, string &nPgSequence);

        static string restoreMatchedPg(string &srcPg, size_t orgSrcLen, const string& destPg,
                istream &pgMapOffSrc, istream &pgMapLenSrc,
                bool revComplMatching, bool plainTextReadMode, bool srcIsDest = false);

        static string restoreMatchedPg(string &srcPg, const string& destPgPrefix,
                                       bool revComplMatching, bool plainTextReadMode);

        static string restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching);

        static void writeMatchingResult(const string &pgPrefix,
                const string &pgMapped, const string &pgMapOff, const string& pgMapLen);

        static char MATCH_MARK;
    };
}

#endif //PGTOOLS_SIMPLEPGMATCHER_H
