#ifndef PGTOOLS_DEFAULTPGMATCHER_H
#define PGTOOLS_DEFAULTPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"
#include "../pseudogenome/readslist/SeparatedExtendedReadsList.h"
#include "TextMatchers.h"

namespace PgTools {

    using namespace PgIndex;

    struct PgMatch;

    class DefaultPgMatcher {
    private:

        const string srcPgPrefix;
        const string targetPgPrefix;
        const bool revComplMatching;

        bool destPgIsSrcPg;

        bool plainTextReadMode = false;

        PseudoGenomeHeader* srcPgh = nullptr;
        ReadsSetProperties* srcRsProp = nullptr;
        uint_read_len_max readLength;
        string srcPg;
        string destPg;

        vector<TextMatch> textMatches;
        vector<PgMatch> pgMatches;

        ExtendedReadsListWithConstantAccessOption* srcRl = nullptr;
        vector<uint_pg_len_max> newRlPos;

        void mapPgMatches2SrcReadsList();
        void reverseDestWithSrcForBetterMatchesMappingInTheSamePg();
        void correctDestPositionDueToRevComplMatching();
        void resolveDestSrcReadsOverlapConflictsInTheSamePg();
        bool resolveDestSrcCollision(PgMatch &destMatch, PgMatch &srcMatch, uint_pg_len_max &collidedCharsCount);
        void resolveMatchesOverlapInSrc();
        bool resolveSrcSrcCollision(PgMatch &leftMatch, PgMatch &rightMatch, uint_pg_len_max &collidedCharsCount);

        void transferMatchesFromSrcToDest(const string &destPgPrefix);

        void fillPgMatches();

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        DefaultPgMatcher(const string& srcPgPrefix, const string& targetPgPrefix, bool revComplMatching);

        void exactMatchPg(uint32_t minMatchLength);

        virtual ~DefaultPgMatcher();

        void writeMatchesInfo(const string &dumpFilePrefix);

        void transferMatchedReads(const string &destPgFilePrefix);

        static void matchPgInPgFile(const string &srcPgPrefix, const string &targetPgPrefix, uint_pg_len_max targetMatchLength,
                             const string &destPgPrefix, bool revComplPg, bool dumpInfo);

    };

    struct PgMatch {
        TextMatch mapping;

        uint_reads_cnt_max startRlIdx = -1;
        uint_reads_cnt_max endRlIdx = -1;

        uint_pg_len_max leftSrcLockedMargin = 0;
        uint_pg_len_max rightSrcLockedMargin = 0;

        PgMatch(TextMatch& match): mapping(match) {}

        void alignToReads(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength, int64_t &i) {
            while(--i > 0 && rlPos[i] >= mapping.posSrcText);
            while(rlPos[++i] < mapping.posSrcText);
            startRlIdx = i--;
            while(rlPos[++i] + readLength <= mapping.posSrcText + mapping.length);
            endRlIdx = i - 1;
        }

        uint_pg_len_max netPosSrcPg(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) const {
            return leftSrcLockedMargin + (startRlIdx > 0?rlPos[startRlIdx - 1] + readLength:0);
        }

        uint_pg_len_max netSrcPosAlignment(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) const {
            return netPosSrcPg(rlPos, readLength) - mapping.posSrcText;
        };

        uint_pg_len_max netEndPosSrcPg(const vector<uint_pg_len_max> &rlPos) const {
            return rlPos[endRlIdx + 1] - rightSrcLockedMargin;
        }

        uint_pg_len_max netSrcLength(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) const {
            int64_t len = (int64_t) netEndPosSrcPg(rlPos) - netPosSrcPg(rlPos, readLength);
            return len > 0?len:0;
        };

        bool inactive(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) {
            return netSrcLength(rlPos, readLength) == 0;
        };

        void trimLeft(uint_pg_len_max overlapLength,
                      const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) {
            const uint_pg_len_max leftPos = leftSrcLockedMargin + (startRlIdx > 0 ? rlPos[startRlIdx - 1] : 0);
            uint_reads_cnt_max rsIdx = startRlIdx - 1;
            while (++rsIdx <= endRlIdx && rlPos[rsIdx] - leftPos <= overlapLength);
            leftSrcLockedMargin = overlapLength - (rlPos[rsIdx - 1] - leftPos);
            startRlIdx = rsIdx;
        }

        void trimRight(uint_pg_len_max overlapLength,
                const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) {
            const uint_pg_len_max rightPos = rlPos[endRlIdx + 1] - rightSrcLockedMargin;
            uint_reads_cnt_max leIdx = endRlIdx + 1;
            while (--leIdx >= startRlIdx && rightPos - rlPos[leIdx] <= overlapLength);
            rightSrcLockedMargin = overlapLength - (rightPos - rlPos[leIdx + 1]);
            endRlIdx = leIdx;
        }

        uint_pg_len_max mapSrcReadToDest(const uint_pg_len_max &srcReadPos, const uint_read_len_max &readLength,
                                         bool revComplMatch) const {
            uint_pg_len_max srcOffset = srcReadPos - mapping.posSrcText;
            return mapping.posDestText + (revComplMatch?mapping.length - srcOffset - readLength:srcOffset);
        }

        uint_pg_len_max netPosDestPg(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength,
                bool revComplMatch) const {
            return revComplMatch?mapSrcReadToDest(rlPos[endRlIdx], readLength, true)
                :mapSrcReadToDest(rlPos[startRlIdx], readLength, false);
        }

        uint_pg_len_max netEndPosDestPg(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength,
                bool revComplMatch) const {
            return readLength + (revComplMatch?mapSrcReadToDest(rlPos[startRlIdx], readLength, true)
                :mapSrcReadToDest(rlPos[endRlIdx], readLength, false));
        }

        uint_pg_len_max netDestLength(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength,
                bool revComplMatch) const {
            return netEndPosDestPg(rlPos, readLength, revComplMatch) + netPosDestPg(rlPos, readLength,revComplMatch);
        };

        void report(ostream& out) {
            out << mapping.length << ": <" << mapping.posSrcText << ", " << mapping.endPosSrcText() << ") in " << mapping.posDestText << endl;
        }
    };

}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
