#ifndef PGTOOLS_DEFAULTPGMATCHER_H
#define PGTOOLS_DEFAULTPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"

namespace PgTools {

    using namespace PgSAIndex;

    struct PgMatch;

    void matchPgInPgFile(const string &srcPgPrefix, const string &targetPgPrefix, uint_pg_len_max targetMatchLength,
                         const string &destPgPrefix, bool revComplPg, bool dumpInfo);

    class DefaultPgMatcher {
    private:

        const string srcPgPrefix;
        const string targetPgPrefix;
        const bool revComplMatching;

        bool destPgIsSrcPg;

        bool plainTextReadMode = false;

        PseudoGenomeHeader* srcPgh = 0;
        uint_read_len_max readLength;
        string srcPg;
        string destPg;


        vector<PgMatch> pgMatches;

        vector<uint_pg_len_max> rlPos;
        vector<uint_pg_len_max> newRlPos;
        vector<uint_reads_cnt_max> rlIdx;
        void fillSrcReadsList();

        void mapPgMatches2SrcReadsList();
        void reverseDestWithSrcForBetterMatchesMappingInTheSamePg();
        void correctDestPositionDueToRevComplMatching();
        void resolveDestOverlapSrcConflictsInTheSamePg();
        bool resolveDestSrcCollision(PgMatch &destMatch, PgMatch &srcMatch, uint_pg_len_max &collidedCharsCount);
        void resolveMatchesOverlapInSrc();
        bool resolveSrcSrcCollision(PgMatch &leftMatch, PgMatch &rightMatch, uint_pg_len_max &collidedCharsCount);

        void transferMatchesFromSrcToDest(const string &destPgPrefix);

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        DefaultPgMatcher(const string& srcPgPrefix, const string& targetPgPrefix, bool revComplMatching);

        void exactMatchPg(uint32_t minMatchLength);

        virtual ~DefaultPgMatcher();

        void writeMatchesInfo(const string &dumpFilePrefix);

        void writeIntoPseudoGenome(const string &destPgFilePrefix);

        int64_t trimOverlap(PgMatch &leftMatch, PgMatch &rightMatch, uint_pg_len_max overlapLength);
    };

    struct PgMatch {
        uint_pg_len_max posGrossSrcPg;
        uint_pg_len_max grossLength;
        uint_pg_len_max posGrossDestPg;

        uint_reads_cnt_max startRlIdx = -1;
        uint_reads_cnt_max endRlIdx = -1;

        uint_pg_len_max leftSrcLockedMargin = 0;
        uint_pg_len_max rightSrcLockedMargin = 0;

        PgMatch(uint_pg_len_max posSrcPg, uint_pg_len_max length, uint_pg_len_max posDestPg) :
            posGrossSrcPg(posSrcPg), grossLength(length), posGrossDestPg(posDestPg) {}

        void alignToReads(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength, int64_t &i) {
            while(--i > 0 && rlPos[i] >= posGrossSrcPg);
            while(rlPos[++i] < posGrossSrcPg);
            startRlIdx = i--;
            while(rlPos[++i] + readLength <= posGrossSrcPg + grossLength);
            endRlIdx = i - 1;
        }

        uint_pg_len_max endGrossPosSrcPg() const {
            return posGrossSrcPg + grossLength;
        }

        uint_pg_len_max endGrossPosDestPg() const {
            return posGrossDestPg + grossLength;
        }

        uint_pg_len_max netPosSrcPg(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) const {
            return leftSrcLockedMargin + (startRlIdx > 0?rlPos[startRlIdx - 1] + readLength:0);
        }

        uint_pg_len_max netSrcPosAlignment(const vector<uint_pg_len_max> &rlPos, const uint_read_len_max &readLength) const {
            return netPosSrcPg(rlPos, readLength) - posGrossSrcPg;
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
            uint_pg_len_max srcOffset = srcReadPos - posGrossSrcPg;
            return posGrossDestPg + (revComplMatch?grossLength - srcOffset - readLength:srcOffset);
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
            out << grossLength << ": <" << posGrossSrcPg << ", " << endGrossPosSrcPg() << ") in " << posGrossDestPg << endl;
        }
    };

}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
