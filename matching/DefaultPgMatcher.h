#ifndef PGTOOLS_DEFAULTPGMATCHER_H
#define PGTOOLS_DEFAULTPGMATCHER_H

#include "../pseudogenome/PseudoGenomeBase.h"

namespace PgTools {

    using namespace PgSAIndex;

    struct PgMatch;

    class DefaultPgMatcher {
    private:
        const string srcPgPrefix;
        const string targetPgPrefix;
        const bool revComplMatching;

        bool textFromSamePg;

        bool plainTextReadMode = false;

        PseudoGenomeHeader* srcPgh = 0;
        string srcPg;
        string destPg;

        vector<PgMatch> pgMatches;

        vector<uint_pg_len_max> rlPos;
        vector<uint_reads_cnt_max> rlIdx;
        void fillSrcReadsList();

        void mapPgMatches2SrcReadsList();
        void reverseDestWithSrcForBetterMatchesMappingInTheSamePg();
        void correctDestPositionDueToRevComplMatching();
        void resolveDestOverlapSrcConflictsInTheSamePg();
        bool resolveCollision(PgMatch &destMatch, PgMatch &srcMatch, uint_pg_len_max& collidedCharsCount);

        string getTotalMatchStat(uint_pg_len_max totalMatchLength);

    public:
        DefaultPgMatcher(const string& srcPgPrefix, const string& targetPgPrefix, bool revComplMatching);

        void exactMatchPg(uint32_t minMatchLength);

        virtual ~DefaultPgMatcher();

        void writeMatchesInfo(const string &dumpFilePrefix);

        void writeIntoPseudoGenome(const string &destPgFilePrefix);
    };


    struct PgMatch{
        uint_pg_len_max posSrcPg;
        uint_pg_len_max length;
        uint_pg_len_max posDestPg;

        uint_read_len_max netPosAlignment = 0;
        uint_pg_len_max netLength = 0;
        uint_reads_cnt_max startRlIdx = -1;
        uint_reads_cnt_max endRlIdx = -1;

        bool inactive = false;

        PgMatch(uint_pg_len_max posSrcPg, uint_pg_len_max length, uint_pg_len_max posDestPg) : posSrcPg(posSrcPg), length(length),
                                                                                          posDestPg(posDestPg) {}
        void reverseMatch() {
            uint_pg_len_max temp = posSrcPg;
            posSrcPg = posDestPg;
            posDestPg = temp;
        }

        uint_pg_len_max endPosSrcPg() const {
            return posSrcPg + length;
        }

        uint_pg_len_max endPosDestPg() const {
            return posDestPg + length;
        }

        uint_pg_len_max netPosSrcPg() const {
            return posSrcPg + netPosAlignment;
        }

        uint_pg_len_max netPosDestPg() const {
            return posDestPg + netPosAlignment;
        }

        uint_pg_len_max netEndPosSrcPg() const {
            return netPosSrcPg() + netLength;
        }

        uint_pg_len_max netEndPosDestPg() const {
            return netPosDestPg() + netLength;
        }

        void report(ostream& out) {
            out << length << ": <" << posSrcPg << ", " << endPosSrcPg() << ") in " << posDestPg << endl;
        }
    };
}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
