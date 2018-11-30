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

        uint_pg_len_max nettoLength = 0;
        uint_reads_cnt_max startRlIdx = -1;

        PgMatch(uint_pg_len_max posSrcPg, uint_pg_len_max length, uint_pg_len_max posDestPg) : posSrcPg(posSrcPg), length(length),
                                                                                          posDestPg(posDestPg) {}
        void reverseMatch() {
            uint_pg_len_max temp = posSrcPg;
            posSrcPg = posDestPg;
            posDestPg = temp;
        }

        uint_pg_len_max endPosSrcPg() {
            return posSrcPg + length;
        }

        uint_pg_len_max endPosDestPg() {
            return posDestPg + length;
        }

        void report(ostream& out) {
            out << length << ": <" << posSrcPg << ", " << endPosSrcPg() << ") in " << posDestPg << endl;
        }
    };
}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
