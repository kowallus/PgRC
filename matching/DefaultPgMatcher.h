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
        uint_pg_len_max posPg;
        uint_pg_len_max length;
        uint_pg_len_max posText;

        PgMatch(uint_pg_len_max posPg, uint_pg_len_max length, uint_pg_len_max posText) : posPg(posPg), length(length),
                                                                                          posText(posText) {}

        uint_pg_len_max endPosPg() {
            return posPg + length;
        }

        uint_pg_len_max endPosText() {
            return posText + length;
        }

        void report(ostream& out) {
            out << length << ": <" << posPg << ", " << endPosPg() << ") in " << posText << endl;
        }
    };
}


#endif //PGTOOLS_DEFAULTPGMATCHER_H
