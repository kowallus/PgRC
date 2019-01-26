#include "SimplePgMatcher.h"

#include "copmem/CopMEMMatcher.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    SimplePgMatcher::SimplePgMatcher(const string& srcPgPrefix, const string& srcPg, uint32_t minMatchLength)
            :srcPgPrefix(srcPgPrefix), srcPg(srcPg), targetMatchLength(minMatchLength), minMatchLength(minMatchLength) {
        cout << "Source pseudogenome length: " << srcPg.length() << endl;

        //matcher = new DefaultTextMatcher(srcPg, minMatchLength);
        matcher = new CopMEMMatcher(srcPg, minMatchLength);

    }

    SimplePgMatcher::~SimplePgMatcher() {
        delete(matcher);
    }

    void SimplePgMatcher::exactMatchPg() {
        clock_checkpoint();
        bool destPgIsSrcPg = srcPgPrefix == targetPgPrefix;
        if (!destPgIsSrcPg)
            cout << "Destination pseudogenome length: " << destPg.length() << endl;
        matcher->matchTexts(textMatches, destPg, destPgIsSrcPg, revComplMatching, targetMatchLength);
        cout << "... found " << textMatches.size() << " exact matchers in " << clock_millis() << " msec. " << endl;
/*        std::sort(textMatches.begin(), textMatches.end(), [](const TextMatch &match1, const TextMatch &match2) -> bool
            { return match1.length > match2.length; });
        cout << "Largest matches:" << endl;
        for (uint32_t i = 0; i < textMatches.size() && i < 10; i++)
            textMatches[i].report(cout);
*/
        if (revComplMatching) {
            correctDestPositionDueToRevComplMatching();
            if (destPgIsSrcPg)
                destPg = srcPg;
            else
                PgSAHelpers::reverseComplementInPlace(destPg);
        }
    }

    using namespace PgTools;

    void SimplePgMatcher::correctDestPositionDueToRevComplMatching() {
        for (TextMatch& match: textMatches)
            match.posDestText = destPg.length() - (match.posDestText + match.length);
    }

    string SimplePgMatcher::getTotalMatchStat(uint_pg_len_max totalMatchLength) {
        return toString(totalMatchLength) + " (" + toString((totalMatchLength * 100.0) / destPg.length(), 1)+ "%)";
    }

    void SimplePgMatcher::markAndRemoveExactMatches(const string &destPgFilePrefix, const string &destPg, bool revComplMatching) {
        this->targetPgPrefix = destPgFilePrefix;
        this->revComplMatching = revComplMatching;
        this->destPg = revComplMatching?reverseComplement(destPg):destPg;

        exactMatchPg();

        if (srcPgPrefix == destPgFilePrefix)
            resolveMappingCollisionsInTheSameText();

        ofstream pgDest = SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(destPgFilePrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX, true);
        ofstream pgMapOffDest = SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(destPgFilePrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX, true);
        ofstream pgMapLenDest = SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(destPgFilePrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX, true);

        sort(textMatches.begin(), textMatches.end(), [this](const TextMatch& match1, const TextMatch& match2) -> bool
        { return match1.posDestText < match2.posDestText; });
        uint_pg_len_max pos = 0;
        uint_pg_len_max totalDestOverlap = 0;
        uint_pg_len_max totalMatched = 0;
        bool isPgLengthStd = srcPg.length() <= UINT32_MAX;
        for(TextMatch& match: textMatches) {
            if (match.posDestText < pos) {
                uint_pg_len_max overflow = pos - match.posDestText;
                if (overflow > match.length) {
                    totalDestOverlap += match.length;
                    match.length = 0;
                    continue;
                }
                totalDestOverlap += overflow;
                match.length -= overflow;
                match.posDestText += overflow;
                if (!revComplMatching)
                    match.posSrcText += overflow;
            }
            if (match.length < minMatchLength) {
                totalDestOverlap += match.length;
                continue;
            }
            totalMatched += match.length;
            PgSAHelpers::writeArray(pgDest, (void*) (destPg.data() + pos), match.posDestText - pos);
            pgDest.put(128);
            if (isPgLengthStd)
                PgSAHelpers::writeValue<uint32_t>(pgMapOffDest, match.posSrcText);
            else
                PgSAHelpers::writeValue<uint64_t>(pgMapOffDest, match.posSrcText);
            PgSAHelpers::writeUIntByteFrugal(pgMapLenDest, match.length - minMatchLength);
            pos = match.endPosDestText();
        }
        PgSAHelpers::writeArray(pgDest, (void*) (destPg.data() + pos), destPg.length() - pos);
        pgDest.close();
        pgMapOffDest.close();
        pgMapLenDest.close();
        SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(destPgFilePrefix, false);

        cout << "Final size of Pg: " << (destPg.length() - totalMatched) << " (removed: " <<
             getTotalMatchStat(totalMatched) << "; " << totalDestOverlap << " chars in overlapped dest symbol)" << endl;
    }

    void SimplePgMatcher::resolveMappingCollisionsInTheSameText() {
        for (TextMatch& match: textMatches) {
            if (match.posSrcText > match.posDestText) {
                uint64_t tmp = match.posSrcText;
                match.posSrcText = match.posDestText;
                match.posDestText = tmp;
            }
            if (revComplMatching && match.endPosSrcText() > match.posDestText) {
                uint64_t margin = (match.endPosSrcText() - match.posDestText + 1) / 2;
                match.length -= margin;
                match.posDestText += margin;
            }
        }
    }

    void SimplePgMatcher::matchPgInPgFiles(string& hqPgSequence, string& lqPgSequence,
            const string &hqPgPrefix, const string &lqPgPrefix, uint_pg_len_max targetMatchLength,
                         bool revComplMatching) {
        clock_t ref_start = clock();
        PgTools::SimplePgMatcher matcher(hqPgPrefix, hqPgSequence, targetMatchLength);
        cout << "Feeding reference pseudogenome finished in " << clock_millis(ref_start) << " msec. " << endl;
        clock_t hq_start = clock();
        matcher.markAndRemoveExactMatches(hqPgPrefix, hqPgSequence, revComplMatching);
        cout << "PgMatching hqPg finished in " << clock_millis(hq_start) << " msec. " << endl;
        clock_t lq_start = clock();
        matcher.markAndRemoveExactMatches(lqPgPrefix, lqPgSequence, revComplMatching);
        cout << "PgMatching lqPg finished in " << clock_millis(lq_start) << " msec. " << endl;
    }
}


