#include "SimplePgMatcher.h"

#include "copmem/CopMEMMatcher.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"

namespace PgTools {

    SimplePgMatcher::SimplePgMatcher(const string& srcPgPrefix, const string& srcPg, uint32_t targetMatchLength,
            uint32_t minMatchLength)
            :srcPgPrefix(srcPgPrefix), srcPg(srcPg), targetMatchLength(targetMatchLength) {
        cout << "Source pseudogenome length: " << srcPg.length() << endl;

        //matcher = new DefaultTextMatcher(srcPg, targetMatchLength);
        matcher = new CopMEMMatcher(srcPg, targetMatchLength, minMatchLength);

    }

    SimplePgMatcher::~SimplePgMatcher() {
        delete(matcher);
    }

    void SimplePgMatcher::exactMatchPg(string& destPg, uint32_t minMatchLength) {
        clock_checkpoint();
        bool destPgIsSrcPg = srcPgPrefix == targetPgPrefix;

        if (!destPgIsSrcPg)
            cout << "Destination pseudogenome length: " << destPgLength << endl;

        if (revComplMatching) {
            if (destPgIsSrcPg) {
                string queryPg = reverseComplement(destPg);
                matcher->matchTexts(textMatches, queryPg, destPgIsSrcPg, revComplMatching, minMatchLength);
            } else {
                reverseComplementInPlace(destPg);
                matcher->matchTexts(textMatches, destPg, destPgIsSrcPg, revComplMatching, minMatchLength);
                reverseComplementInPlace(destPg);
            }
        } else
            matcher->matchTexts(textMatches, destPg, destPgIsSrcPg, revComplMatching, minMatchLength);

        cout << "... found " << textMatches.size() << " exact matches in " << clock_millis() << " msec. " << endl;

        /*        std::sort(textMatches.begin(), textMatches.end(), [](const TextMatch &match1, const TextMatch &match2) -> bool
            { return match1.length > match2.length; });
        cout << "Largest matches:" << endl;
        for (uint32_t i = 0; i < textMatches.size() && i < 10; i++)
            textMatches[i].report(cout);/**/

        if (revComplMatching)
            correctDestPositionDueToRevComplMatching();
    }

    using namespace PgTools;

    void SimplePgMatcher::correctDestPositionDueToRevComplMatching() {
        for (TextMatch& match: textMatches)
            match.posDestText = destPgLength - (match.posDestText + match.length);
    }

    string SimplePgMatcher::getTotalMatchStat(uint_pg_len_max totalMatchLength) {
        return toString(totalMatchLength) + " (" + toString((totalMatchLength * 100.0) / destPgLength, 1)+ "%)";
    }

    static const char MATCH_MARK = 128;

    void SimplePgMatcher::markAndRemoveExactMatches(const string &destPgFilePrefix, string &destPg, bool revComplMatching,
                                                    uint32_t minMatchLength) {
        this->targetPgPrefix = destPgFilePrefix;
        this->revComplMatching = revComplMatching;
        this->destPgLength = destPg.length();

        if (minMatchLength == UINT32_MAX)
            minMatchLength = targetMatchLength;
        exactMatchPg(destPg, minMatchLength);

        clock_t post_start = clock();
        if (srcPgPrefix == destPgFilePrefix)
            resolveMappingCollisionsInTheSameText();

        ofstream pgDest = SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(destPgFilePrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX, true);
        ofstream pgMapOffDest = SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(destPgFilePrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX, true);
        ofstream pgMapLenDest = SeparatedPseudoGenomePersistence::getPseudoGenomeElementDest(destPgFilePrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX, true);

        PgSAHelpers::writeUIntByteFrugal(pgMapLenDest, minMatchLength);

        sort( textMatches.begin(), textMatches.end() );
        textMatches.erase( unique( textMatches.begin(), textMatches.end() ), textMatches.end() );
        cout << "Unique exact matches: " << textMatches.size() << endl;

        uint_pg_len_max pos = 0;
        uint_pg_len_max nPos = 0;
        uint_pg_len_max totalDestOverlap = 0;
        uint_pg_len_max totalMatched = 0;
        bool isPgLengthStd = srcPg.length() <= UINT32_MAX;
        for(TextMatch& match: textMatches) {
            if (match.posDestText < pos) {
                uint_pg_len_max overflow = pos - match.posDestText;
                if (overflow >= match.length) {
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
            uint64_t length = match.posDestText - pos;
            destPg.replace(nPos, length, destPg, pos, length);
            nPos += length;
            destPg[nPos++] = MATCH_MARK;
            if (isPgLengthStd)
                PgSAHelpers::writeValue<uint32_t>(pgMapOffDest, match.posSrcText);
            else
                PgSAHelpers::writeValue<uint64_t>(pgMapOffDest, match.posSrcText);
            PgSAHelpers::writeUIntByteFrugal(pgMapLenDest, match.length - minMatchLength);
            pos = match.endPosDestText();
        }
        uint64_t length = destPg.length() - pos;
        destPg.replace(nPos, length, destPg, pos, length);
        nPos += length;
        destPg.resize(nPos);
        PgSAHelpers::writeArray(pgDest, (void*) (destPg.data()), destPg.length());
        pgDest.close();
        pgMapOffDest.close();
        pgMapLenDest.close();
        SeparatedPseudoGenomePersistence::acceptTemporaryPseudoGenomeElements(destPgFilePrefix, false);
        cout << "Writing files time: " << clock_millis(post_start) << endl;
        cout << "Final size of Pg: " << nPos << " (removed: " <<
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
            const string &hqPgPrefix, const string &lqPgPrefix, uint_pg_len_max targetMatchLength
            , uint32_t minMatchLength) {
        clock_t ref_start = clock();
        PgTools::SimplePgMatcher matcher(hqPgPrefix, hqPgSequence, targetMatchLength, minMatchLength);
        cout << "Feeding reference pseudogenome finished in " << clock_millis(ref_start) << " msec. " << endl;
        clock_t lq_start = clock();
        matcher.markAndRemoveExactMatches(lqPgPrefix, lqPgSequence, true, minMatchLength);
        cout << "PgMatching lqPg finished in " << clock_millis(lq_start) << " msec. " << endl;
        clock_t hq_start = clock();
        matcher.markAndRemoveExactMatches(hqPgPrefix, hqPgSequence, true, minMatchLength);
        cout << "PgMatching hqPg finished in " << clock_millis(hq_start) << " msec. " << endl;

    }

    string
    SimplePgMatcher::restoreMatchedPg(string &srcPg, const string &destPgPrefix, bool revComplMatching,
                                      bool plainTextReadMode, bool srcIsDest) {
        string destPg = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(destPgPrefix);
        bool isPgLengthStd = srcPg.length() <= UINT32_MAX;
        if (srcIsDest)
            srcPg.resize(0);
        string tmp;
        string& resPg = srcIsDest?srcPg:tmp;
        uint64_t posDest = 0;
        uint32_t minMatchLength = 0;
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        PgSAHelpers::readUIntByteFrugal(pgMapLenSrc, minMatchLength);
        uint64_t markPos = 0;
        while ((markPos = destPg.find(MATCH_MARK, posDest)) != std::string::npos) {
            resPg.append(destPg, posDest, markPos - posDest);
            posDest = markPos + 1;
            uint64_t matchSrcPos = 0;
            if (isPgLengthStd) {
                uint32_t tmp;
                PgSAHelpers::readValue<uint32_t>(pgMapOffSrc, tmp, plainTextReadMode);
                matchSrcPos = tmp;
            } else
                PgSAHelpers::readValue<uint64_t>(pgMapOffSrc, matchSrcPos, plainTextReadMode);
            uint16_t matchLength = 0;
            PgSAHelpers::readUIntByteFrugal(pgMapLenSrc, matchLength);
            matchLength += minMatchLength;
            if (revComplMatching)
                resPg.append(reverseComplement(srcPg.substr(matchSrcPos, matchLength)));
            else
                resPg.append(srcPg.substr(matchSrcPos, matchLength));
        }
        resPg.append(destPg, posDest, destPg.length() - posDest);

        cout << "Restored Pg sequence of length: " << resPg.length() << endl;

        return resPg;
    }

    string
    SimplePgMatcher::restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching) {
        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string tmpPg;
        tmpPg.resize(pgh->getPseudoGenomeLength());
        tmpPg = SimplePgMatcher::restoreMatchedPg(tmpPg, pgPrefix, true, plainTextReadMode, true);
        delete(pgh);
        delete(prop);
        return tmpPg;
    }
}


