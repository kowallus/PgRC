#include "SimplePgMatcher.h"

#include "copmem/CopMEMMatcher.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../utils/LzmaLib.h"

namespace PgTools {

    SimplePgMatcher::SimplePgMatcher(const string& srcPg, uint32_t targetMatchLength,
            uint32_t minMatchLength)
            :srcPg(srcPg), targetMatchLength(targetMatchLength) {
        cout << "Source pseudogenome length: " << srcPg.length() << endl;

        //matcher = new DefaultTextMatcher(srcPg, targetMatchLength);
        matcher = new CopMEMMatcher(srcPg, targetMatchLength, minMatchLength);

    }

    SimplePgMatcher::~SimplePgMatcher() {
        delete(matcher);
    }

    void SimplePgMatcher::exactMatchPg(string& destPg, bool destPgIsSrcPg, uint32_t minMatchLength) {
        clock_checkpoint();

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

    void SimplePgMatcher::markAndRemoveExactMatches(
            bool destPgIsSrcPg, string &destPg, bool revComplMatching, uint32_t minMatchLength) {
        this->revComplMatching = revComplMatching;
        this->destPgLength = destPg.length();

        if (minMatchLength == UINT32_MAX)
            minMatchLength = targetMatchLength;
        exactMatchPg(destPg, destPgIsSrcPg, minMatchLength);

        clock_t post_start = clock();
        if (destPgIsSrcPg)
            resolveMappingCollisionsInTheSameText();

        stringstream pgDest;
        stringstream pgMapOffDest;
        stringstream pgMapLenDest;

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

        pgMapped = std::move(pgDest.str());
        pgMapOff = std::move(pgMapOffDest.str());
        pgMapLen = std::move(pgMapLenDest.str());

        cout << "Preparing output time: " << clock_millis(post_start) << endl;
        cout << "Final size of Pg: " << nPos << " (removed: " <<
             getTotalMatchStat(totalMatched) << "; " << totalDestOverlap << " chars in overlapped dest symbol)" << endl;
    }

    void SimplePgMatcher::writeMatchingResult(const string &pgPrefix) {
        PgSAHelpers::writeArrayToFile(pgPrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_FILE_SUFFIX,
                                      (void*) pgMapped.data(), pgMapped.length());
        PgSAHelpers::writeArrayToFile(pgPrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX,
                                      (void*) pgMapOff.data(), pgMapOff.length());
        PgSAHelpers::writeArrayToFile(pgPrefix + SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX,
                                      (void*) pgMapLen.data(), pgMapLen.length());
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

    void SimplePgMatcher::matchPgsInPg(string &hqPgSequence, string &lqPgSequence, ostream &pgrcOut,
                                       const string &hqPgPrefix, const string &lqPgPrefix,
                                       uint_pg_len_max targetMatchLength, uint32_t minMatchLength) {
        clock_t ref_start = clock();
        PgTools::SimplePgMatcher matcher(hqPgSequence, targetMatchLength, minMatchLength);
        cout << "Feeding reference pseudogenome finished in " << clock_millis(ref_start) << " msec. " << endl;
        clock_t lq_start = clock();
        matcher.markAndRemoveExactMatches(false, lqPgSequence, true, minMatchLength);
        if (!lqPgPrefix.empty())
            matcher.writeMatchingResult(lqPgPrefix);
        string lqPgMapped = std::move(matcher.pgMapped);
        string lqPgMapOff = std::move(matcher.pgMapOff);
        string lqPgMapLen = std::move(matcher.pgMapLen);
        cout << "PgMatching lqPg finished in " << clock_millis(lq_start) << " msec. " << endl;
        clock_t hq_start = clock();
        matcher.markAndRemoveExactMatches(true, hqPgSequence, true, minMatchLength);
        if (!hqPgPrefix.empty())
            matcher.writeMatchingResult(hqPgPrefix);
        string pgSeq = std::move(matcher.pgMapped);
        string hqPgMapOff = std::move(matcher.pgMapOff);
        string hqPgMapLen = std::move(matcher.pgMapLen);
        cout << "PgMatching hqPg finished in " << clock_millis(hq_start) << " msec. " << endl;

        PgSAHelpers::writeValue<uint_pg_len_max>(pgrcOut, pgSeq.length(), false);

        pgSeq.append(std::move(lqPgMapped));
        lqPgMapped.clear();
        writeCompressed(pgrcOut, pgSeq.data(), pgSeq.size(), LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM);
        writeCompressed(pgrcOut, hqPgMapOff.data(), hqPgMapOff.size(), LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM);
        writeCompressed(pgrcOut, hqPgMapLen.data(), hqPgMapLen.size(), LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM);
        writeCompressed(pgrcOut, lqPgMapOff.data(), lqPgMapOff.size(), LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM);
        writeCompressed(pgrcOut, lqPgMapLen.data(), lqPgMapLen.size(), LZMA_CODER, PGRC_CODER_LEVEL_MAXIMUM);
    }

    void SimplePgMatcher::restoreMatchedPgs(istream &pgrcIn, string &hqPgSequence, string &lqPgSequence) {
        uint_pg_len_max hqPgSeqLen = 0;
        PgSAHelpers::readValue<uint_pg_len_max>(pgrcIn, hqPgSeqLen, false);
        string comboPgSeq;
        readCompressed(pgrcIn, comboPgSeq);
        size_t lqPgSeqLen = comboPgSeq.length() - hqPgSeqLen;
        string lqPgMapped(comboPgSeq, hqPgSeqLen);
        comboPgSeq.resize(hqPgSeqLen);
        string hqPgMapped = std::move(comboPgSeq);
        string pgMapOff, pgMapLen;
        readCompressed(pgrcIn, pgMapOff);
        readCompressed(pgrcIn, pgMapLen);
        istringstream pgMapOffSrc, pgMapLenSrc;
        pgMapOffSrc.str(pgMapOff);
        pgMapLenSrc.str(pgMapLen);
        hqPgSequence.clear();
        hqPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, hqPgSeqLen, hqPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                  true, false, true);
        readCompressed(pgrcIn, pgMapOff);
        readCompressed(pgrcIn, pgMapLen);
        pgMapOffSrc.str(pgMapOff);
        pgMapLenSrc.str(pgMapLen);
        lqPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, lqPgSeqLen, lqPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                         true, false);
    }

    string
    SimplePgMatcher::restoreMatchedPg(string &srcPg, size_t srcLen, const string& destPg, istream &pgMapOffSrc, istream &pgMapLenSrc,
            bool revComplMatching, bool plainTextReadMode, bool srcIsDest) {
        bool isPgLengthStd = srcLen <= UINT32_MAX;
        if (srcIsDest)
            srcPg.resize(0);
        string tmp;
        string& resPg = srcIsDest?srcPg:tmp;
        uint64_t posDest = 0;
        uint32_t minMatchLength = 0;

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

    string SimplePgMatcher::restoreMatchedPg(string &srcPg, const string& destPgPrefix,
                                   bool revComplMatching, bool plainTextReadMode) {
        string destPg = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(destPgPrefix);
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        return SimplePgMatcher::restoreMatchedPg(srcPg, srcPg.length(), destPg, pgMapOffSrc, pgMapLenSrc, revComplMatching, plainTextReadMode, false);
    }

    string
    SimplePgMatcher::restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching) {
        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string destPg = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(pgPrefix);
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(pgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(pgPrefix, SeparatedPseudoGenomePersistence::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        string resPg;
        resPg = SimplePgMatcher::restoreMatchedPg(resPg, pgh->getPseudoGenomeLength(), destPg, pgMapOffSrc, pgMapLenSrc,
                revComplMatching, plainTextReadMode, true);
        delete(pgh);
        delete(prop);
        return resPg;
    }
}


