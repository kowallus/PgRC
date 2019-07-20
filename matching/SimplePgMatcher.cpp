#include "SimplePgMatcher.h"

#include "copmem/CopMEMMatcher.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../utils/LzmaLib.h"

namespace PgTools {

    SimplePgMatcher::SimplePgMatcher(const string& srcPg, uint32_t targetMatchLength,
            uint32_t minMatchLength)
            :srcPg(srcPg), targetMatchLength(targetMatchLength) {
        cout << "Source pseudogenome length: " << srcPg.length() << endl;
        if (srcPg.size() >= targetMatchLength)
            matcher = new CopMEMMatcher(srcPg.data(), srcPg.length(), targetMatchLength, minMatchLength);
        //matcher = new DefaultTextMatcher(srcPg, targetMatchLength);
    }

    SimplePgMatcher::~SimplePgMatcher() {
        if (matcher)
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

        *logout << "... found " << textMatches.size() << " exact matches in " << clock_millis() << " msec. " << endl;

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
            bool destPgIsSrcPg, string &destPg, string &resPgMapOff, string& resPgMapLen,
            bool revComplMatching, uint32_t minMatchLength) {
        if (!matcher) {
            resPgMapOff.clear();
            resPgMapLen.clear();
            return;
        }

        this->revComplMatching = revComplMatching;
        this->destPgLength = destPg.length();

        if (minMatchLength == UINT32_MAX)
            minMatchLength = targetMatchLength;
        exactMatchPg(destPg, destPgIsSrcPg, minMatchLength);

        clock_t post_start = clock();
        if (destPgIsSrcPg)
            resolveMappingCollisionsInTheSameText();

        ostringstream pgMapOffDest;
        ostringstream pgMapLenDest;

        PgSAHelpers::writeUIntByteFrugal(pgMapLenDest, minMatchLength);

        sort( textMatches.begin(), textMatches.end() );
        textMatches.erase( unique( textMatches.begin(), textMatches.end() ), textMatches.end() );
        *logout << "Unique exact matches: " << textMatches.size() << endl;

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

        textMatches.clear();
        resPgMapOff = pgMapOffDest.str();
        pgMapOffDest.clear();
        resPgMapLen = pgMapLenDest.str();
        pgMapLenDest.clear();

        *logout << "Preparing output time: " << clock_millis(post_start) << endl;
        cout << "Final size of Pg: " << nPos << " (removed: " <<
             getTotalMatchStat(totalMatched) << "; " << totalDestOverlap << " chars in overlapped dest symbol)" << endl;
    }

    void SimplePgMatcher::writeMatchingResult(const string &pgPrefix,
                                              const string &pgMapped, const string &pgMapOff, const string& pgMapLen) {
        PgSAHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX,
                                      pgMapped);
        PgSAHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX,
                                      pgMapOff);
        PgSAHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX,
                                      pgMapLen);
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

    void SimplePgMatcher::matchPgsInPg(string &hqPgSequence, string &lqPgSequence, string &nPgSequence,
                                        bool separateNReads, ostream &pgrcOut, uint8_t coder_level,
                                        const string &hqPgPrefix, const string &lqPgPrefix, const string &nPgPrefix,
                                        uint_pg_len_max targetMatchLength, uint32_t minMatchLength) {
        clock_t ref_start = clock();
        const unsigned long refSequenceLength = hqPgSequence.length();
        bool isPgLengthStd = refSequenceLength <= UINT32_MAX;
        
        PgTools::SimplePgMatcher* matcher = new PgTools::SimplePgMatcher(hqPgSequence, targetMatchLength, minMatchLength);
        *logout << "Feeding reference pseudogenome finished in " << clock_millis(ref_start) << " msec. " << endl;
        clock_t lq_start = clock();
        string lqPgMapOff, lqPgMapLen;
        matcher->markAndRemoveExactMatches(false, lqPgSequence, lqPgMapOff, lqPgMapLen, true, minMatchLength);
        if (!lqPgPrefix.empty())
            matcher->writeMatchingResult(lqPgPrefix, lqPgSequence, lqPgMapOff, lqPgMapLen);
        *logout << "PgMatching lqPg finished in " << clock_millis(lq_start) << " msec. " << endl;

        clock_t n_start = clock();
        string nPgMapOff, nPgMapLen;
        matcher->markAndRemoveExactMatches(false, nPgSequence, nPgMapOff, nPgMapLen, true, minMatchLength);
        if (!nPgPrefix.empty())
            matcher->writeMatchingResult(nPgPrefix, nPgSequence, nPgMapOff, nPgMapLen);
        if (!nPgSequence.empty())
            *logout << "PgMatching nPg finished in " << clock_millis(n_start) << " msec. " << endl;

        clock_t hq_start = clock();
        string hqPgMapOff, hqPgMapLen;
        matcher->markAndRemoveExactMatches(true, hqPgSequence, hqPgMapOff, hqPgMapLen, true, minMatchLength);
        if (!hqPgPrefix.empty())
            matcher->writeMatchingResult(hqPgPrefix, hqPgSequence, hqPgMapOff, hqPgMapLen);
        *logout << "PgMatching hqPg finished in " << clock_millis(hq_start) << " msec. " << endl;
        delete(matcher);

        PgSAHelpers::writeValue<uint_pg_len_max>(pgrcOut, hqPgSequence.length(), false);
        PgSAHelpers::writeValue<uint_pg_len_max>(pgrcOut, lqPgSequence.length(), false);
        PgSAHelpers::writeValue<uint_pg_len_max>(pgrcOut, nPgSequence.length(), false);
        string comboPgSeq = std::move(hqPgSequence);
        comboPgSeq.reserve(comboPgSeq.size() + lqPgSequence.size() + nPgSequence.size());
        comboPgSeq.append(lqPgSequence);
        lqPgSequence.clear();
        lqPgSequence.shrink_to_fit();
        comboPgSeq.append(nPgSequence);
        nPgSequence.clear();
        nPgSequence.shrink_to_fit();

        *logout << "Joined mapped sequences (good&bad" << (nPgSequence.empty()?"":"&N") << ")... ";
        writeCompressed(pgrcOut, comboPgSeq.data(), comboPgSeq.size(), LZMA_CODER, coder_level,
                PGRC_DATAPERIODCODE_8_t, COMPRESSION_ESTIMATION_BASIC_DNA);
        *logout << "Good sequence mapping - offsets... ";
        double estimated_pg_offset_ratio = simpleUintCompressionEstimate(refSequenceLength, isPgLengthStd?UINT32_MAX:UINT64_MAX);
        const int pgrc_pg_offset_dataperiodcode = isPgLengthStd ? PGRC_DATAPERIODCODE_32_t : PGRC_DATAPERIODCODE_64_t;
        writeCompressed(pgrcOut, hqPgMapOff.data(), hqPgMapOff.size(), LZMA_CODER, coder_level,
                        pgrc_pg_offset_dataperiodcode, estimated_pg_offset_ratio);
        *logout << "lengths... ";
        writeCompressed(pgrcOut, hqPgMapLen.data(), hqPgMapLen.size(), LZMA_CODER, coder_level,
                        PGRC_DATAPERIODCODE_8_t);
        *logout << "Bad sequence mapping - offsets... ";
        writeCompressed(pgrcOut, lqPgMapOff.data(), lqPgMapOff.size(), LZMA_CODER, coder_level,
                        pgrc_pg_offset_dataperiodcode, estimated_pg_offset_ratio);
        *logout << "lengths... ";
        writeCompressed(pgrcOut, lqPgMapLen.data(), lqPgMapLen.size(), LZMA_CODER, coder_level,
                        PGRC_DATAPERIODCODE_8_t);
        if (separateNReads) {
            *logout << "N sequence mapping - offsets... ";
            writeCompressed(pgrcOut, nPgMapOff.data(), nPgMapOff.size(), LZMA_CODER, coder_level,
                            pgrc_pg_offset_dataperiodcode, estimated_pg_offset_ratio);
            *logout << "lengths... ";
            writeCompressed(pgrcOut, nPgMapLen.data(), nPgMapLen.size(), LZMA_CODER, coder_level,
                            PGRC_DATAPERIODCODE_8_t);
        }
    }

    void SimplePgMatcher::restoreMatchedPgs(istream &pgrcIn, uint_pg_len_max orgHqPgLen, string &hqPgSequence, string &lqPgSequence,
            string &nPgSequence) {
        uint_pg_len_max hqPgMappedLen, lqPgMappedLen, nPgMappedLen;
        string pgMapOff, pgMapLen;
        istringstream pgMapOffSrc, pgMapLenSrc;
        PgSAHelpers::readValue<uint_pg_len_max>(pgrcIn, hqPgMappedLen, false);
        PgSAHelpers::readValue<uint_pg_len_max>(pgrcIn, lqPgMappedLen, false);
        PgSAHelpers::readValue<uint_pg_len_max>(pgrcIn, nPgMappedLen, false);
        string comboPgMapped;
        readCompressed(pgrcIn, comboPgMapped);
        string nPgMapped(comboPgMapped, hqPgMappedLen + lqPgMappedLen);
        comboPgMapped.resize(hqPgMappedLen + lqPgMappedLen);
        {
            string lqPgMapped(comboPgMapped, hqPgMappedLen);
            comboPgMapped.resize(hqPgMappedLen);
            comboPgMapped.shrink_to_fit();
            {
                string hqPgMapped = std::move(comboPgMapped);
                readCompressed(pgrcIn, pgMapOff);
                readCompressed(pgrcIn, pgMapLen);
                pgMapOffSrc.str(pgMapOff);
                pgMapLenSrc.str(pgMapLen);
                hqPgSequence.clear();
                hqPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, orgHqPgLen, hqPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                                 true, false, true);
            }
            readCompressed(pgrcIn, pgMapOff);
            readCompressed(pgrcIn, pgMapLen);
            pgMapOffSrc.str(pgMapOff);
            pgMapLenSrc.str(pgMapLen);
            lqPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, orgHqPgLen, lqPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                             true, false);
        }
        if (nPgMappedLen) {
            readCompressed(pgrcIn, pgMapOff);
            readCompressed(pgrcIn, pgMapLen);
            pgMapOffSrc.str(pgMapOff);
            pgMapLenSrc.str(pgMapLen);
            nPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, orgHqPgLen, nPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                            true, false);
        }
    }

    string
    SimplePgMatcher::restoreMatchedPg(string &srcPg, size_t orgSrcLen, const string& destPg, istream &pgMapOffSrc, istream &pgMapLenSrc,
            bool revComplMatching, bool plainTextReadMode, bool srcIsDest) {
        bool isPgLengthStd = orgSrcLen <= UINT32_MAX;
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
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        return SimplePgMatcher::restoreMatchedPg(srcPg, srcPg.length(), destPg, pgMapOffSrc, pgMapLenSrc, revComplMatching, plainTextReadMode, false);
    }

    string
    SimplePgMatcher::restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching) {
        PseudoGenomeHeader* pgh = 0;
        ReadsSetProperties* prop = 0;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string destPg = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(pgPrefix);
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        string resPg;
        resPg = SimplePgMatcher::restoreMatchedPg(resPg, pgh->getPseudoGenomeLength(), destPg, pgMapOffSrc, pgMapLenSrc,
                revComplMatching, plainTextReadMode, true);
        delete(pgh);
        delete(prop);
        return resPg;
    }
}


