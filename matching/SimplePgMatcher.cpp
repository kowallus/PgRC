#include "SimplePgMatcher.h"

#include "copmem/CopMEMMatcher.h"

#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "../coders/CodersLib.h"
#include "../coders/VarLenDNACoder.h"
#include "../coders/PropsLibrary.h"

namespace PgTools {

    SimplePgMatcher::SimplePgMatcher(const string &srcPg, uint32_t targetMatchLength,
                                     uint32_t minMatchLength)
            : srcPg(srcPg), targetMatchLength(targetMatchLength) {
        cout << "Source pseudogenome length: " << srcPg.length() << endl;
        if (srcPg.size() >= targetMatchLength)
            matcher = new CopMEMMatcher(srcPg.data(), srcPg.length(), targetMatchLength, minMatchLength);
        //matcher = new DefaultTextMatcher(srcPg, targetMatchLength);
    }

    SimplePgMatcher::~SimplePgMatcher() {
        if (matcher)
            delete (matcher);
    }

    void SimplePgMatcher::exactMatchPg(string &destPg, bool destPgIsSrcPg, uint32_t minMatchLength) {
        time_checkpoint();

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

        *logout << "... found " << textMatches.size() << " exact matches in " << time_millis() << " msec. " << endl;

        /*        std::sort(textMatches.begin(), textMatches.end(), [](const TextMatch &match1, const TextMatch &match2) -> bool
            { return match1.length > match2.length; });
        cout << "Largest matches:" << endl;
        for (uint32_t i = 0; i < textMatches.size() && i < 10; i++)
            textMatches[i].report(cout);*/

        if (revComplMatching)
            correctDestPositionDueToRevComplMatching();
    }

    using namespace PgTools;

    void SimplePgMatcher::correctDestPositionDueToRevComplMatching() {
        for (TextMatch &match: textMatches)
            match.posDestText = destPgLength - (match.posDestText + match.length);
    }

    string SimplePgMatcher::getTotalMatchStat(uint_pg_len_max totalMatchLength) {
        return toString(totalMatchLength) + " (" + toString((totalMatchLength * 100.0) / destPgLength, 1) + "%)";
    }

    char SimplePgMatcher::MATCH_MARK = '%';

    void SimplePgMatcher::markAndRemoveExactMatches(
            bool destPgIsSrcPg, string &destPg, string &resPgMapOff, string &resPgMapLen,
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

        chrono::steady_clock::time_point post_start = chrono::steady_clock::now();
        if (destPgIsSrcPg)
            resolveMappingCollisionsInTheSameText();

        ostringstream pgMapOffDest;
        ostringstream pgMapLenDest;

        PgHelpers::writeUIntByteFrugal(pgMapLenDest, minMatchLength);

        sort(textMatches.begin(), textMatches.end());
        textMatches.erase(unique(textMatches.begin(), textMatches.end()), textMatches.end());
        *logout << "Unique exact matches: " << textMatches.size() << endl;

        char *destPtr = (char *) destPg.data();
        uint_pg_len_max pos = 0;
        uint_pg_len_max nPos = 0;
        uint_pg_len_max totalDestOverlap = 0;
        uint_pg_len_max totalMatched = 0;
        bool isPgLengthStd = srcPg.length() <= UINT32_MAX;
        for (TextMatch &match: textMatches) {
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
            memmove(destPtr + nPos, destPtr + pos, length);
            nPos += length;
            destPg[nPos++] = MATCH_MARK;
            if (isPgLengthStd)
                PgHelpers::writeValue<uint32_t>(pgMapOffDest, match.posSrcText);
            else
                PgHelpers::writeValue<uint64_t>(pgMapOffDest, match.posSrcText);
            PgHelpers::writeUIntByteFrugal(pgMapLenDest, match.length - minMatchLength);
            pos = match.endPosDestText();
        }
        uint64_t length = destPg.length() - pos;
        memmove(destPtr + nPos, destPtr + pos, length);
        nPos += length;
        destPg.resize(nPos);

        textMatches.clear();
        resPgMapOff = pgMapOffDest.str();
        pgMapOffDest.clear();
        resPgMapLen = pgMapLenDest.str();
        pgMapLenDest.clear();

        *logout << "Preparing output time: " << time_millis(post_start) << " msec." << endl;
        cout << "Final size of Pg: " << nPos << " (removed: " <<
             getTotalMatchStat(totalMatched) << "; " << totalDestOverlap << " chars in overlapped dest symbol)" << endl;
    }

    void SimplePgMatcher::writeMatchingResult(const string &pgPrefix,
                                              const string &pgMapped, const string &pgMapOff, const string &pgMapLen) {
        PgHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_FILE_SUFFIX,
                                       pgMapped);
        PgHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX,
                                       pgMapOff);
        PgHelpers::writeStringToFile(pgPrefix + SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX,
                                       pgMapLen);
    }

    void SimplePgMatcher::resolveMappingCollisionsInTheSameText() {
        for (TextMatch &match: textMatches) {
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
        chrono::steady_clock::time_point ref_start = chrono::steady_clock::now();
        const unsigned long refSequenceLength = hqPgSequence.length();
        bool isPgLengthStd = refSequenceLength <= UINT32_MAX;
        
        PgTools::SimplePgMatcher* matcher = new PgTools::SimplePgMatcher(hqPgSequence, targetMatchLength, minMatchLength);
        *logout << "Feeding reference pseudogenome finished in " << time_millis(ref_start) << " msec. " << endl;
        chrono::steady_clock::time_point lq_start = chrono::steady_clock::now();
        string lqPgMapOff, lqPgMapLen;
        matcher->markAndRemoveExactMatches(false, lqPgSequence, lqPgMapOff, lqPgMapLen, true, minMatchLength);
        if (!lqPgPrefix.empty())
            matcher->writeMatchingResult(lqPgPrefix, lqPgSequence, lqPgMapOff, lqPgMapLen);
        *logout << "PgMatching lqPg finished in " << time_millis(lq_start) << " msec. " << endl;

        chrono::steady_clock::time_point n_start = chrono::steady_clock::now();
        string nPgMapOff, nPgMapLen;
        matcher->markAndRemoveExactMatches(false, nPgSequence, nPgMapOff, nPgMapLen, true, minMatchLength);
        if (!nPgPrefix.empty())
            matcher->writeMatchingResult(nPgPrefix, nPgSequence, nPgMapOff, nPgMapLen);
        if (!nPgSequence.empty())
            *logout << "PgMatching nPg finished in " << time_millis(n_start) << " msec. " << endl;

        chrono::steady_clock::time_point hq_start = chrono::steady_clock::now();
        string hqPgMapOff, hqPgMapLen;
        matcher->markAndRemoveExactMatches(true, hqPgSequence, hqPgMapOff, hqPgMapLen, true, minMatchLength);
        if (!hqPgPrefix.empty())
            matcher->writeMatchingResult(hqPgPrefix, hqPgSequence, hqPgMapOff, hqPgMapLen);
        *logout << "PgMatching hqPg finished in " << time_millis(hq_start) << " msec. " << endl;
        delete(matcher);

        ostringstream pgsLen;
        PgHelpers::writeValue<uint_pg_len_max>(pgsLen, hqPgSequence.length(), false);
        PgHelpers::writeValue<uint_pg_len_max>(pgsLen, lqPgSequence.length(), false);
        PgHelpers::writeValue<uint_pg_len_max>(pgsLen, nPgSequence.length(), false);
        string comboPgSeq = std::move(hqPgSequence);
        comboPgSeq.reserve(comboPgSeq.size() + lqPgSequence.size() + nPgSequence.size());
        comboPgSeq.append(lqPgSequence);
        lqPgSequence.clear();
        lqPgSequence.shrink_to_fit();
        comboPgSeq.append(nPgSequence);
        nPgSequence.clear();
        nPgSequence.shrink_to_fit();

        string pgsLenString = pgsLen.str();
        vector<CompressionJob> cJobs;
        auto fseCoderProps = getDefaultFSECoderProps();
        cJobs.emplace_back("Sequences length info... ", pgsLenString, fseCoderProps.get());
        bool noNPgSequence = nPgSequence.empty();
        auto dnaCoderProps = getDefaultCoderProps(VARLEN_DNA_CODER, coder_level);
        auto pgSeqCoderProps = getVarLenEncodedPgCoderProps(coder_level);
        auto compoundCoderProps = getCompoundCoderProps(dnaCoderProps.get(), pgSeqCoderProps.get());
        cJobs.emplace_back(string("Joined mapped sequences (good&bad") + (noNPgSequence ? "" : "&N") + ")... ",
                           (unsigned char*) comboPgSeq.data(), comboPgSeq.size(),
                           compoundCoderProps.get(), COMPRESSION_ESTIMATION_VAR_LEN_DNA);
        double estimated_pg_offset_ratio = simpleUintCompressionEstimate(refSequenceLength, isPgLengthStd ? UINT32_MAX : UINT64_MAX);
        const int pgrc_pg_offset_dataperiodcode = isPgLengthStd ? LZMA_DATAPERIODCODE_32_t : LZMA_DATAPERIODCODE_64_t;
        auto hqMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, pgrc_pg_offset_dataperiodcode);
        cJobs.emplace_back("Good sequence mapping - offsets... ", (unsigned char*) hqPgMapOff.data(),
                           hqPgMapOff.size(), hqMapOffCoderProps.get(), estimated_pg_offset_ratio);
        auto hqMapLenCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, LZMA_DATAPERIODCODE_8_t);
        cJobs.emplace_back("lengths... ", (unsigned char*) hqPgMapLen.data(), hqPgMapLen.size(),
                           hqMapLenCoderProps.get());
        auto lqMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, pgrc_pg_offset_dataperiodcode);
        cJobs.emplace_back("Bad sequence mapping - offsets... ", (unsigned char*) lqPgMapOff.data(),
                           lqPgMapOff.size(), lqMapOffCoderProps.get(), estimated_pg_offset_ratio);
        auto lqMapLenCoderProps = getDefaultFSECoderProps(12);
        cJobs.emplace_back("lengths... ", (unsigned char*) lqPgMapLen.data(), lqPgMapLen.size(),
                           lqMapLenCoderProps.get());
        auto nMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, coder_level, pgrc_pg_offset_dataperiodcode);
        auto nMapLenCoderProps =  getDefaultFSECoderProps();
        if (separateNReads) {
            cJobs.emplace_back("N sequence mapping - offsets... ", (unsigned char*) nPgMapOff.data(),
                               nPgMapOff.size(), nMapOffCoderProps.get(),estimated_pg_offset_ratio);


            cJobs.emplace_back("lengths... ", (unsigned char*) nPgMapLen.data(), nPgMapLen.size(),
                               nMapLenCoderProps.get());
        }
        CompressionJob::writeCompressedCollectiveParallel(pgrcOut, cJobs);
    }

    void SimplePgMatcher::restoreMatchedPgs(istream &pgrcIn, uint_pg_len_max orgHqPgLen, string &hqPgSequence, string &lqPgSequence,
                                            string &nPgSequence, PgRCParams* params) {
        uint_pg_len_max hqPgMappedLen, lqPgMappedLen, nPgMappedLen;
        string hqPgMapOff, hqPgMapLen, lqPgMapOff, lqPgMapLen, nPgMapOff, nPgMapLen;
        istringstream pgMapOffSrc, pgMapLenSrc;
        istream* propsIn = &pgrcIn;
        string propsString;
        if (params->isVersionAtLeast(1, 3)) {
            readCompressed(pgrcIn, propsString);
            propsIn = new istringstream(propsString);
        }
        PgHelpers::readValue<uint_pg_len_max>(*propsIn, hqPgMappedLen, false);
        PgHelpers::readValue<uint_pg_len_max>(*propsIn, lqPgMappedLen, false);
        PgHelpers::readValue<uint_pg_len_max>(*propsIn, nPgMappedLen, false);
        if (params->isVersionAtLeast(1, 3))
            delete propsIn;
        string comboPgMapped;
        vector<string*> destStrings;
        destStrings.push_back(&comboPgMapped);
        destStrings.push_back(&hqPgMapOff);
        destStrings.push_back(&hqPgMapLen);
        destStrings.push_back(&lqPgMapOff);
        destStrings.push_back(&lqPgMapLen);
        if (nPgMappedLen) {
            destStrings.push_back(&nPgMapOff);
            destStrings.push_back(&nPgMapLen);
        }
        readCompressedCollectiveParallel(pgrcIn, destStrings, params->isVersion(1, 2) ? 0 : -1);
        string nPgMapped(comboPgMapped, hqPgMappedLen + lqPgMappedLen);
        comboPgMapped.resize(hqPgMappedLen + lqPgMappedLen);
        {
            string lqPgMapped(comboPgMapped, hqPgMappedLen);
            comboPgMapped.resize(hqPgMappedLen);
            comboPgMapped.shrink_to_fit();
            {
                string hqPgMapped = std::move(comboPgMapped);
                pgMapOffSrc.str(hqPgMapOff);
                pgMapLenSrc.str(hqPgMapLen);
                hqPgSequence.clear();
                hqPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, orgHqPgLen, hqPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                                 true, false, true);
            }

            pgMapOffSrc.str(lqPgMapOff);
            pgMapLenSrc.str(lqPgMapLen);
            lqPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, orgHqPgLen, lqPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                             true, false);
        }
        if (nPgMappedLen) {
            pgMapOffSrc.str(nPgMapOff);
            pgMapLenSrc.str(nPgMapLen);
            nPgSequence = SimplePgMatcher::restoreMatchedPg(hqPgSequence, orgHqPgLen, nPgMapped, pgMapOffSrc, pgMapLenSrc,
                                                            true, false);
        }
    }

    string
    SimplePgMatcher::restoreMatchedPg(string &srcPg, size_t orgSrcLen, const string &destPg, istream &pgMapOffSrc, istream &pgMapLenSrc,
                                      bool revComplMatching, bool plainTextReadMode, bool srcIsDest) {
        bool isPgLengthStd = orgSrcLen <= UINT32_MAX;
        if (srcIsDest)
            srcPg.resize(0);
        string tmp;
        string &resPg = srcIsDest ? srcPg : tmp;
        uint64_t posDest = 0;
        uint32_t minMatchLength = 0;

        PgHelpers::readUIntByteFrugal(pgMapLenSrc, minMatchLength);
        uint64_t markPos = 0;
        while ((markPos = destPg.find(MATCH_MARK, posDest)) != std::string::npos) {
            resPg.append(destPg, posDest, markPos - posDest);
            posDest = markPos + 1;
            uint64_t matchSrcPos = 0;
            if (isPgLengthStd) {
                uint32_t tmp;
                PgHelpers::readValue<uint32_t>(pgMapOffSrc, tmp, plainTextReadMode);
                matchSrcPos = tmp;
            } else
                PgHelpers::readValue<uint64_t>(pgMapOffSrc, matchSrcPos, plainTextReadMode);
            uint64_t matchLength = 0;
            PgHelpers::readUIntByteFrugal(pgMapLenSrc, matchLength);
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

    string SimplePgMatcher::restoreMatchedPg(string &srcPg, const string &destPgPrefix,
                                             bool revComplMatching, bool plainTextReadMode) {
        string destPg = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(destPgPrefix);
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(destPgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        return SimplePgMatcher::restoreMatchedPg(srcPg, srcPg.length(), destPg, pgMapOffSrc, pgMapLenSrc, revComplMatching, plainTextReadMode, false);
    }

    string
    SimplePgMatcher::restoreAutoMatchedPg(const string &pgPrefix, bool revComplMatching) {
        PseudoGenomeHeader *pgh = nullptr;
        ReadsSetProperties *prop = nullptr;
        bool plainTextReadMode = false;
        SeparatedPseudoGenomeBase::getPseudoGenomeProperties(pgPrefix, pgh, prop, plainTextReadMode);
        string destPg = SeparatedPseudoGenomePersistence::loadPseudoGenomeSequence(pgPrefix);
        ifstream pgMapOffSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_OFFSETS_FILE_SUFFIX);
        ifstream pgMapLenSrc = SeparatedPseudoGenomePersistence::getPseudoGenomeElementSrc(pgPrefix, SeparatedPseudoGenomeBase::PSEUDOGENOME_MAPPING_LENGTHS_FILE_SUFFIX);
        string resPg;
        resPg = SimplePgMatcher::restoreMatchedPg(resPg, pgh->getPseudoGenomeLength(), destPg, pgMapOffSrc, pgMapLenSrc,
                                                  revComplMatching, plainTextReadMode, true);
        delete (pgh);
        delete (prop);
        return resPg;
    }
}