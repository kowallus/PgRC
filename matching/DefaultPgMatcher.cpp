#include "DefaultPgMatcher.h"

#include "../pseudogenome/DefaultPseudoGenome.h"
#include "../pseudogenome/PackedPseudoGenome.h"
#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include <list>

namespace PgTools {

    DefaultPgMatcher::DefaultPgMatcher(const string& srcPgPrefix)
            :srcPgPrefix(srcPgPrefix) {

        PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(srcPgPrefix, pgh, plainTextReadMode);
    }

    DefaultPgMatcher::~DefaultPgMatcher() {
    }

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

    void backMatchExpand(const string& text, uint64_t& textPos, const string& pattern, uint64_t& patternPos, uint64_t& length) {
        const char* textGuardPtr = text.data();
        const char* textPtr = text.data() + textPos;
        const char* patternGuardPtr = pattern.data();
        const char* patternPtr = pattern.data() + patternPos;
        uint64_t i = 0;
        while (textPtr-- != textGuardPtr && patternPtr-- != patternGuardPtr && *textPtr == *patternPtr)
            i++;
        textPos -= i;
        patternPos -= i;
        length += i;
    }

    void forwardMatchExpand(const string& text, uint64_t textPos, const string& pattern, uint64_t patternPos, uint64_t& length) {
        const char* textGuardPtr = text.data() + text.length();
        const char* textPtr = text.data() + textPos + length;
        const char* patternGuardPtr = pattern.data() + pattern.length();
        const char* patternPtr = pattern.data() + patternPos + length;
        while (*textPtr == *patternPtr && textPtr++ != textGuardPtr && patternPtr++ != patternGuardPtr)
            length++;
    }

    void DefaultPgMatcher::exactMatchPg(const string& text,
            ofstream &offsetsDest, uint32_t minMatchLength, bool textFromSamePg, bool textIsRevComplOfPg) {
        clock_checkpoint();
        const uint_read_len_max readLength  = pgh->getMaxReadLength();
        const char* textPtr = text.data();

        cout << "Feeding pattern pseudogenome parts...\n" << endl;
        const uint_pg_len_max matchingLength = minMatchLength / 2;
        DefaultConstantLengthPatternsOnTextHashMatcher hashMatcher(matchingLength);
        string pgPattern = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(srcPgPrefix);
        const char* pgPtr = pgPattern.data();
        uint32_t partsCount = pgPattern.length() / matchingLength;
        for (uint32_t i = 0; i < partsCount; i++)
            hashMatcher.addPattern(pgPtr + i * matchingLength, i);

        cout << "... finished in " << clock_millis() << " msec. " << endl;
        clock_checkpoint();
        cout << "Matching...\n" << endl;
        hashMatcher.iterateOver(textPtr, text.length());

        vector<PgMatch> pgMatches;
        pgMatches.push_back(PgMatch(0,0,0));
        list<PgMatch> currentMatches;

        int i = 0;
        uint64_t matchCount = 0;
        uint64_t matchCharsCount = 0;
        uint64_t matchCharsWithOverlapCount = 0;
        uint64_t shorterMatchCount = 0;
        uint64_t falseMatchCount = 0;
        while (hashMatcher.moveNext()) {
            uint64_t matchTextPos = hashMatcher.getHashMatchTextPosition();
            const uint32_t matchPatternIndex = hashMatcher.getHashMatchPatternIndex();
            uint64_t matchPatternPos = matchPatternIndex * matchingLength;
            if(textFromSamePg && textIsRevComplOfPg?text.length() - matchPatternPos < matchTextPos:matchTextPos >= matchPatternPos)
                continue;
            auto cmIt = currentMatches.begin();
            bool continueMatch = false;
            while (cmIt != currentMatches.end()) {
                if (matchTextPos + matchingLength > (*cmIt).endPosText())  {
                    currentMatches.erase(cmIt++);
                } else {
                    if ((*cmIt).posPg - (*cmIt).posText == matchPatternPos - matchTextPos) {
                        continueMatch = true;
                        break;
                    }
                    cmIt++;
                }
            }
            if (continueMatch)
                continue;

            bool confirmPatternMatch = strncmp(textPtr + matchTextPos, pgPtr + matchPatternPos, matchingLength) == 0;
            uint64_t matchLength = matchingLength;
            if (confirmPatternMatch) {
                backMatchExpand(text, matchTextPos, pgPattern, matchPatternPos, matchLength);
                forwardMatchExpand(text, matchTextPos, pgPattern, matchPatternPos, matchLength);
                if (matchLength >= minMatchLength) {
                    matchCount++;
                    matchCharsWithOverlapCount += matchLength;
                    matchCharsCount += matchLength - (matchTextPos<pgMatches.back().endPosText()?pgMatches.back().endPosText() - matchTextPos:0);
                    const PgMatch &matchInfo = PgMatch(matchPatternPos, matchLength, matchTextPos);
                    pgMatches.push_back(matchInfo);
                    currentMatches.push_back(matchInfo);
                } else
                    shorterMatchCount++ ;
            } else
                falseMatchCount++;

            if ((i++ % 100000) == 0 || matchLength > minMatchLength * 15) {
                 cout << (confirmPatternMatch?"Matched ":"False-matched ") << matchLength << " chars: <" << matchPatternPos << "; "
                     << (matchPatternPos + matchLength) << ") " << " in " << matchTextPos << " (" << matchCharsCount << " chars matched total)" << endl;
                 cout << "Elapsed time: " << ((double) text.length() / matchTextPos) * clock_millis() / 1000 << "[s]" << endl;
                if (matchLength < matchingLength * 10) {
                    uint8_t beforeChars = matchingLength / 4;
                    const uint64_t &minPos = min<uint64_t>(matchPatternIndex, matchTextPos);
                    beforeChars = minPos < beforeChars ? minPos : beforeChars;
                    uint8_t afterChars = matchingLength / 4;
                    const uint64_t &maxLength =
                            min<uint64_t>(pgPattern.length() - matchPatternPos, text.length() - matchTextPos) -
                            matchLength;
                    afterChars = maxLength < afterChars ? maxLength : afterChars;
                    cout << pgPattern.substr(matchPatternPos, matchLength) << endl;
                    cout << text.substr(matchTextPos, matchLength) << endl;
                    cout << pgPattern.substr(matchPatternPos - beforeChars, beforeChars) << " ... " << pgPattern.substr(matchPatternPos + matchLength, afterChars) << endl;
                    cout << text.substr(matchTextPos - beforeChars, beforeChars) << " ... " << text.substr(matchTextPos + matchLength, afterChars) << endl;
                }
            }
        }
        cout << "Finished matching in  " << clock_millis() << " msec. " << endl;
        cout << "Exact matched " << matchCount << " parts (too short matches: " << shorterMatchCount << "). False matches reported: " << falseMatchCount << "." << endl;
        cout << "Matched " << matchCharsCount << " Pg characters. Sum length of all matches is " << matchCharsWithOverlapCount << "." << endl;

        std::sort(pgMatches.begin(), pgMatches.end(), [](const PgMatch& match1, const PgMatch& match2) -> bool
        { return match1.length > match2.length; });

        cout << "Largest matches:" << endl;
        for(uint32_t i = 0; i < pgMatches.size() && i < 10; i++)
            pgMatches[i].report(cout);
    }
}

