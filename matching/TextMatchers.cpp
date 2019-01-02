
#include <list>
#include "TextMatchers.h"

using namespace PgSAHelpers;

namespace PgTools {

    void backMatchExpand(const string &text, uint64_t &textPos, const string &pattern, uint64_t &patternPos, uint64_t &length) {
        const char *textGuardPtr = text.data();
        const char *textPtr = text.data() + textPos;
        const char *patternGuardPtr = pattern.data();
        const char *patternPtr = pattern.data() + patternPos;
        uint64_t i = 0;
        while (textPtr-- != textGuardPtr && patternPtr-- != patternGuardPtr && *textPtr == *patternPtr)
            i++;
        textPos -= i;
        patternPos -= i;
        length += i;
    }

    void forwardMatchExpand(const string &text, uint64_t textPos, const string &pattern, uint64_t patternPos, uint64_t &length) {
        const char *textGuardPtr = text.data() + text.length();
        const char *textPtr = text.data() + textPos + length;
        const char *patternGuardPtr = pattern.data() + pattern.length();
        const char *patternPtr = pattern.data() + patternPos + length;
        while (*textPtr == *patternPtr && textPtr++ != textGuardPtr && patternPtr++ != patternGuardPtr)
            length++;
    }

    DefaultTextMatcher::DefaultTextMatcher(const string &srcText, const uint32_t targetMatchLength) :
            srcText(srcText), targetMatchLength(targetMatchLength), matchingLength(targetMatchLength / 2), hashMatcher(matchingLength) {
        clock_checkpoint();

        cout << "Feeding pattern pseudogenome parts...\n" << endl;

        const char *srcTextPtr = srcText.data();
        uint32_t partsCount = srcText.length() / matchingLength;
        for (uint32_t i = 0; i < partsCount; i++)
            hashMatcher.addPattern(srcTextPtr + i * matchingLength, i);

        cout << "... finished in " << clock_millis() << " msec. " << endl;
    }

    void DefaultTextMatcher::matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                                        uint32_t minMatchLength) {
        clock_checkpoint();
        cout << "Matching...\n" << endl;

        const char *srcTextPtr = srcText.data();
        const char *destTextPtr = destText.data();
        hashMatcher.iterateOver(destTextPtr, destText.length());

        resMatches.clear();
        uint64_t furthestMatchEndPos = 0;
        list<TextMatch> currentMatches;

        int i = 0;
        uint64_t matchCount = 0;
        uint64_t matchCharsCount = 0;
        uint64_t matchCharsWithOverlapCount = 0;
        uint64_t shorterMatchCount = 0;
        uint64_t falseMatchCount = 0;
        while (hashMatcher.moveNext()) {
            uint64_t matchDestPos = hashMatcher.getHashMatchTextPosition();
            const uint32_t matchPatternIndex = hashMatcher.getHashMatchPatternIndex();
            uint64_t matchSrcPos = matchPatternIndex * matchingLength;
            if (destIsSrc && revComplMatching ? destText.length() - matchSrcPos < matchDestPos : matchDestPos >= matchSrcPos)
                continue;
            auto cmIt = currentMatches.begin();
            bool continueMatch = false;
            while (cmIt != currentMatches.end()) {
                if (matchDestPos + matchingLength > (*cmIt).endPosDestText()) {
                    currentMatches.erase(cmIt++);
                } else {
                    if ((*cmIt).posSrcText - (*cmIt).posDestText == matchSrcPos - matchDestPos) {
                        continueMatch = true;
                        break;
                    }
                    cmIt++;
                }
            }
            if (continueMatch)
                continue;

            bool confirmPatternMatch = strncmp(destTextPtr + matchDestPos, srcTextPtr + matchSrcPos, matchingLength) == 0;
            uint64_t matchLength = matchingLength;
            if (confirmPatternMatch) {
                backMatchExpand(destText, matchDestPos, srcText, matchSrcPos, matchLength);
                forwardMatchExpand(destText, matchDestPos, srcText, matchSrcPos, matchLength);
                if (matchLength >= minMatchLength) {
                    matchCount++;
                    matchCharsWithOverlapCount += matchLength;
                    matchCharsCount += matchLength - (matchDestPos < furthestMatchEndPos ?
                                                      furthestMatchEndPos - matchDestPos : 0);
                    const TextMatch &matchInfo = TextMatch(matchSrcPos, matchLength, matchDestPos);
                    if (furthestMatchEndPos < matchDestPos + matchLength)
                        furthestMatchEndPos = matchDestPos + matchLength;
                    resMatches.push_back(matchInfo);
                    currentMatches.push_back(matchInfo);
                } else
                    shorterMatchCount++;
            } else
                falseMatchCount++;

            /*            if ((i++ % 1000000) == 0 || matchLength > minMatchLength * 15) {
                             cout << (confirmPatternMatch?"Matched ":"False-matched ") << matchLength << " chars: <" << matchSrcPos << "; "
                                 << (matchSrcPos + matchLength) << ") " << " in " << matchDestPos << " (" << matchCharsCount << " chars matched total)" << endl;
                             cout << "Elapsed time: " << ((double) destPg.length() / matchDestPos) * clock_millis() / 1000 << "[s]" << endl;
                            if (matchLength < matchingLength * 10) {
                                uint8_t beforeChars = matchingLength / 4;
                                const uint64_t &minPos = min<uint64_t>(matchPatternIndex, matchDestPos);
                                beforeChars = minPos < beforeChars ? minPos : beforeChars;
                                uint8_t afterChars = matchingLength / 4;
                                const uint64_t &maxLength =
                                        min<uint64_t>(srcPg.length() - matchSrcPos, destPg.length() - matchDestPos) -
                                        matchLength;
                                afterChars = maxLength < afterChars ? maxLength : afterChars;
                                cout << srcPg.substr(matchSrcPos, matchLength) << endl;
                                cout << destPg.substr(matchDestPos, matchLength) << endl;
                                cout << srcPg.substr(matchSrcPos - beforeChars, beforeChars) << " ... " << srcPg.substr(matchSrcPos + matchLength, afterChars) << endl;
                                cout << destPg.substr(matchDestPos - beforeChars, beforeChars) << " ... " << destPg.substr(matchDestPos + matchLength, afterChars) << endl;
                            }
                        }*/
        }
        cout << "Finished matching in  " << clock_millis() << " msec. " << endl;
        cout << "Exact matched " << matchCount << " parts (too short matches: " << shorterMatchCount << "). False matches reported: " << falseMatchCount << "." << endl;
        cout << "Stage 1: Matched " << toString(matchCharsCount) << " (" << (toString((matchCharsCount * 100.0) / srcText.length(), 1)) << "%)" <<
             " Pg characters. Sum length of all matches is " << matchCharsWithOverlapCount << "." << endl;

        std::sort(resMatches.begin(), resMatches.end(), [](const TextMatch &match1, const TextMatch &match2) -> bool { return match1.length > match2.length; });

        cout << "Largest matches:" << endl;
        for (uint32_t i = 0; i < resMatches.size() && i < 10; i++)
            resMatches[i].report(cout);
    }
}
