#ifndef PGTOOLS_TEXTMATCHERS_H
#define PGTOOLS_TEXTMATCHERS_H

#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include <vector>
#include "../utils/helper.h"

namespace PgTools {

    struct TextMatch {
        uint64_t posSrcText;
        uint64_t length;
        uint64_t posDestText;

        TextMatch(uint64_t posSrcText, uint64_t length, uint64_t posDestText) : posSrcText(posSrcText), length(length),
                                                                                posDestText(posDestText) {}

        uint64_t endPosSrcText() const {
            return posSrcText + length;
        }

        uint64_t endPosDestText() const {
            return posDestText + length;
        }

        void report(ostream &out) {
            out << length << ": <" << posSrcText << ", " << endPosSrcText() << ") in " << posDestText << endl;
        }
    };

    class TextMatcher {

    public:
        virtual void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                                uint32_t minMatchLength) = 0;

    };

    class DefaultTextMatcher : public TextMatcher {
    private:
        const string &srcText;
        const uint32_t targetMatchLength;
        const uint64_t matchingLength;
        DefaultConstantLengthPatternsOnTextHashMatcher hashMatcher;
    public:
        DefaultTextMatcher(const string &srcText, const uint32_t targetMatchLength);

    public:
        void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                        uint32_t minMatchLength);

    };
}

#endif //PGTOOLS_TEXTMATCHERS_H
