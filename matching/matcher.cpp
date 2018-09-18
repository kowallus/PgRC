#include "matcher.h"

#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include "../readsset/PackedReadsSet.h"


void exactMatchConstantLengthPatterns(string text, string readsFile, ofstream& offsetsDest,
                                      ofstream& missedPatternsDest) {
    cout << "Reading reads set\n";
    PackedReadsSet* readsSet = PackedReadsSet::readReadsSet(readsFile);
    readsSet->printout();

    cout << "Feeding patterns...\n" << endl;
    ConstantLengthPatternsOnTextHashMatcher hashMatcher(readsSet->readLength(0));
    for(uint_reads_cnt_max i = 0; i < readsSet->readsCount(); i++)
        hashMatcher.addPattern(readsSet->getRead(i).data(), i);

    cout << "Matching...\n" << endl;
    hashMatcher.iterateOver(text.data(), text.length());

    for(int i = 0; i < 50; i++) {
        hashMatcher.moveNext();
        cout << "Matched: " << hashMatcher.getHashMatchPatternIndex() << "; " << hashMatcher.getHashMatchTextPosition() << endl;
        cout << readsSet->getRead(hashMatcher.getHashMatchPatternIndex()) << endl;
    }

    delete(readsSet);
}
