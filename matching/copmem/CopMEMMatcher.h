#ifndef PGTOOLS_COPMEMMATCHER_H
#define PGTOOLS_COPMEMMATCHER_H

#include "../TextMatchers.h"

using namespace PgTools;

enum verbosity { v0, v1, v2 };
enum reverseMode { no, yes, both };

typedef std::pair<std::string, size_t> SequenceItem;
typedef std::vector<SequenceItem> SequenceVector;

template<class MyUINT1, class MyUINT2>
using HashBuffer =  std::pair<MyUINT1*, MyUINT2* >;

class CopMEMMatcher: public TextMatcher {
private:
    const string &srcText;
    const uint32_t targetMatchLength;
    int bigRef;

    verbosity isVerbose;
    reverseMode isRC;
    bool isFast;

    int K, H, L, k1, k2;
    std::uint32_t(*hashFunc)(const char*);
    std::uint32_t(*hashFuncMatrix[64][6])(const char*);

    std::string R_FN;
    std::string Q_FN;

    void initHashFuncMatrix();
    void initGlobals();
    void initParams();
    void calcCoprimes();
    void displayParams();

    template <class MyUINT>
    void genCumm(size_t N, const char* gen, MyUINT* cumm);

    void dumpMEM(SequenceItem& item1, SequenceItem& item2, size_t* match);
    void dumpMEMTight(SequenceItem& item1, size_t* match, size_t counter);

    std::pair<std::uint64_t*, std::uint64_t*> buffer2;
    std::pair<std::uint64_t*, std::uint32_t*> buffer1;
    std::pair<std::uint32_t*, std::uint32_t*> buffer0;

    template<typename MyUINT1, typename MyUINT2>
    HashBuffer<MyUINT1, MyUINT2> processRef();

    template <class MyUINT1, class MyUINT2>
    void deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf);

    template<typename MyUINT1, typename MyUINT2>
    void processQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector <TextMatch> &resMatches, const string &destText,
            bool destIsSrc, bool revComplMatching, uint32_t minMatchLength);

public:
    CopMEMMatcher(const string &srcText, const uint32_t targetMatchLength);

    virtual ~CopMEMMatcher();

    void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                    uint32_t minMatchLength) override;
};

#endif //PGTOOLS_COPMEMMATCHER_H
