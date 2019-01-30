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

class processApproxMatchQueryTight;

class CopMEMMatcher: public TextMatcher {
private:
    const string &srcText;
    const char* start1;
    size_t N;
    int bigRef;
    const int L;
    int K, k1, k2;
    std::uint32_t(*hashFunc)(const char*);
    std::uint32_t(*hashFuncMatrix[64][6])(const char*);
    const int H = 1;

    int LK2 = (L - K) / 2;
    int LK2_MINUS_4 = LK2 - 4;
    int K_PLUS_LK24 = K + LK2_MINUS_4;

    void initHashFuncMatrix();
    void initParams(uint32_t minMatchLength);
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
    void processExactMatchQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches,
                                     const string &destText,
                                     bool destIsSrc, bool revComplMatching, uint32_t minMatchLength);

    template<typename MyUINT1, typename MyUINT2>
    uint64_t processApproxMatchQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, const char *pattern, const uint_read_len_max length,
                                 uint8_t maxMismatches, uint8_t minMismatches, uint8_t &mismatchesCount,
                                 uint64_t& multiMatchCount, uint64_t& falseMatchCount);

public:
    CopMEMMatcher(const string &srcText, const uint32_t targetMatchLength, uint32_t minMatchLength = UINT32_MAX);

    virtual ~CopMEMMatcher();

    void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                    uint32_t minMatchLength) override;

    uint64_t approxMatchPattern(const char *pattern, const uint_read_len_max length, uint8_t maxMismatches, uint8_t minMismatches,
            uint8_t &mismatchesCount, uint64_t& multiMatchCount, uint64_t& falseMatchCount);

};

#endif //PGTOOLS_COPMEMMATCHER_H
