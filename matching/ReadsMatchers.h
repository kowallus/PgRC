#ifndef PGTOOLS_READSMATCHERS_H
#define PGTOOLS_READSMATCHERS_H

#include <iostream>
#include <fstream>
#include <string>

#include "../readsset/PackedConstantLengthReadsSet.h"
#include "ConstantLengthPatternsOnTextHashMatcher.h"
#include "../pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "copmem/CopMEMMatcher.h"

using namespace std;

namespace PgTools {

    static const int NOT_MATCHED_COUNT = UINT8_MAX;

    uint8_t
    countMismatches(const char *pattern, const char *text, uint64_t length, uint8_t maxMismatches = NOT_MATCHED_COUNT - 1);

    class AbstractReadsApproxMatcher;

    class DefaultReadsMatcher {
    protected:
        SeparatedPseudoGenome* sPg;
        const char* pgPtr;
        const uint_pg_len_max pgLength;
        bool revComplPg;
        PackedConstantLengthReadsSet *readsSet;
        uint32_t matchPrefixLength;

        vector<uint64_t> readMatchPos;
        vector<bool> readMatchRC;

        uint_reads_cnt_max matchedReadsCount = 0;

        const uint_reads_cnt_max readsCount;
        const uint_read_len_max readLength;
        const uint_read_len_max matchingLength;

        // stats
        uint64_t betterMatchCount = 0;
        uint64_t falseMatchCount = 0;

        virtual void initMatching();
        virtual void executeMatching(bool revCompMode = false) = 0;
        virtual void writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest) = 0;

        virtual SeparatedPseudoGenomeOutputBuilder *createSeparatedPseudoGenomeOutputBuilder(
                bool enableRevComp, bool enableMismatches) = 0;

        virtual void initEntryUpdating() = 0;
        virtual void updateEntry(DefaultReadsListEntry &entry, uint_reads_cnt_max matchIdx) = 0;
        virtual void closeEntryUpdating() = 0;

    public:
        static const uint_read_len_max DISABLED_PREFIX_MODE;
        static const uint64_t NOT_MATCHED_POSITION;

        DefaultReadsMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                            uint32_t matchPrefixLength);

        virtual ~DefaultReadsMatcher();

        void matchConstantLengthReads();

        static const string OFFSETS_SUFFIX;
        static const string SUFFIXES_SUFFIX;
        static const string MISSED_READS_SUFFIX;

        void writeMatchesInfo(const string &outPrefix);

        void writeIntoPseudoGenome(ostream& pgrcOut, uint8_t compressionLevel,
                const string &outPgPrefix, IndexesMapping* orgIndexesMapping);

        virtual const vector<bool> getMatchedReadsBitmap(uint8_t maxMismatches = NOT_MATCHED_COUNT - 1);

        virtual void transferMatchingResults(AbstractReadsApproxMatcher* approxMatcher);
    };

    class DefaultReadsExactMatcher: public DefaultReadsMatcher {
    private:
        DefaultConstantLengthPatternsOnTextHashMatcher* hashMatcher = 0;

    protected:
        void initMatching();
        void executeMatching(bool revCompMode = false);
        void writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest);

        SeparatedPseudoGenomeOutputBuilder *createSeparatedPseudoGenomeOutputBuilder(
                bool enableRevComp, bool enableMismatches) override;

        void initEntryUpdating() override {};
        void updateEntry(DefaultReadsListEntry &entry, uint_reads_cnt_max matchIdx) override {};
        void closeEntryUpdating() override {};
    public:
        DefaultReadsExactMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                                 uint32_t matchPrefixLength);

        virtual ~DefaultReadsExactMatcher();

        void transferMatchingResults(AbstractReadsApproxMatcher *approxMatcher) override;
    };

    class AbstractReadsApproxMatcher: public DefaultReadsMatcher {
    protected:
        uint8_t maxMismatches;
        uint8_t targetMismatches = 0;
        uint8_t minMismatches = 0;

        vector<uint8_t> readMismatchesCount;
        uint_reads_cnt_max matchedCountPerMismatches[NOT_MATCHED_COUNT + 1] = {};

        void writeMatchesInfo(ofstream &offsetsDest, ofstream &missedPatternsDest, ofstream &suffixesDest);
        void printApproxMatchingStats();

        SeparatedPseudoGenomeOutputBuilder *createSeparatedPseudoGenomeOutputBuilder(
                                                                                     bool enableRevComp, bool enableMismatches) override;
        void initEntryUpdating() override;
        string currentRead;
        void updateEntry(DefaultReadsListEntry &entry, uint_reads_cnt_max matchIdx) override;
        void closeEntryUpdating() override;

        virtual void initMatchingContinuation(DefaultReadsMatcher *pMatcher);

    public:
        AbstractReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                                   uint32_t matchPrefixLength, uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches = 0);

        virtual ~AbstractReadsApproxMatcher();

        void continueMatchingConstantLengthReads(DefaultReadsMatcher *pMatcher);

        void transferMatchingResults(AbstractReadsApproxMatcher* approxMatcher) override;

        const vector<bool> getMatchedReadsBitmap(uint8_t maxMismatches = NOT_MATCHED_COUNT - 1) override;

        friend class DefaultReadsExactMatcher;
    };

    class DefaultReadsApproxMatcher: public AbstractReadsApproxMatcher {
    private:
        DefaultConstantLengthPatternsOnTextHashMatcher* hashMatcher = 0;

    protected:
        uint_read_len_max partLength;

        void initMatching();

        void initMatchingContinuation(DefaultReadsMatcher *pMatcher) override;

        void executeMatching(bool revCompMode = false);

    public:
        DefaultReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                                  uint32_t matchPrefixLength, uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches = 0);

        virtual ~DefaultReadsApproxMatcher();
    };

    class InterleavedReadsApproxMatcher: public AbstractReadsApproxMatcher {
    private:
        InterleavedConstantLengthPatternsOnTextHashMatcher* hashMatcher = 0;

    protected:
        uint_read_len_max partLength;

        void initMatching();
        void executeMatching(bool revCompMode = false);

        void initMatchingContinuation(DefaultReadsMatcher *pMatcher) override;

    public:
        InterleavedReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                                  uint32_t matchPrefixLength, uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches = 0);

        virtual ~InterleavedReadsApproxMatcher();
    };

    class CopMEMReadsApproxMatcher: public AbstractReadsApproxMatcher {
    private:

    protected:
        uint_read_len_max partLength;

        void initMatching();
        void executeMatching(bool revCompMode = false);

        void initMatchingContinuation(DefaultReadsMatcher *pMatcher) override;

    public:
        CopMEMReadsApproxMatcher(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                                  uint32_t matchPrefixLength, uint16_t readsExactMatchingChars, uint8_t maxMismatches, uint8_t minMismatches = 0);

        virtual ~CopMEMReadsApproxMatcher();
    };

    const vector<bool> mapReadsIntoPg(SeparatedPseudoGenome* sPg, bool revComplPg, PackedConstantLengthReadsSet *readsSet,
                        uint_read_len_max matchPrefixLength, uint16_t preReadsExactMatchingChars,
                        uint16_t readsExactMatchingChars, uint16_t minCharsPerMismatch, char preMatchingMode,
                        char matchingMode, bool dumpInfo, ostream& pgrcOut, uint8_t compressionLevel,
                        const string &pgDestFilePrefix, IndexesMapping* orgIndexesMapping);
}

#endif //PGTOOLS_READSMATCHERS_H
