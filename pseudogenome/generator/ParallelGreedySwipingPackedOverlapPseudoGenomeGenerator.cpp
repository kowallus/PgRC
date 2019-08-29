#include "ParallelGreedySwipingPackedOverlapPseudoGenomeGenerator.h"
#include "../../readsset/persistance/ReadsSetPersistence.h"
#include "AbstractOverlapPseudoGenomeGenerator.h"

#include <parallel/algorithm>
#include <cassert>

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

// GENERATOR

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt>
    ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::ParallelGreedySwipingPackedOverlapGeneratorTemplate(
            PackedConstantLengthReadsSet* orgReadsSet, bool ownReadsSet):
        packedReadsSet(orgReadsSet), ownReadsSet(ownReadsSet)
    {
        if (!orgReadsSet->isReadLengthConstant())
            cout << "Unsupported: variable length reads :(";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::~ParallelGreedySwipingPackedOverlapGeneratorTemplate() {
        if (ownReadsSet && this->packedReadsSet)
            delete(this->packedReadsSet);
    }
    
        template<typename uint_read_len, typename uint_reads_cnt>
    string ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::getReadUpToOverlap(uint_reads_cnt incIdx) {
        return packedReadsSet->getReadPrefix(incIdx - 1, this->overlap[incIdx]);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    uint_read_len ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::readLength(uint_reads_cnt incIdx) {
        return packedReadsSet->readLength(incIdx - 1);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::readsTotal() {
        return packedReadsSet->readsCount();
    }  
    template<typename uint_read_len, typename uint_reads_cnt>
    ReadsSetProperties* ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::getReadsSetProperties() {
        return packedReadsSet->getReadsSetProperties();
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    int ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::compareReads(uint_reads_cnt lIncIdx, uint_reads_cnt rIncIdx) {
        return packedReadsSet->comparePackedReads(lIncIdx - 1, rIncIdx - 1);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uchar ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::getSymbolOrderFromRead(uint_reads_cnt incIdx, uint_read_len offset) {
        char tmp;
        packedReadsSet->sPacker->reverseSequence(packedReadsSet->getPackedRead(incIdx - 1), offset, 1, &tmp);
        return (uchar) getReadsSetProperties()->symbolOrder[(uchar) tmp];
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    int ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::compareSuffixes(
            uchar lSymOrder, uchar rSymOrder, uint_read_len offset, uint_reads_cnt* ssiSymbolIdx) {
        uint_reads_cnt lIncIdx = sortedSuffixIdxsPtr[ssiSymbolIdx[lSymOrder]];
        uint_reads_cnt rIncIdx = sortedSuffixIdxsPtr[ssiSymbolIdx[rSymOrder]];
        return packedReadsSet->comparePackedReads(lIncIdx - 1, rIncIdx - 1, offset);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    int ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::compareSuffixWithPrefix(uint_reads_cnt sufIncIdx, uint_reads_cnt preIncIdx, uint_read_len sufOffset) {
        return packedReadsSet->compareSuffixWithPrefix(sufIncIdx - 1, preIncIdx - 1, sufOffset);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::updateSuffixQueue(
            uchar suffixGroup, uint_read_len suffixOffset,
            uint_reads_cnt* ssiSymbolIdx, uint_reads_cnt* ssiSymbolEnd, deque<uchar> &ssiOrder) {
        if (ssiSymbolIdx[suffixGroup] < ssiSymbolEnd[suffixGroup]) {
            deque<uchar>::reverse_iterator it = ssiOrder.rbegin();
            while (true) {
                if (it == ssiOrder.rend() || compareSuffixes(suffixGroup, *it, suffixOffset, ssiSymbolIdx) >= 0) {
                    ssiOrder.insert(it.base(), suffixGroup);
                    break;
                }
                it++;
            }
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    template<bool pgGenerationMode>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::setDuplicateSuccessor(
            uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght) {
        this->nextRead[curIdx] = nextIdx;
        this->overlap[curIdx] = overlapLenght;
        if (pgGenerationMode) {
            if (this->headRead[curIdx] == 0)
                this->headRead[nextIdx] = curIdx;
            else
                this->headRead[nextIdx] = this->headRead[curIdx];
        }
    }


    template<typename uint_read_len, typename uint_reads_cnt>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::prepareSortedReadsBlocks() {
        for (uint_reads_cnt i = 1; i <= packedReadsSet->readsCount(); i++)
            sortedReadsIdxs.push_back(i);
        PackedReadsComparator comparePacked = PackedReadsComparator(this);
        __gnu_parallel::sort(sortedReadsIdxs.begin(), sortedReadsIdxs.end(), comparePacked);

        const uint_symbols_cnt symbolsCount = packedReadsSet->getReadsSetProperties()->symbolsCount;
        blocksCount = pow(symbolsCount, blockPrefixLength);
        #pragma omp parallel for
        for (uint16_t b = 0; b < blocksCount; b++) {
            char *prefix = blockPrefixes[b];
            uint16_t val = b;
            for (uint8_t j = 0; j < blockPrefixLength; j++) {
                uint8_t div = val % symbolsCount;
                prefix[blockPrefixLength - j - 1] = packedReadsSet->getReadsSetProperties()->symbolsList[div];
                val = val / symbolsCount;
            }
            PackedReadVsPatternComparator comparePackedWithPattern = PackedReadVsPatternComparator(this, prefix);
            sortedReadsBlockPos[b] = lower_bound(sortedReadsIdxs.begin(), sortedReadsIdxs.end(),
                                                 comparePackedWithPattern.PATTERN_INDEX,
                                                 comparePackedWithPattern) - sortedReadsIdxs.begin();
        }
        sortedReadsBlockPos[blocksCount] = packedReadsSet->readsCount();
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    template<bool avoidCyclesMode>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::initAndFindDuplicates() {
        auto start_t = chrono::steady_clock::now();
        prepareSortedReadsBlocks();
        uint16_t b = 0;
        threadStartBlock[UINT8_MAX] = 0;
        for(uint8_t t = 1; t <= numberOfThreads; t++) {
            uint_reads_cnt threshold = ((double) t / numberOfThreads) * this->readsLeft;
            while (sortedReadsBlockPos[b] < threshold) {
                b++;
            }
            threadStartBlock[t] = b;
        }
        const uint_symbols_cnt symbolsCount = packedReadsSet->getReadsSetProperties()->symbolsCount;
        uint_reads_cnt sortedSuffixesLeftCount[MAX_BLOCKS_COUNT] = { 0 };
        uint_reads_cnt duplicatesCount = 0;
        #pragma omp parallel for reduction(+:sortedSuffixesLeftCount[0:MAX_BLOCKS_COUNT]) reduction(+:duplicatesCount)
        for(uint8_t t = 0; t < numberOfThreads; t++)
        {
            for (uint16_t b = threadStartBlock[t]; b < threadStartBlock[t + 1]; b++)
            {
                sortedReadsCount[b] = sortedReadsBlockPos[b + 1] - sortedReadsBlockPos[b];
                uchar curSymOrder = 0;
                sortedSuffixBlockPlusSymbolPos[b][curSymOrder] = sortedReadsBlockPos[b];
                if (sortedReadsCount[b]) {
                    uint16_t youngestNextSuffixBlock = (b % (blocksCount / symbolsCount)) * symbolsCount;
                    const auto &blockEnd = sortedReadsIdxs.begin() + sortedReadsBlockPos[b + 1];
                    for (auto srIt = sortedReadsIdxs.begin() + sortedReadsBlockPos[b]; srIt != blockEnd;) {
                        srIt++;
                        if (srIt != blockEnd && compareReads(*(srIt - 1), *srIt) == 0) {
                            if (avoidCyclesMode)
                                this->setDuplicateSuccessor<true>(*(srIt - 1), *srIt, packedReadsSet->maxReadLength());
                            else
                                this->setDuplicateSuccessor<false>(*(srIt - 1), *srIt, packedReadsSet->maxReadLength());
                            duplicatesCount++;
                        } else {
                            uchar nextSuffixSymbolOrder = getSymbolOrderFromRead(*(srIt - 1), blockPrefixLength);
                            sortedSuffixesLeftCount[youngestNextSuffixBlock + nextSuffixSymbolOrder]++;
                            while (curSymOrder != nextSuffixSymbolOrder)
                                sortedSuffixBlockPlusSymbolPos[b][++curSymOrder] = (srIt - 1) - sortedReadsIdxs.begin();
                        }
                    }
                }
                while (curSymOrder < symbolsCount)
                    sortedSuffixBlockPlusSymbolPos[b][++curSymOrder] = sortedReadsBlockPos[b + 1];
            }
        }
        this->readsLeft -= duplicatesCount;
        assert(accumulate(sortedSuffixesLeftCount, sortedSuffixesLeftCount + blocksCount, 0) == this->readsLeft);
        cout << "Found " << (readsTotal() - this->readsLeft) << " duplicates (..." << time_millis() << " msec)" << endl;
        sortedSuffixIdxs.resize(this->readsLeft, 0);
        mergeSortOfLeftSuffixes(1, sortedSuffixesLeftCount, sortedSuffixIdxs.data(), sortedReadsIdxs.data());

        #pragma omp parallel for schedule(guided)
        for(uint16_t b = 0; b < blocksCount; b++) {
            if (!sortedReadsCount[b])
                continue;
            uint_reads_cnt i = sortedReadsBlockPos[b];
            for (uint_reads_cnt j = i; j < sortedReadsBlockPos[b] + sortedReadsCount[b]; j++) {
                if (j == 0 || this->nextRead[sortedReadsIdxs[j - 1]] == 0)
                    sortedReadsIdxs[i++] = sortedReadsIdxs[j];
            }
            sortedReadsCount[b] = i - sortedReadsBlockPos[b];
        }
        assert(accumulate(sortedReadsCount, sortedReadsCount + blocksCount, 0) == this->readsLeft);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::mergeSortOfLeftSuffixes(
            uint_read_len offset, const uint_reads_cnt *sortedSuffixesLeftCount, uint_reads_cnt *sortedSuffixLeftIdxsPtr,
            const uint_reads_cnt *sortedSuffixIdxsPtr) {
        sortedSuffixBlockPos[0] = 0;
        for(uint16_t b = 1; b <= blocksCount; b++)
            sortedSuffixBlockPos[b] = sortedSuffixBlockPos[b - 1] + sortedSuffixesLeftCount[b - 1];

        const uint_symbols_cnt symbolsCount = packedReadsSet->getReadsSetProperties()->symbolsCount;
        this->sortedSuffixIdxsPtr = sortedSuffixIdxsPtr;
        #pragma omp parallel for schedule(guided)
        for(uint16_t b = 0; b < blocksCount; b++)
        {
            uint16_t prevYoungestBlock = b / symbolsCount;
            uint8_t lastPrefixSymbolOrder = b % symbolsCount;
            uint_reads_cnt ssiSymbolIdx[MAX_SYMBOLS_COUNT];
            uint_reads_cnt ssiSymbolEnd[MAX_SYMBOLS_COUNT];
            for (uint8_t j = 0; j < symbolsCount; j++) {
                ssiSymbolIdx[j] = sortedSuffixBlockPlusSymbolPos[prevYoungestBlock + (blocksCount / symbolsCount) * j]
                [lastPrefixSymbolOrder];
                ssiSymbolEnd[j] = sortedSuffixBlockPlusSymbolPos[prevYoungestBlock + (blocksCount / symbolsCount) * j]
                [lastPrefixSymbolOrder + 1];
            }
            deque<uchar> ssiOrder;
            for (uint8_t j = 0; j < symbolsCount; j++) {
                while (ssiSymbolIdx[j] < ssiSymbolEnd[j] &&
                       this->nextRead[sortedSuffixIdxsPtr[ssiSymbolIdx[j]]] != 0)
                    ssiSymbolIdx[j]++;
                updateSuffixQueue(j, offset + blockPrefixLength, ssiSymbolIdx, ssiSymbolEnd, ssiOrder);
            }
            uint_reads_cnt curPos = sortedSuffixBlockPos[b];
            while (!ssiOrder.empty()) {
                uchar j = ssiOrder.front();
                uint_reads_cnt sufIdx = sortedSuffixIdxsPtr[ssiSymbolIdx[j]];
                sortedSuffixLeftIdxsPtr[curPos++] = sufIdx;
                ssiOrder.pop_front();
                while (++ssiSymbolIdx[j] < ssiSymbolEnd[j] &&
                        this->nextRead[sortedSuffixIdxsPtr[ssiSymbolIdx[j]]] != 0);
                updateSuffixQueue(j, offset + blockPrefixLength, ssiSymbolIdx, ssiSymbolEnd, ssiOrder);
            }
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::findOverlappingReads(
            double overlappedReadsCountStopCoef, bool pgGenerationMode) {

        initAndFindDuplicates<false>();
        *logout << "Start overlapping.\n";

        uint_read_len overlapIterations = packedReadsSet->maxReadLength() * overlappedReadsCountStopCoef;
        if (overlapIterations > packedReadsSet->maxReadLength() - blockPrefixLength)
            overlapIterations = packedReadsSet->maxReadLength() - blockPrefixLength;
        for (int i = 1; i < overlapIterations; i++) {

            uint_reads_cnt sortedSuffixesLeftCount[MAX_BLOCKS_COUNT] = { 0 };
            overlapSortedReadsAndSuffixes<false>(i, sortedSuffixesLeftCount);
            vector<uint_reads_cnt> sortedSuffixLeft(this->readsLeft, 0);
            mergeSortOfLeftSuffixes(i + 1, sortedSuffixesLeftCount, sortedSuffixLeft.data(), sortedSuffixIdxs.data());
            sortedSuffixIdxs.swap(sortedSuffixLeft);

            *logout << this->readsLeft << " reads left after "
                    << (uint_read_len_max) (packedReadsSet->maxReadLength() - i) << " overlap (..." << time_millis() << ") msec" << endl;
        }

        sortedReadsIdxs.clear();
        sortedReadsIdxs.shrink_to_fit();
        sortedSuffixIdxs.clear();
        sortedSuffixIdxs.shrink_to_fit();

        if (pgGenerationMode) {
            this->removeCyclesAndPrepareComponents();
            *logout << this->countComponents() << " pseudo-genome components" << endl;
        }

        cout << "Overlapping done in " << time_millis() << " msec" << endl;
        *logout << endl;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    template<bool avoidCyclesMode>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::
            overlapSortedReadsAndSuffixes(uint_read_len suffixesOffset, uint_reads_cnt *sortedSuffixesLeftCount) {
        const uint_symbols_cnt symbolsCount = packedReadsSet->getReadsSetProperties()->symbolsCount;

        uint8_t threadsInIteration = suffixesOffset<25?2:(suffixesOffset<40?4:numberOfThreads);
        if (!avoidCyclesMode || threadsInIteration > numberOfThreads)
            threadsInIteration = numberOfThreads;

        uint16_t b = 0;
        uint_reads_cnt readsAndSuffixesCount = 0;
        threadStartBlock[UINT8_MAX] = 0;
        for(uint8_t t = 1; t <= threadsInIteration; t++) {
            uint_reads_cnt threshold = ((double) t / threadsInIteration) * this->readsLeft * 2;
            while (readsAndSuffixesCount < threshold) {
                readsAndSuffixesCount += sortedReadsCount[b] + sortedSuffixBlockPos[b + 1] - sortedSuffixBlockPos[b];
                b++;
            }
            threadStartBlock[t] = b;
        }

        uint_reads_cnt overlapsCount = 0;
#pragma omp parallel for reduction(+:sortedSuffixesLeftCount[0:MAX_BLOCKS_COUNT]) num_threads(threadsInIteration) \
                        reduction(+:overlapsCount)
        for(uint8_t t = 0; t < threadsInIteration; t++)
        {
            for (uint16_t b = threadStartBlock[t]; b < threadStartBlock[t + 1]; b++)
            {
                uchar curSymOrder = 0;
                sortedSuffixBlockPlusSymbolPos[b][curSymOrder] = sortedSuffixBlockPos[b];
                uint16_t youngestNextSuffixBlock = (b % (blocksCount / symbolsCount)) * symbolsCount;
                auto preIt = sortedReadsIdxs.begin() + sortedReadsBlockPos[b];
                auto sufIt = sortedSuffixIdxs.begin() + sortedSuffixBlockPos[b];
                const auto &preEnd = preIt + sortedReadsCount[b];
                const auto &sufEnd = sortedSuffixIdxs.begin() + sortedSuffixBlockPos[b + 1];
                sortedReadsCount[b] = 0;
                while (sufIt != sufEnd || preIt != preEnd) {
                    if (sufIt == sufEnd)
                        sortedReadsIdxs[sortedReadsBlockPos[b] + sortedReadsCount[b]++] = (*(preIt++));
                    else {
                        int cmpRes = -1;
                        auto curPreIt = preIt;
                        while (preIt != preEnd) {
                            if ((cmpRes = compareSuffixWithPrefix(*sufIt, *preIt, suffixesOffset)) != 0)
                                break;
                            if ((*sufIt != *preIt) && (!avoidCyclesMode || !this->isHeadOf(*sufIt, *preIt)))
                                break;
                            cmpRes = -1;
                            preIt++;
                        }

                        if (cmpRes)
                            preIt = curPreIt;
                        else {
                            uint_reads_cnt preIdx = *preIt;
                            while (preIt > curPreIt) {
                                *preIt = *(preIt - 1);
                                preIt--;
                            }
                            *preIt = preIdx;
                        }

                        if (cmpRes == 0) {
                            if (avoidCyclesMode) {
                                bool cycleCheck = false;
#pragma omp critical
                                {
                                    if (!(cycleCheck = this->isHeadOf(*sufIt, *preIt))) {
                                        if (this->headRead[*sufIt] == 0)
                                            this->headRead[*preIt] = *sufIt;
                                        else
                                            this->headRead[*preIt] = this->headRead[*sufIt];
                                    }
                                }
                                if (cycleCheck)
                                    continue;
                            }
                            this->setReadSuccessor(*sufIt, *preIt,
                                                       packedReadsSet->maxReadLength() - suffixesOffset);
                            overlapsCount++;
                            preIt++;
                        } else if (cmpRes > 0) {
                            sortedReadsIdxs[sortedReadsBlockPos[b] + sortedReadsCount[b]++] = (*(preIt++));
                            continue;
                        } else {
                            uchar nextSuffixSymbolOrder = getSymbolOrderFromRead(*sufIt, suffixesOffset +
                                                                                         blockPrefixLength);
                            sortedSuffixesLeftCount[youngestNextSuffixBlock + nextSuffixSymbolOrder]++;
                            while (curSymOrder != nextSuffixSymbolOrder)
                                sortedSuffixBlockPlusSymbolPos[b][++curSymOrder] = sufIt - sortedSuffixIdxs.begin();
                        }
                        sufIt++;
                    }
                }
                while (curSymOrder < symbolsCount)
                    sortedSuffixBlockPlusSymbolPos[b][++curSymOrder] = sortedSuffixBlockPos[b + 1];
            }
        }
        this->readsLeft -= overlapsCount;
        assert(accumulate(sortedSuffixesLeftCount, sortedSuffixesLeftCount + blocksCount, 0) == this->readsLeft);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>::validateSortedSuffixes(uint_read_len offset) const {
        for(uint_reads_cnt i = 1; i < this->readsLeft; i++) {
            uint_reads_cnt readA = sortedSuffixIdxs[i - 1] - 1;
            uint_reads_cnt readB = sortedSuffixIdxs[i] - 1;
            if (packedReadsSet->comparePackedReads(readA, readB, offset) > 0) {
                cout << i << "\t" << endl;
                cout << packedReadsSet->getRead(readA) << endl;
                cout << packedReadsSet->getRead(readB) << endl;
                exit(0);
            }
        }
    }

    // FACTORY

    template<typename uint_read_len, typename uint_reads_cnt>
    PseudoGenomeGeneratorBase* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getGeneratorFullTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet) {
        return new ParallelGreedySwipingPackedOverlapGeneratorTemplate<uint_read_len, uint_reads_cnt>(readsSet, ownReadsSet);
    }

    template<typename uint_read_len>
    PseudoGenomeGeneratorBase* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getGeneratorPartialTemplate(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet) {

        if (isReadsCountStd(readsSet->readsCount()))
            return getGeneratorFullTemplate<uint_read_len, uint_reads_cnt_std>(readsSet, ownReadsSet);
        else
            cout << "UNSUPPORTED READS COUNT!!!???";

        return 0;
    }

    PseudoGenomeGeneratorBase* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getGenerator(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator) {
        cout << "Reading reads set\n";
        // readsSet will be freed during generator destruction.
        PackedConstantLengthReadsSet *readsSet = PackedConstantLengthReadsSet::loadReadsSet(readsIterator);
        readsSet->printout();
        return this->getGenerator(readsSet, true);
    }

    PseudoGenomeGeneratorBase* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getGenerator(PackedConstantLengthReadsSet* readsSet, bool ownReadsSet) {
        PseudoGenomeGeneratorBase* generatorBase = 0;
        if (isReadLengthMin(readsSet->maxReadLength()))
            generatorBase = getGeneratorPartialTemplate<uint_read_len_min>(readsSet, ownReadsSet);
        else if (isReadLengthStd(readsSet->maxReadLength()))
            generatorBase = getGeneratorPartialTemplate<uint_read_len_std>(readsSet, ownReadsSet);
        else cout << "UNSUPPORTED READS LENGTH!!!";
        
        return generatorBase;
    }

    PseudoGenomeBase* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(
            ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator) {
        PseudoGenomeGeneratorFactory* pggf = new ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
        PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsIterator);
        PseudoGenomeBase* pgb = pggb->generatePseudoGenomeBase();
        delete(pggb);
        delete(pggf);

        if (pgb == 0) {
            fprintf(stderr, "Failed generating Pg\n");
            exit(EXIT_FAILURE);
        }

        return pgb;
    }

    PseudoGenomeBase* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generatePg(
            PackedConstantLengthReadsSet *readsSet) {
        ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory* pggf = new ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
        PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsSet, false);
        PseudoGenomeBase* pgb = pggb->generatePseudoGenomeBase();
        delete(pggb);
        delete(pggf);

        if (pgb == 0) {
            fprintf(stderr, "Failed generating Pg\n");
            exit(EXIT_FAILURE);
        }

        return pgb;
    }

    SeparatedPseudoGenome* ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::generateSeparatedPg(
            PackedConstantLengthReadsSet *readsSet) {
        ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory* pggf = new ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
        PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsSet, false);
        SeparatedPseudoGenome* pg = pggb->generateSeparatedPseudoGenome();
        delete(pggb);
        delete(pggf);

        if (pg == 0) {
            fprintf(stderr, "Failed generating Pg\n");
            exit(EXIT_FAILURE);
        }

        return pg;
    }

    const vector<bool> ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getHQReads(
            ReadsSourceIteratorTemplate<uint_read_len_max> *readsIterator, double qualityCoef) {
        PseudoGenomeGeneratorFactory* pggf = new ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
        PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsIterator);
        const vector<bool> res = pggb->getBothSidesOverlappedReads(qualityCoef);
        delete(pggb);
        delete(pggf);

        return res;
    }

    const vector<bool> ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory::getHQReads(
            PackedConstantLengthReadsSet *readsSet, double qualityCoef) {
        ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory* pggf = new ParallelGreedySwipingPackedOverlapPseudoGenomeGeneratorFactory();
        PseudoGenomeGeneratorBase* pggb = pggf->getGenerator(readsSet, false);
        const vector<bool> res = pggb->getBothSidesOverlappedReads(qualityCoef);
        delete(pggb);
        delete(pggf);

        return res;
    }


}
