#include "DividedPCLReadsSets.h"

#include "iterator/DivisionReadsSetDecorators.h"

namespace PgTools {
    static const char *const DNA_SYMBOLS = "ACGT";

    static const char *const DNA_AND_N_SYMBOLS = "ACGNT";

    DividedPCLReadsSets::DividedPCLReadsSets(uint_read_len_max readLength, bool separateNReadsSet, bool nReadsLQ) :
                                       nReadsLQ(nReadsLQ), separateNReadsSet(separateNReadsSet) {
        if (separateNReadsSet || nReadsLQ)
            hqReadsSet = new PackedConstantLengthReadsSet(readLength, DNA_SYMBOLS, 4);
        else
            hqReadsSet = new PackedConstantLengthReadsSet(readLength, DNA_AND_N_SYMBOLS, 5);
        if (separateNReadsSet) {
            lqReadsSet = new PackedConstantLengthReadsSet(readLength, DNA_SYMBOLS, 4);
            nReadsSet = new PackedConstantLengthReadsSet(readLength, DNA_AND_N_SYMBOLS, 5);
        } else
            lqReadsSet = new PackedConstantLengthReadsSet(readLength, DNA_AND_N_SYMBOLS, 5);
    }

    DividedPCLReadsSets::~DividedPCLReadsSets() {
        disposeHqReadsSet();
        disposeLqReadsSet();
        disposeNReadsSet();
    }

    void DividedPCLReadsSets::disposeHqReadsSet() {
        if (hqReadsSet) {
            delete (hqReadsSet);
            hqReadsSet = 0;
        }
    }

    void DividedPCLReadsSets::disposeLqReadsSet() {
        if (lqReadsSet) {
            delete (lqReadsSet);
            lqReadsSet = 0;
        }
        if (lqMapping) {
            delete (lqMapping);
            lqMapping = 0;
        }
    }

    void DividedPCLReadsSets::disposeNReadsSet() {
        if (nReadsSet) {
            delete (nReadsSet);
            nReadsSet = 0;
        }
        if (nMapping) {
            delete (nMapping);
            nMapping = 0;
        }
    }

    DividedPCLReadsSets*
    DividedPCLReadsSets::getQualityDivisionBasedReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt,
                                                          uint_read_len_max readLength,
                                                          double error_limit, bool simplified_suffix_mode,
                                                          bool separateNReadsSet, bool nReadsLQ) {
        DividedPCLReadsSets* readsSets = new DividedPCLReadsSets(readLength, separateNReadsSet, nReadsLQ);
        time_checkpoint();
        QualityDividingReadsSetIterator<uint_read_len_max> *divReadsIt =
                new QualityDividingReadsSetIterator<uint_read_len_max>(readsIt, error_limit,
                                                                       simplified_suffix_mode);
        vector<uint_reads_cnt_max> lqMapping, nMapping;
        while (divReadsIt->moveNext()) {
            if (separateNReadsSet || nReadsLQ) {
                if (divReadsIt->containsN()) {
                    (separateNReadsSet ? readsSets->nReadsSet : readsSets->lqReadsSet)->addRead(
                            divReadsIt->getRead().data(), divReadsIt->getReadLength());
                    (separateNReadsSet ? nMapping : lqMapping).push_back(divReadsIt->getReadOriginalIndex());
                    continue;
                }
            }
            if (error_limit < 1 && !divReadsIt->isQualityHigh()) {
                readsSets->lqReadsSet->addRead(
                        divReadsIt->getRead().data(), divReadsIt->getReadLength());
                lqMapping.push_back(divReadsIt->getReadOriginalIndex());
            } else
            {
                readsSets->hqReadsSet->addRead(divReadsIt->getRead().data(), divReadsIt->getReadLength());
            }
        };
        cout << "Filtered " << (lqMapping.size() + nMapping.size());
        if (separateNReadsSet)
            cout << " (including " << nMapping.size() << " containing N)";
        cout << " reads (out of " <<
            (divReadsIt->getReadOriginalIndex()) << ") in " << time_millis() << " msec." << endl;
        *logout << endl;

        readsSets->lqMapping = new VectorMapping(std::move(lqMapping), divReadsIt->getReadOriginalIndex());
        if (separateNReadsSet)
            readsSets->nMapping = new VectorMapping(std::move(nMapping), divReadsIt->getReadOriginalIndex());

        delete(divReadsIt);
        return readsSets;
    }

    DividedPCLReadsSets *
    DividedPCLReadsSets::getSimpleDividedPCLReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt,
                                                      uint_read_len_max readLength, bool separateNReadsSet,
                                                      bool nReadsLQ) {
        if (separateNReadsSet || nReadsLQ)
            return DividedPCLReadsSets::getQualityDivisionBasedReadsSets(readsIt, readLength, 1, false, separateNReadsSet, nReadsLQ);
        DividedPCLReadsSets* readsSets = new DividedPCLReadsSets(readLength, false, false);
        while (readsIt->moveNext()) {
            readsSets->hqReadsSet->addRead(readsIt->getRead().data(), readsIt->getReadLength());
        }
        readsSets->lqMapping = new VectorMapping({}, readsSets->hqReadsSet->readsCount());
        return readsSets;
    }

    DividedPCLReadsSets *
    DividedPCLReadsSets::loadDivisionReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt,
                                               uint_read_len_max readLength, string lqDivisionFile, bool nReadsLQ,
                                               string nDivisionFile, bool skipHQReadsSet) {
        bool separateNReadsSet = nDivisionFile != "";
        DividedPCLReadsSets* readsSets = new DividedPCLReadsSets(readLength, separateNReadsSet, nReadsLQ);
        readsSets->lqMapping = VectorMapping::loadMapping(lqDivisionFile);
        if (separateNReadsSet)
            readsSets->nMapping = VectorMapping::loadMapping(nDivisionFile);

        uint_reads_cnt_max allCounter = 0;
        uint_reads_cnt_max lqCounter = 0;
        uint_reads_cnt_max nCounter = 0;
        while (readsIt->moveNext()) {
            if (readsSets->lqMapping->getReadOriginalIndex(lqCounter) == allCounter) {
                readsSets->lqReadsSet->addRead(readsIt->getRead().data(), readsIt->getReadLength());
                lqCounter++;
            } else if (separateNReadsSet &&
                    readsSets->nMapping->getReadOriginalIndex(nCounter) == allCounter) {
                readsSets->nReadsSet->addRead(readsIt->getRead().data(), readsIt->getReadLength());
                nCounter++;
            } else if (!skipHQReadsSet)
                readsSets->hqReadsSet->addRead(readsIt->getRead().data(), readsIt->getReadLength());
            allCounter++;
        }

        return readsSets;
    }

    void DividedPCLReadsSets::moveLqReadsFromHqReadsSetsToLqReadsSets(const vector<bool> &isReadHqInHqReadsSet) {
        if (nReadsLQ) {
            fprintf(stdout, "Unimplemented transferring reads between reads sets packed with different alphabet.\n");
            exit(EXIT_FAILURE);
        }
        uint_reads_cnt_max readsToMoveCount = 0;
        for(const bool hqFlag: isReadHqInHqReadsSet)
            if (!hqFlag) readsToMoveCount++;

        uint_reads_cnt_max newLqCounter = lqReadsSet->readsCount() + readsToMoveCount;
        uint_reads_cnt_max lqCounter = lqReadsSet->readsCount();
        bool ignoreLqSet = (lqCounter-- == 0);
        lqReadsSet->resize(newLqCounter);
        vector<uint_reads_cnt_max> &lqReadIdx = lqMapping->getMappingVector();
        lqReadIdx.resize(newLqCounter + 1);
        uint_reads_cnt_max allCounter = lqMapping->getReadsTotalCount();
        lqReadIdx[newLqCounter--] = allCounter;
        uint_reads_cnt_max nCounter = separateNReadsSet?nReadsSet->readsCount():0;
        bool ignoreNSet = (nCounter-- == 0);
        uint_reads_cnt_max hqCounter = hqReadsSet->readsCount();

        while(allCounter-- > 0 && hqCounter != 0) {
            if(!ignoreNSet) {
                if(nMapping->getReadOriginalIndex(nCounter) == allCounter) {
                    ignoreNSet = (nCounter-- == 0);
                    continue;
                }
            }
            if (!ignoreLqSet) {
                if(lqMapping->getReadOriginalIndex(lqCounter) == allCounter) {
                    lqReadsSet->copyRead(lqCounter, newLqCounter);
                    lqReadIdx[newLqCounter] = allCounter;
                    ignoreLqSet = (lqCounter-- == 0);
                    if (newLqCounter-- == 0)
                        break;
                    continue;
                }
            }
            if (!isReadHqInHqReadsSet[--hqCounter]) {
                lqReadsSet->copyPackedRead(hqReadsSet->getPackedRead(hqCounter), newLqCounter);
                lqReadIdx[newLqCounter] = allCounter;
                if (newLqCounter-- == 0)
                    break;
            }
        }

        uint_reads_cnt_max newHqCounter = 0;
        for(hqCounter = 0; hqCounter < hqReadsSet->readsCount(); hqCounter++) {
            if (isReadHqInHqReadsSet[hqCounter])
                hqReadsSet->copyRead(hqCounter, newHqCounter++);
        }
        hqReadsSet->resize(hqReadsSet->readsCount() - readsToMoveCount);
    }

    IndexesMapping* DividedPCLReadsSets::generateHqReadsIndexesMapping() {
        vector<uint_reads_cnt_max> hqReadIdx;
//        hqReadIdx.reserve(lqMapping->getReadsTotalCount() -
//            lqReadsSet->readsCount() - (separateNReadsSet?nReadsSet->readsCount():0));
        int64_t allCounter = -1;
        uint_reads_cnt_max lqCounter = 0;
        uint_reads_cnt_max nCounter = 0;
        while (++allCounter < lqMapping->getReadsTotalCount()) {
            if (lqMapping->getReadOriginalIndex(lqCounter) == allCounter) {
                lqCounter++;
            } else if (separateNReadsSet &&
                       nMapping->getReadOriginalIndex(nCounter) == allCounter) {
                nCounter++;
            } else
                hqReadIdx.push_back(allCounter);
        }
        return new VectorMapping(std::move(hqReadIdx), lqMapping->getReadsTotalCount());
    }

    void DividedPCLReadsSets::removeReadsFromLqReadsSet(const vector<bool> &isLqReadMappedIntoHqPg) {
        vector<uint_reads_cnt_max> &lqReadIdx = lqMapping->getMappingVector();
        uint_reads_cnt_max newLqCounter = 0;
        for(uint_reads_cnt_max lqCounter = 0; lqCounter < lqReadsSet->readsCount(); lqCounter++) {
            if (!isLqReadMappedIntoHqPg[lqCounter]) {
                lqReadIdx[newLqCounter] = lqReadIdx[lqCounter];
                lqReadsSet->copyRead(lqCounter, newLqCounter++);
            }
        }
        lqReadsSet->resize(newLqCounter);
        lqReadIdx[newLqCounter++] = lqMapping->getReadsTotalCount();
        lqReadIdx.resize(newLqCounter);
    }


    void DividedPCLReadsSets::removeReadsFromNReadsSet(const vector<bool> &isReadMappedIntoHqPg,
            uint_reads_cnt_max nBegIdx) {
        vector<uint_reads_cnt_max> &nReadIdx = nMapping->getMappingVector();
        uint_reads_cnt_max newNCounter = 0;
        for(uint_reads_cnt_max nCounter = 0; nCounter < nReadsSet->readsCount(); nCounter++) {
            if (!isReadMappedIntoHqPg[nCounter + nBegIdx]) {
                nReadIdx[newNCounter] = nReadIdx[nCounter];
                nReadsSet->copyRead(nCounter, newNCounter++);
            }
        }
        nReadsSet->resize(newNCounter);
        nReadIdx[newNCounter++] = nMapping->getReadsTotalCount();
        nReadIdx.resize(newNCounter);
    }

}
