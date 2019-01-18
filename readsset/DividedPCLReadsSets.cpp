#include "DividedPCLReadsSets.h"

#include "iterator/DivisionReadsSetDecorators.h"

namespace PgTools {
    DividedPCLReadsSets::DividedPCLReadsSets(uint_read_len_max readLength, bool separateNReadsSet, bool nReadsLQ) :
                                       nReadsLQ(nReadsLQ), separateNReadsSet(separateNReadsSet) {
        if (separateNReadsSet || nReadsLQ)
            hqReadsSet = new PackedConstantLengthReadsSet(readLength, "ACGT", 4);
        else
            hqReadsSet = new PackedConstantLengthReadsSet(readLength, "ACGTN", 5);
        if (separateNReadsSet) {
            lqReadsSet = new PackedConstantLengthReadsSet(readLength, "ACGT", 4);
            nReadsSet = new PackedConstantLengthReadsSet(readLength, "ACGTN", 5);
        } else
            lqReadsSet = new PackedConstantLengthReadsSet(readLength, "ACGTN", 5);
    }

    DividedPCLReadsSets::~DividedPCLReadsSets() {
        if (hqReadsSet)
            delete(hqReadsSet);
        if (lqReadsSet)
            delete(lqReadsSet);
        if (nReadsSet)
            delete(nReadsSet);
        if (lqMapping)
            delete(lqMapping);
        if (nMapping)
            delete(nMapping);
    }

    DividedPCLReadsSets*
    DividedPCLReadsSets::getQualityDivisionBasedReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt,
                                                          uint_read_len_max readLength,
                                                          double error_limit, bool separateNReadsSet, bool nReadsLQ) {
        DividedPCLReadsSets* readsSets = new DividedPCLReadsSets(readLength, separateNReadsSet, nReadsLQ);
        clock_checkpoint();
        QualityDividingReadsSetIterator<uint_read_len_max> *divReadsIt =
                new QualityDividingReadsSetIterator<uint_read_len_max>(readsIt, error_limit);
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
            if (!divReadsIt->isQualityHigh()) {
                readsSets->lqReadsSet->addRead(
                        divReadsIt->getRead().data(), divReadsIt->getReadLength());
                lqMapping.push_back(divReadsIt->getReadOriginalIndex());
            } else {
                readsSets->hqReadsSet->addRead(divReadsIt->getRead().data(), divReadsIt->getReadLength());
            }
        };
        cout << "Filtered " << (lqMapping.size() + nMapping.size());
        if (separateNReadsSet)
            cout << " (including " << nMapping.size() << " containing N)";
        cout << " reads (out of " <<
            (divReadsIt->getReadOriginalIndex()) << ") in " << clock_millis() << " msec." << endl << endl;

        readsSets->lqMapping = new VectorMapping(std::move(lqMapping), divReadsIt->getReadOriginalIndex());
        if (separateNReadsSet)
            readsSets->nMapping = new VectorMapping(std::move(nMapping), divReadsIt->getReadOriginalIndex());

        return readsSets;
    }

    DividedPCLReadsSets *
    DividedPCLReadsSets::getSimpleDividedPCLReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt,
                                                      uint_read_len_max readLength, bool separateNReadsSet,
                                                      bool nReadsLQ) {
        DividedPCLReadsSets* readsSets = new DividedPCLReadsSets(readLength, separateNReadsSet, nReadsLQ);
        while (readsIt->moveNext()) {
            readsSets->hqReadsSet->addRead(readsIt->getRead().data(), readsIt->getReadLength());
        }
        return readsSets;
    }

    DividedPCLReadsSets *
    DividedPCLReadsSets::loadDivisionReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt,
                                               uint_read_len_max readLength, string lqDivisionFile, bool nReadsLQ,
                                               string nDivisionFile) {
        bool separateNReadsSet = nDivisionFile != "";
        DividedPCLReadsSets* readsSets = new DividedPCLReadsSets(readLength, separateNReadsSet, nReadsLQ);
        readsSets->lqMapping = VectorMapping::loadMapping(lqDivisionFile);
        if (separateNReadsSet)
            readsSets->nMapping = VectorMapping::loadMapping(nDivisionFile);

        uint_reads_cnt_max allCounter = 0;
        uint_reads_cnt_max lqCounter = 0;
        uint_reads_cnt_max nCounter = 0;
        while (readsIt->moveNext()) {
            PackedConstantLengthReadsSet* targetSet = readsSets->hqReadsSet;
            if (readsSets->lqMapping->getReadOriginalIndex(lqCounter) == allCounter) {
                targetSet = readsSets->lqReadsSet;
                lqCounter++;
            } else if (separateNReadsSet &&
                    readsSets->nMapping->getReadOriginalIndex(nCounter) == allCounter) {
                targetSet = readsSets->nReadsSet;
                nCounter++;
            }
            targetSet->addRead(readsIt->getRead().data(), readsIt->getReadLength());
            allCounter++;
        }

        return readsSets;
    }
}
