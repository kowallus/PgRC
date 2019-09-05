#include <cassert>
#include "AbstractOverlapPseudoGenomeGenerator.h"

template<typename uint_read_len, typename uint_reads_cnt>
void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::removeCyclesAndPrepareComponents() {
    this->headRead = (uint_reads_cnt *) calloc(this->readsTotal() + 1, sizeof(uint_reads_cnt));
    uint_reads_cnt cyclesCount = 0;
    uint_reads_cnt overlapLost = 0;
    uint_reads_cnt nextIdx;
    for(uint_reads_cnt curIdx = 1; curIdx <= this->readsTotal(); curIdx++) {
        if (nextIdx = this->nextRead[curIdx]) {
            if (this->isHeadOf(curIdx, nextIdx)) {
                cyclesCount++;
                uint_read_len minOverlap = this->overlap[curIdx];
                uint_reads_cnt minOverlapIdx = curIdx;
                nextIdx = curIdx;
                while ((nextIdx = this->nextRead[nextIdx]) != curIdx) {
                    if (minOverlap > this->overlap[nextIdx]) {
                        uint_read_len minOverlap = this->overlap[nextIdx];
                        uint_reads_cnt minOverlapIdx = nextIdx;
                    }
                }
                overlapLost += minOverlap;
                uint_reads_cnt headIdx = this->nextRead[minOverlapIdx];
                this->nextRead[minOverlapIdx] = 0;
                this->overlap[minOverlapIdx] = 0;
                nextIdx = headIdx;
                this->headRead[headIdx] = 0;
                while ((nextIdx = this->nextRead[nextIdx]))
                    this->headRead[nextIdx] = headIdx;
            } else {
                if (this->headRead[curIdx] == 0)
                    this->headRead[nextIdx] = curIdx;
                else
                    this->headRead[nextIdx] = this->headRead[curIdx];
            }
        }
    }
    *logout << "Removed " << cyclesCount << " cycles (lost " << overlapLost << " symbols)" << endl;
}

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

// GENERATOR

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::init(bool pgGenerationMode) {
        nextRead = (uint_reads_cnt*) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        overlap = (uint_read_len *) calloc(readsTotal() + 1, sizeof(uint_read_len));
        if (isGenerationCyclesAware(pgGenerationMode))
            headRead = (uint_reads_cnt *) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        readsLeft = readsTotal();
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::dispose() {   
        free(nextRead);
        free(overlap);
        free(headRead);
    }


    template<typename uint_read_len, typename uint_reads_cnt>
    const vector<bool> AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::getBothSidesOverlappedReads(
            double overlappedReadsCountStopCoef) {

        clock_checkpoint();
        init(false);
        performOverlapping(overlappedReadsCountStopCoef, false);

        vector<uint_read_len> prevOverlap(readsTotal() + 1, 0);
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (hasSuccessor(i))
                prevOverlap[nextRead[i]] = overlap[i];

        uint_reads_cnt resCount = readsTotal();
        vector<bool> res(readsTotal(), true);
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++) {
            if (prevOverlap[i] && hasSuccessor(i))
                continue;
            if (hasSuccessor(i) && overlap[i] == readLength(i))
                continue;
            if (prevOverlap[i] == readLength(i))
                continue;
            res[i-1] = false;
            resCount--;
        }
        dispose();

        cout << "Found " << resCount << " both-side overlapped reads in " << clock_millis() << " msec\n";
        *logout << endl;

        return res;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::performOverlapping(
            double overlappedReadsCountStopCoef, bool pgGenerationMode) {
        clock_checkpoint();
        findOverlappingReads(overlappedReadsCountStopCoef, pgGenerationMode);
        if (pgGenerationMode) {
            pseudoGenomeLength = countPseudoGenomeLength();
            quick_stats();
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::unionOverlappedReads(
            uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght) {
        assert(headRead != 0);
        /* validation
        if (prevRead[nextIdx] > 0) 
            cout << curIdx << " cannot have successor " << nextIdx << " since it already success " << prevRead[nextIdx] << " by " << (int) overlap[prevRead[nextIdx]] << " symbols.\n";
        if (nextRead[curIdx] > 0)
            cout << nextIdx << " cannot have predecessor " << curIdx << " since it already precedes " << nextRead[curIdx] << " by " << (int) overlap[curIdx] << " symbols.\n";
        */
        setReadSuccessor(curIdx, nextIdx, overlapLenght);
        if (headRead[curIdx] == 0)
            headRead[nextIdx] = curIdx;
        else
            headRead[nextIdx] = headRead[curIdx];
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::setReadSuccessor(
            uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght) {
        nextRead[curIdx] = nextIdx;
        overlap[curIdx] = overlapLenght;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::getHead(uint_reads_cnt idx) {
        assert(headRead != 0);
        if (headRead[idx] == 0)
            return idx;

        headRead[idx] = getHead(headRead[idx]);
        return headRead[idx];
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_pg_len_max AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::countPseudoGenomeLength() {
        uint_pg_len_max len = 0;
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (overlap[i] < readLength(i))
                len += (readLength(i) - overlap[i]);

        return len;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::countComponents() {
        uint_reads_cnt count = 0;
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (!hasPredecessor(i) && hasSuccessor(i))
                count++;
        return count;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::countSingles() {
        uint_reads_cnt count = 0;
        for(uint_reads_cnt i = 1; i <= readsTotal(); i++)
            if (!hasPredecessor(i) && !hasSuccessor(i))
                count++;
        return count;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::quick_stats() {

        *logout << pseudoGenomeLength << " pseudo-genome length after overlapping\n";
        *logout << countComponents() << " pseudo-genome components\n";
        *logout << countSingles() << " single reads\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    template<class GeneratedPseudoGenome>
    GeneratedPseudoGenome* AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::assemblePseudoGenomeTemplate() {
        clock_checkpoint();
        if (!getReadsSetProperties()->constantReadLength) {
            cout << "ERROR: Unsuported variable reads length!";
            exit(EXIT_FAILURE);
        }

        GeneratedPseudoGenome* genPG =
                new GeneratedPseudoGenome(this->pseudoGenomeLength, getReadsSetProperties());

        for (uint_reads_cnt i = 1; i <= readsTotal(); i++) {
            uint_reads_cnt idx = i;
            if (!hasPredecessor(idx))
                do {
                    genPG->append(getReadUpToOverlap(idx), readLength(idx), overlap[idx], idx - 1);
                    idx = nextRead[idx];
                } while (idx != 0);
        }

        genPG->validate();
        *logout << "Pseudogenome assembled in " << clock_millis() << " msec\n\n";

        return genPG;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    bool AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::hasPredecessor(uint_reads_cnt incIdx) {
        assert(headRead != 0);
        return headRead[incIdx];
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    bool AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::hasSuccessor(uint_reads_cnt incIdx) {
        return nextRead[incIdx];
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    bool AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::isHeadOf(uint_reads_cnt idx, uint_reads_cnt head) {
        return getHead(idx) == head;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    bool AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::isPseudoGenomeLengthMaximalVirtual() {
        return isPGLengthMax(pseudoGenomeLength);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    bool AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::isPseudoGenomeLengthStandardVirtual() {
        return isPGLengthStd(pseudoGenomeLength);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    SeparatedPseudoGenome *
    AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::generateSeparatedPseudoGenome() {
        init();
        performOverlapping();

        SeparatedPseudoGenome* pg = assemblePseudoGenomeTemplate<GeneratedSeparatedPseudoGenome>();

        dispose();
        return pg;

    }

    template<typename uint_read_len, typename uint_reads_cnt>
    PseudoGenomeBase *
    AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::generatePseudoGenomeBase() {
        init();
        performOverlapping();

        PseudoGenomeBase* pg = 0;

        if (isPseudoGenomeLengthStandardVirtual())
            pg = assemblePseudoGenomeTemplate<GeneratedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len_std>>();
        else if (isPseudoGenomeLengthMaximalVirtual())
            pg = assemblePseudoGenomeTemplate<GeneratedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len_max>>();
        else
            cout << "Unsupported: pseudo genome length :(";

        dispose();
        return pg;
    }

    template class AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len_min, uint_reads_cnt_std>;
    template class AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len_std, uint_reads_cnt_std>;
    
}
