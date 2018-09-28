#include "AbstractOverlapPseudoGenomeGenerator.h"

using namespace PgSAReadsSet;
using namespace PgSAHelpers;

// GENERATOR

namespace PgSAIndex {

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::init() {
        nextRead = (uint_reads_cnt*) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        overlap = (uint_read_len*) calloc(readsTotal() + 1, sizeof(uint_read_len));
        headRead = (uint_reads_cnt*) calloc(readsTotal() + 1, sizeof(uint_reads_cnt));
        readsLeft = readsTotal();
    }
    
    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::dispose() {   
        free(nextRead);
        free(overlap);
        free(headRead);
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    PseudoGenomeBase* AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::generatePseudoGenomeBase() {

        clock_checkpoint();
      
        init();  
        
        findOverlappingReads();
        
        pseudoGenomeLength = countPseudoGenomeLength();
        quick_stats();

        PseudoGenomeBase* pgb = 0;
        
        if (isPseudoGenomeLengthStandardVirtual())
            pgb = assemblePseudoGenomeTemplate<uint_pg_len_std>();
        else if (isPseudoGenomeLengthMaximalVirtual())
            pgb = assemblePseudoGenomeTemplate<uint_pg_len_max>();
        else
            cout << "Unsupported: pseudo genome length :(";
        
        dispose();
        
        return pgb;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    void AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::setReadSuccessor(uint_reads_cnt curIdx, uint_reads_cnt nextIdx, uint_read_len overlapLenght) {
        
        /* validation
        if (prevRead[nextIdx] > 0) 
            cout << curIdx << " cannot have successor " << nextIdx << " since it already success " << prevRead[nextIdx] << " by " << (int) overlap[prevRead[nextIdx]] << " symbols.\n";
        if (nextRead[curIdx] > 0)
            cout << nextIdx << " cannot have predecessor " << curIdx << " since it already precedes " << nextRead[curIdx] << " by " << (int) overlap[curIdx] << " symbols.\n";
        */
        
        nextRead[curIdx] = nextIdx;
        overlap[curIdx] = overlapLenght;
        if (headRead[curIdx] == 0)
            headRead[nextIdx] = curIdx;
        else
            headRead[nextIdx] = headRead[curIdx];
        readsLeft--;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    uint_reads_cnt AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::getHead(uint_reads_cnt idx) {
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

        cout << pseudoGenomeLength << " bytes after overlapping\n";
        cout << countComponents() << " pseudo-genome components\n";
        cout << countSingles() << " single reads\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    template<typename uint_pg_len, class GeneratedPseudoGenome>
    PseudoGenomeBase* AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::assemblePseudoGenomeFullTemplate() {
        clock_checkpoint();

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
        cout << "Pseudogenome assembled in " << clock_millis() << " msec\n\n";

        return genPG;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    template<typename uint_pg_len>
    PseudoGenomeBase* AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::assemblePseudoGenomeTemplate() {
        if (getReadsSetProperties()->constantReadLength)
            return assemblePseudoGenomeFullTemplate<uint_pg_len, GeneratedPseudoGenomeOfConstantLengthReadsType < uint_read_len, uint_reads_cnt, uint_pg_len >> ();

        cout << "ERROR: Unsuported variable reads length!";
        return 0;
    }

    template<typename uint_read_len, typename uint_reads_cnt>
    bool AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len, uint_reads_cnt>::hasPredecessor(uint_reads_cnt incIdx) {
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

    template class AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len_min, uint_reads_cnt_std>;
    template class AbstractOverlapPseudoGenomeGeneratorTemplate<uint_read_len_std, uint_reads_cnt_std>;
    
}
