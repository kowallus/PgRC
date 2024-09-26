#include "ListOfConstantLengthReads.h"
#include "ReadsListTypes.h"

namespace PgIndex {

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::ListOfConstantLengthReads(uint_read_len readLength, uint_reads_cnt readsCount, uint_pg_len pseudoGenomeLength)
    : readLength(readLength) {
        this->readsCount = readsCount;
        this->pseudoGenomeLength = pseudoGenomeLength;
        this->pgReadsList = new uchar[(readsCount + 1) * (uint_max) LIST_ELEMENT_SIZE]();
        this->curRawIdx = 0;
        this->maxPos = 0;
        this->isSortRequired = false;
        this->duplicateFilterK = -1;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::ListOfConstantLengthReads(uint_read_len readLength, std::istream& src)
    : readLength(readLength) {
        src >> readsCount;
        src >> pseudoGenomeLength;
        int srchelper;
        src >> srchelper;
        duplicateFilterK = srchelper;
        src.get(); //"\n";
        this->pgReadsList = (uchar*) PgHelpers::readArray(src, sizeof (uchar) * (readsCount + 1) * (uint_max) LIST_ELEMENT_SIZE);

        this->isSortRequired = false;
        generateReverseIndex();
        this->pgReadsListEnd = this->pgReadsList + (uint_max) readsCount * LIST_ELEMENT_SIZE;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::~ListOfConstantLengthReads() {
        delete[] (this->pgReadsList);
        if (this->readsListIdx != 0)
            delete[] (this->readsListIdx);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::addImpl(uint_pg_len pos, uint_read_len len, uint_reads_cnt idx) {
        if (maxPos <= pos)
            maxPos = pos;
        else
            isSortRequired = true;
        *((uint_pg_len*) (this->pgReadsList + curRawIdx)) = pos;
        *((uint_reads_cnt*) (this->pgReadsList + curRawIdx + ORIGINAL_INDEX_OFFSET)) = idx;
        curRawIdx += LIST_ELEMENT_SIZE;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::validateImpl() {
        if (curRawIdx != readsCount * (uint_max) LIST_ELEMENT_SIZE)
            cout << "WARNING: Reads list generation failed: " << curRawIdx / LIST_ELEMENT_SIZE << " elements instead of " << readsCount << "\n";

        if (isSortRequired) {
            qsort(this->pgReadsList, readsCount, sizeof (uchar) * LIST_ELEMENT_SIZE, elementsCompare);
            isSortRequired = false;
        }

        // add guard
        *((uint_pg_len*) (this->pgReadsList + curRawIdx)) = this->pseudoGenomeLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::writeImpl(std::ostream& dest) {
        if (isSortRequired)
            cout << "WARNING: Reads list invalid!";
        dest << readsCount << "\n";
        dest << pseudoGenomeLength << "\n";
        dest << (int) duplicateFilterK << "\n";
        PgHelpers::writeArray(dest, this->pgReadsList, sizeof (uchar) * (readsCount + 1) * (uint_max) LIST_ELEMENT_SIZE);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::generateReverseIndex() {
        readsListIdx = new uint_reads_cnt[readsCount];
        for (uint_reads_cnt i = 0; i < readsCount; i++)
            readsListIdx[this->getReadOriginalIndexImpl(i)] = i;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::buildLUTImpl() {
        while (pseudoGenomeLength >> lookupStepShift++ > readsCount);
        
        lookup.resize((pseudoGenomeLength >> lookupStepShift) + 2);
//        cout << ((pseudoGenomeLength >> lookupStepShift) + 2) << " " << lookup.max_size() << "\n";
        uint_pg_len j = 0;
        uint_reads_cnt i = 0;
        while (i < readsCount) {

            uint_pg_len nextReadPos = getReadPositionImpl(i + 1);           
            while ((j << lookupStepShift) < nextReadPos)
                lookup[j++] = i;
            
            i++;
        }

        lookup[j] = readsCount;
    }
        
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    uint_reads_cnt ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::findFurthestReadContainingImpl(uint_pg_len pos) {
      
        uint_reads_cnt lIdx = lookup[pos >> lookupStepShift];
        uint_reads_cnt rIdx = lookup[(pos >> lookupStepShift) + 1] + 1;
        if (rIdx > readsCount)
            rIdx = readsCount;
        
        if ((getReadPositionImpl(lIdx) > pos) || (getReadPositionImpl(rIdx) <= pos))
            cout << "whoops\n";
        
        uint_reads_cnt mIdx;
        
        long long int cmpRes = 0;
        
        while (lIdx < rIdx) {
            mIdx = (lIdx + rIdx) / 2;
            cmpRes = (long long int) pos - (long long int) getReadPositionImpl(mIdx);

            if (cmpRes >= 0) 
                lIdx = mIdx + 1;
            else if (cmpRes < 0)
                rIdx = mIdx;
        }
 
        return lIdx - 1;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    inline uint_reads_cnt ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::clearAllOccurFlagsAndCountSingleOccurrencesImpl() {
        uint_reads_cnt count = 0;
        for (typename vector<uchar*>::iterator it = occurFlagsReadsList.begin(); it != occurFlagsReadsList.end(); ++it) {
            if (hasOccurOnceFlagByAddress(*it))
                count++;
            clearOccurFlagsByAddress(*it);
        }

        occurFlagsReadsList.clear();

        return count;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    inline void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::clearAllOccurFlagsAndPushSingleOccurrencesFlattenImpl(vector<uint_flatten_occurrence_max>& flattenOccurrences) {
        for (typename vector<pair < uchar*, uint_read_len>>::iterator it = occurFlagsOccurrencesList.begin(); it != occurFlagsOccurrencesList.end(); ++it) {
            if (this->hasOccurOnceFlagByAddress((*it).first))
                flattenOccurrences.push_back(((uint_flatten_occurrence_max) this->getReadOriginalIndexByAddress((*it).first)) * this->getMaxReadLength() + (*it).second);
            this->clearOccurFlagsByAddress((*it).first);
        }

        occurFlagsOccurrencesList.clear();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    inline uint_reads_cnt ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::clearAllOccurFlagsImpl() {
        uint_reads_cnt readsWithOccurFlagCount = occurFlagsReadsList.size();
        for (typename vector<uchar*>::iterator it = occurFlagsReadsList.begin(); it != occurFlagsReadsList.end(); ++it)
            clearOccurFlagsByAddress(*it);

        occurFlagsReadsList.clear();

        return readsWithOccurFlagCount;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    inline void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::clearOccurFlagsByAddress(uchar* address) {
        flagsVariableByAddress(address) &= ~(OCCURFLAGS_MASK);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, unsigned char LIST_ELEMENT_SIZE, uchar FLAGS_OFFSET>
    inline void ListOfConstantLengthReads<uint_read_len, uint_reads_cnt, uint_pg_len, LIST_ELEMENT_SIZE, FLAGS_OFFSET>::clearOccurFlagsImpl(uint_reads_cnt idx) {
        clearOccurFlagsByAddress(idx2Address(idx));
    }
    
    template class ListOfConstantLengthReads<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, 9, 4>;
    template class ListOfConstantLengthReads<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, 9, 4>;
    template class ListOfConstantLengthReads<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, 13, 8>;
    template class ListOfConstantLengthReads<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, 13, 8>;
        
}