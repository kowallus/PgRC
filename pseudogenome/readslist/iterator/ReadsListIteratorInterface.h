/* 
 * File:   ReadsListIteratorInterface.h
 * Author: Tomek
 *
 * Created on 4 luty 2014, 12:06
 */

#ifndef READSLISTITERATORINTERFACE_H
#define	READSLISTITERATORINTERFACE_H

#include "../../../pgsaconfig.h"
#include "../../../utils/helper.h"
#include <vector>

using namespace std;

namespace PgSAIndex {

    template < typename uint_read_len, typename uint_reads_cnt, class ReadsIteratorClass >
    class ReadsIteratorInterface
    {
        public:
            
            virtual ~ReadsIteratorInterface() {};

            inline bool moveNext() { return static_cast<ReadsIteratorClass*>(this)->moveNextImpl(); };
            
            inline uint_reads_cnt getReadOriginalIndex() { return static_cast<ReadsIteratorClass*>(this)->getReadOriginalIndexImpl(); };
            inline uint_read_len getOccurrenceOffset() { return static_cast<ReadsIteratorClass*>(this)->getOccurrenceOffsetImpl(); };
            
            inline uint_flatten_occurrence_max getFlattenOccurrence() { return static_cast<ReadsIteratorClass*>(this)->getFlattenOccurrenceImpl(); };
            
            inline uint_reads_cnt getReadIndex() { return static_cast<ReadsIteratorClass*>(this)->getReadIndexImpl(); };
            
            //true if read does not contain duplicate kmers (of filter kmer length)
            inline bool hasDuplicateFilterFlag() { return static_cast<ReadsIteratorClass*>(this)->hasDuplicateFilterFlagImpl(); };

            inline bool hasOccurFlag() { return static_cast<ReadsIteratorClass*>(this)->hasOccurFlagImpl(); };
            inline void setOccurFlag() { static_cast<ReadsIteratorClass*>(this)->setOccurFlagImpl(); };
            inline uint_reads_cnt clearAllOccurFlags() { return static_cast<ReadsIteratorClass*>(this)->clearAllOccurFlagsImpl(); };

            inline bool hasOccurOnceFlag() { return static_cast<ReadsIteratorClass*>(this)->hasOccurOnceFlagImpl(); };
            inline void setOccurOnceFlagAndPushRead() { static_cast<ReadsIteratorClass*>(this)->setOccurOnceFlagAndPushReadImpl(); };
            inline void setOccurOnceFlagAndPushOccurrence() { static_cast<ReadsIteratorClass*>(this)->setOccurOnceFlagAndPushOccurrenceImpl(); };
            
            template<typename api_uint_reads_cnt>
            inline void clearAllOccurFlagsAndPushReadsWithSingleOccurrence(vector<api_uint_reads_cnt>& reads) { static_cast<ReadsIteratorClass*>(this)->template clearAllOccurFlagsAndPushReadsWithSingleOccurrenceImpl<api_uint_reads_cnt>(reads);  };
            
            inline uint_reads_cnt clearAllOccurFlagsAndCountSingleOccurrences() { return static_cast<ReadsIteratorClass*>(this)->clearAllOccurFlagsAndCountSingleOccurrencesImpl(); };
            
            template<typename api_uint_reads_cnt, typename api_uint_read_len>
            inline void clearAllOccurFlagsAndPushSingleOccurrences(vector<pair<api_uint_reads_cnt, api_uint_read_len>>& occurrences) { static_cast<ReadsIteratorClass*>(this)->template clearAllOccurFlagsAndPushSingleOccurrencesImpl<api_uint_reads_cnt, api_uint_read_len>(occurrences);  };
                        
            inline void clearAllOccurFlagsAndPushSingleOccurrencesFlatten(vector<uint_flatten_occurrence_max>& flattenOccurrences) { static_cast<ReadsIteratorClass*>(this)->clearAllOccurFlagsAndPushSingleOccurrencesFlattenImpl(flattenOccurrences);  };
            
    };
    
    template < typename uint_read_len, typename uint_reads_cnt, class ReadsListIteratorClass >
    class ReadsListIteratorInterface: public ReadsIteratorInterface<uint_read_len, uint_reads_cnt, ReadsListIteratorClass> {
        public:
            
            virtual ~ReadsListIteratorInterface() {};

            inline void initIteration(const uint_read_len& kmerLength) { static_cast<ReadsListIteratorClass*>(this)->initIterationImpl(kmerLength); };
            inline void setIterationPosition(const RPGOffset<uint_read_len, uint_reads_cnt>& saPos) { static_cast<ReadsListIteratorClass*>(this)->setIterationPositionImpl(saPos); };
            
            // for sparse SA
            inline void setIterationPosition(const RPGOffset<uint_read_len, uint_reads_cnt>& saPos, uchar shift) { static_cast<ReadsListIteratorClass*>(this)->setIterationPositionImpl(saPos, shift); };
  
    };

}


#endif	/* READSLISTITERATORINTERFACE_H */

