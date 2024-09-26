#ifndef READSLISTINTERFACE_H_INCLUDED
#define READSLISTINTERFACE_H_INCLUDED

#define DUPINREADS_FLAG 0x10
#define OCCUR_FLAG 0x01
#define OCCURFLAGS_MASK 0x03

#include "../../pgrc/pg-config.h"
#include "iterator/ReadsListIteratorInterface.h"

namespace PgIndex {
        
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    class ReadsListInterface
    {
        public:

            virtual ~ReadsListInterface() {};
            
            inline uint_reads_cnt getReadsCount() { return static_cast<ReadsListClass*>(this)->getReadsCountImpl(); };
            inline uint_read_len getMaxReadLength() { return static_cast<ReadsListClass*>(this)->getMaxReadLengthImpl(); };
            
            inline uint_pg_len getReadPosition(uint_reads_cnt idx) { return static_cast<ReadsListClass*>(this)->getReadPositionImpl(idx); };
            inline uint_read_len getReadLength(uint_reads_cnt idx) { return static_cast<ReadsListClass*>(this)->getReadLengthImpl(idx); };
            inline uint_reads_cnt getReadOriginalIndex(uint_reads_cnt idx) { return static_cast<ReadsListClass*>(this)->getReadOriginalIndexImpl(idx); };

            inline uint_reads_cnt getReadsListIndexOfOriginalIndex(uint_reads_cnt originalIdx) { return static_cast<ReadsListClass*>(this)->getReadsListIndexOfOriginalIndexImpl(originalIdx); };
            
            inline void buildLUT() { static_cast<ReadsListClass*>(this)->buildLUTImpl(); };
            inline uint_reads_cnt findFurthestReadContaining(uint_pg_len pos) { return static_cast<ReadsListClass*>(this)->findFurthestReadContainingImpl(pos); };
            
            inline uint_read_len getDuplicateFilterKmerLength() { return static_cast<ReadsListClass*>(this)->getDuplicateFilterKmerLengthImpl(); };
            inline void setDuplicateFilterKmerLength(uint_read_len kLength) { static_cast<ReadsListClass*>(this)->setDuplicateFilterKmerLengthImpl(kLength); };
            
            //true if read does not contain duplicate kmers (of filter kmer length)
            inline bool hasDuplicateFilterFlag(uint_reads_cnt idx) { return static_cast<ReadsListClass*>(this)->hasDuplicateFilterFlagImpl(idx); };
            inline void setDuplicateFilterFlag(uint_reads_cnt idx) { static_cast<ReadsListClass*>(this)->setDuplicateFilterFlagImpl(idx); };

            inline bool hasOccurFlag(uint_reads_cnt idx) { return static_cast<ReadsListClass*>(this)->hasOccurFlagImpl(idx); };
            inline bool hasOccurOnceFlag(uint_reads_cnt idx) { return static_cast<ReadsListClass*>(this)->hasOccurOnceFlagImpl(idx); };

            inline void setOccurFlag(uint_reads_cnt idx) { static_cast<ReadsListClass*>(this)->setOccurFlagImpl(idx); };
            inline void setOccurOnceFlag(uint_reads_cnt idx) { static_cast<ReadsListClass*>(this)->setOccurOnceFlagImpl(idx); };

            inline void clearOccurFlags(uint_reads_cnt idx) { static_cast<ReadsListClass*>(this)->clearOccurFlagsImpl(idx); };

            inline void write(ostream& dest) { static_cast<ReadsListClass*>(this)->writeImpl(dest); };
            
            inline const uchar getListElementSize() { return static_cast<ReadsListClass*>(this)->getListElementSizeImpl(); };
    };

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class GeneratedReadsListInterface: public ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>
    {
        public:

            virtual ~GeneratedReadsListInterface() {};

            virtual void add(uint_pg_len pos, uint_read_len len, uint_reads_cnt idx) { static_cast<ReadsListClass*>(this)->addImpl(pos, len, idx); };
            virtual void validate() { static_cast<ReadsListClass*>(this)->validateImpl(); };

    };

    template <class ReadsListClass> 
    struct ReadsListIteratorFactoryTemplate{
        typedef int ReadsListIteratorClass; // wrong type !!!
        static ReadsListIteratorClass* getReadsListIterator() {
            return 0;
        };
    };
    
}

#endif // READSLISTINTERFACE_H_INCLUDED
