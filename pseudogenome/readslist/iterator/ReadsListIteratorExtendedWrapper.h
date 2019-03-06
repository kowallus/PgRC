#ifndef PGTOOLS_READSLISTITERATOREXTENDEDWRAPPER_H
#define PGTOOLS_READSLISTITERATOREXTENDEDWRAPPER_H

#include "../ReadsListInterface.h"

#include "ExtendedReadsListIteratorInterface.h"
#include "../../../readsset/persistance/ReadsSetPersistence.h"

namespace PgTools {

    using namespace PgSAIndex;

    class ReadsListIteratorExtendedWrapperBase: public DefaultReadsListIteratorInterface {
    public:
        virtual void applyIndexesMapping(IndexesMapping* orgIndexesMapping) = 0;

        virtual void applyRevComplPairFileFlag() = 0;
    };

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class ReadsListIteratorExtendedWrapper: public ReadsListIteratorExtendedWrapperBase {
    private:
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList;
        uint_reads_cnt currentIdx = 0;
        DefaultReadsListEntry entry;
        bool revComplPairFile = false;
        IndexesMapping* orgIndexesMapping;;

    public:
        ReadsListIteratorExtendedWrapper(ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList)
        : readsList(readsList), orgIndexesMapping(new DirectMapping(readsList->getReadsCount()))  {}

        ~ReadsListIteratorExtendedWrapper() { }

        inline uint_reads_cnt_max mapIndex(uint_reads_cnt idx) {
            return orgIndexesMapping->getReadOriginalIndex(idx);
        }

        void applyIndexesMapping(IndexesMapping* orgIndexesMapping) {
            delete(this->orgIndexesMapping);
            this->orgIndexesMapping = orgIndexesMapping;
        }

        void applyRevComplPairFileFlag() {
            revComplPairFile = true;
        }

        bool moveNext() override {
            const uint_reads_cnt_max orgIdx = mapIndex(readsList->getReadOriginalIndex(currentIdx));
            entry.advanceEntryByPosition(readsList->getReadPosition(currentIdx), orgIdx, revComplPairFile?orgIdx % 2 == 1:false);
            return currentIdx++ < readsList->getReadsCount();
        }

        void rewind() override {
            currentIdx = 0;
        }

        DefaultReadsListEntry &peekReadEntry() override {
            return entry;
        }
    };

}

#endif //PGTOOLS_READSLISTITERATOREXTENDEDWRAPPER_H
