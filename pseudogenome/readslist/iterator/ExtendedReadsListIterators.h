#ifndef PGTOOLS_EXTENDEDREADSLISTITERATORS_H
#define PGTOOLS_EXTENDEDREADSLISTITERATORS_H

#include "../ReadsListInterface.h"

#include "ExtendedReadsListIteratorInterface.h"
#include "../../../readsset/persistance/ReadsSetPersistence.h"

namespace PgTools {

    using namespace PgSAIndex;

    class ReadsListIteratorExtendedWrapperBase: public DefaultReadsListIteratorInterface {
    public:
        virtual void applyDivision(string divisionFile, bool divisionComplement) = 0;
    };

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    class ReadsListIteratorExtendedWrapper: public ReadsListIteratorExtendedWrapperBase {
    private:
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList;
        uint_reads_cnt currentIdx = 0;
        DefaultReadsListEntry entry;
        bool mapping = false;
        vector<uint_reads_cnt_max> orgIndexesMapping;

    public:
        ReadsListIteratorExtendedWrapper(ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass> *readsList)
        : readsList(readsList) {}

        ~ReadsListIteratorExtendedWrapper() { }

        inline uint_reads_cnt_max mapIndex(uint_reads_cnt idx) {
            return mapping?orgIndexesMapping[idx]:idx;
        }

        void applyDivision(string divisionFile, bool divisionComplement) {
            mapping = true;
            orgIndexesMapping = ReadsSetPersistence::getReadsOriginalIndexes(divisionFile, divisionComplement, readsList->getReadsCount());
        }

        bool moveNext() override {
            entry.advanceEntryByPosition(readsList->getReadPosition(currentIdx), mapIndex(readsList->getReadOriginalIndex(currentIdx)));
            return currentIdx++ < readsList->getReadsCount();
        }

        const ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &peekReadEntry() override {
            return entry;
        }
    };

    class SeparatedExtendedReadsListIterator: public DefaultReadsListIteratorInterface {
    private:
        const string &pseudoGenomePrefix;
        ifstream* rlPosSrc = 0;
        ifstream* rlOrgIdxSrc = 0;
        ifstream* rlRevCompSrc = 0;
        ifstream* rlMisCntSrc = 0;
        ifstream* rlMisSymSrc = 0;
        ifstream* rlMisOffSrc = 0;

        DefaultReadsListEntry entry;

        void initSrc(ifstream* &src);
        void initSrcs();
        void freeSrc(ifstream* &src);
        void freeSrcs();
    public:
        SeparatedExtendedReadsListIterator(const string &pseudoGenomePrefix);

        ~SeparatedExtendedReadsListIterator() override;

        bool moveNext() override;

        const ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &peekReadEntry() override;
    };
}

#endif //PGTOOLS_EXTENDEDREADSLISTITERATORS_H
