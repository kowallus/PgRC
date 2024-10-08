#ifndef DEFAULTREADSSET_H_INCLUDED
#define DEFAULTREADSSET_H_INCLUDED

#include <algorithm>
#include "ReadsSetBase.h"
#include "iterator/ReadsSetIterator.h"
#include "ReadsSetInterface.h"

namespace PgReadsSet {

    class DefaultReadsSet: public ReadsSetBase, public ReadsSetInterface< uint_read_len_max, uint_reads_cnt_max >
    {
        private:

            vector<string> reads;

        public:

            template<class ReadsSourceIterator>
            DefaultReadsSet(ReadsSourceIterator*);

            ~DefaultReadsSet() override {};

            inline uint_read_len_max maxReadLength() override { return properties->maxReadLength; };
            inline uint_reads_cnt_max readsCount() override { return properties->readsCount; };

            inline bool isReadLengthConstant() override { return properties->constantReadLength; };

            inline const string getRead(uint_reads_cnt_max i) override { return reads[i];};
            inline uint_read_len_max readLength(uint_reads_cnt_max i) override { return reads[i].length(); };

            uint_read_len_max maxReadLengthVirtual() { return maxReadLength(); };
            uint_reads_cnt_max readsCountVirtual() { return readsCount(); };
            bool isReadLengthConstantVirtual() { return isReadLengthConstant(); };
            const string getReadVirtual(uint_reads_cnt_max i) { return getRead(i); };
            uint_read_len_max readLengthVirtual(uint_reads_cnt_max i) { return readLength(i); };

            void printout() {
                properties->printout();

//                for (auto readIt = reads.begin(); readIt != reads.end(); readIt++)
//                    cout << *readIt << "\n";
            }

    };

    //typedef DefaultReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_max>> ConcatenatedReadsSet;
}

#endif // DEFAULTREADSSET_H_INCLUDED
