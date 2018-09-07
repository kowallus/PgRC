#ifndef DEFAULTREADSSET_H_INCLUDED
#define DEFAULTREADSSET_H_INCLUDED

#include <algorithm>
#include "ReadsSetBase.h"
#include "iterator/ReadsSetIterator.h"
#include "ReadsSetInterface.h"

namespace PgSAReadsSet {

    class DefaultReadsSet: public ReadsSetBase, public ReadsSetInterface< uint_read_len_max, uint_reads_cnt_max >
    {
        private:

            vector<string> reads;

        public:

            template<class ReadsSourceIterator>
            DefaultReadsSet(ReadsSourceIterator*);

            virtual ~DefaultReadsSet() {};

            inline uint_read_len_max maxReadLength() { return properties->maxReadLength; };
            inline uint_reads_cnt_max readsCount() { return properties->readsCount; };

            inline bool isReadLengthConstant() { return properties->constantReadLength; };

            inline const string getRead(uint_reads_cnt_max i) { return reads[i];};
            inline uint_read_len_max readLength(uint_reads_cnt_max i) { return reads[i].length(); };

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
           
            static DefaultReadsSet* readReadsSet(std::string filename, std::string pairfile = "") {
                istream* streamSource = new std::ifstream(filename, std::ios::in | std::ios::binary);
                istream* pairSource = 0;
                if (pairfile != "")
                    pairSource = new std::ifstream(pairfile, std::ios::in | std::ios::binary);
                        
                DefaultReadsSet* readsSet;
                if (filename.substr(filename.length() - 6) == ".fasta") {
                    FASTAReadsSourceIterator<uint_read_len_max>* readsSource = new FASTAReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
                    readsSet = new DefaultReadsSet(readsSource);
                } else if (filename.substr(filename.length() - 6) == ".fastq") {
                    FASTQReadsSourceIterator<uint_read_len_max>* readsSource = new FASTQReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
                    readsSet = new DefaultReadsSet(readsSource);
                } else {
                    ConcatenatedReadsSourceIterator<uint_read_len_max>* readsSource = new ConcatenatedReadsSourceIterator<uint_read_len_max>(streamSource);
                    readsSet = new DefaultReadsSet(readsSource);
                }
                delete(streamSource);
                return readsSet;
            }
    };

    //typedef DefaultReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_max>> ConcatenatedReadsSet;
}

#endif // DEFAULTREADSSET_H_INCLUDED
