#ifndef PACKEDREADSSET_H_INCLUDED
#define PACKEDREADSSET_H_INCLUDED

#include <algorithm>
#include "ReadsSetBase.h"
#include "iterator/ReadsSetIterator.h"
#include "ReadsSetInterface.h"
#include "../pseudogenome/packing/SymbolsPackingFacility.h"

using namespace PgSAIndex;

namespace PgSAReadsSet {

    class PackedReadsSet: public ReadsSetBase, public ReadsSetInterface< uint_read_len_max, uint_reads_cnt_max >
    {
        private:

            uint_ps_element_min* packedReads;
            uchar packedLength;
            vector<uint_read_len_max> lengths;

        public:

            SymbolsPackingFacility<uint_ps_element_min>* sPacker;
            
            template<class ReadsSourceIterator>
            PackedReadsSet(ReadsSourceIterator*);

            virtual ~PackedReadsSet();

            inline uint_read_len_max maxReadLength() { return properties->maxReadLength; };
            inline uint_reads_cnt_max readsCount() { return properties->readsCount; };

            inline bool isReadLengthConstant() { return properties->constantReadLength; };

            inline const uint_ps_element_min* getPackedRead(uint_reads_cnt_max i) { return packedReads + i * (size_t) packedLength;};
            inline const string getReadPrefix(uint_reads_cnt_max i, uint_read_len_max skipSuffix) { return sPacker->reverseSequence(packedReads + i * (size_t) packedLength, 0, readLength(i) - skipSuffix);};
            inline const string getRead(uint_reads_cnt_max i) { return sPacker->reverseSequence(packedReads + (size_t) packedLength * i, 0, readLength(i));};
            inline uint_read_len_max readLength(uint_reads_cnt_max i) { return lengths[i]; };

            int comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx);
            int comparePackedReads(uint_reads_cnt_max lIdx, uint_reads_cnt_max rIdx, uint_read_len_max offset);
            int compareSuffixWithPrefix(uint_reads_cnt_max sufIdx, uint_reads_cnt_max preIdx, uint_read_len_max sufOffset);         
            
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
           
            static PackedReadsSet* readReadsSet(std::string filename, std::string pairfile = "") {
                istream* streamSource = new std::ifstream(filename, std::ios::in | std::ios::binary);
                istream* pairSource = 0;
                if (pairfile != "")
                    pairSource = new std::ifstream(pairfile, std::ios::in | std::ios::binary);
                        
                PackedReadsSet* readsSet;
                if (filename.substr(filename.length() - 6) == ".fasta") {
                    FASTAReadsSourceIterator<uint_read_len_max>* readsSource = new FASTAReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
                    readsSet = new PackedReadsSet(readsSource);
                    delete(readsSource);
                } else if (filename.substr(filename.length() - 6) == ".fastq") {
                    FASTQReadsSourceIterator<uint_read_len_max>* readsSource = new FASTQReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
                    readsSet = new PackedReadsSet(readsSource);
                    delete(readsSource);
                } else {
                    ConcatenatedReadsSourceIterator<uint_read_len_max>* readsSource = new ConcatenatedReadsSourceIterator<uint_read_len_max>(streamSource);
                    readsSet = new PackedReadsSet(readsSource);
                    delete(readsSource);
                }
                delete(streamSource);
                if (pairSource)
                    delete(pairSource);
                return readsSet;
            }
    };

    //typedef PackedReadsSet<ConcatenatedReadsSourceIterator<uint_read_len_max>> ConcatenatedReadsSet;
}

#endif // PACKEDREADSSET_H_INCLUDED
