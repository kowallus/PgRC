#ifndef PGTOOLS_READSSETPERSISTENCE_H
#define PGTOOLS_READSSETPERSISTENCE_H

#include "../DefaultReadsSet.h"
#include "../PackedReadsSet.h"
#include "../iterator/ReadsSetIterator.h"

namespace PgSAReadsSet {

    class ReadsSetPersistence {
    private:

        class ManagedReadsSetIterator: public ReadsSourceIteratorTemplate<uint_read_len_max> {
        private:
            vector<ReadsSourceIteratorTemplate< uint_read_len_max>*> coreIterators;
            ReadsSourceIteratorTemplate< uint_read_len_max>* readsIterator = 0;

            ifstream* srcSource = 0;
            ifstream* pairSource = 0;
            ifstream* divSource = 0;
        public:
            ManagedReadsSetIterator(const string &srcFile, const string &pairFile = "", const string &divisionFile = "",
                    bool divisionComplement = false, bool revComplPairFile = false);

        public:
            bool moveNextVirtual();
            string getReadVirtual();
            string getQualityInfoVirtual();
            uint_read_len_max getReadLengthVirtual();
            void rewindVirtual();

            virtual ~ManagedReadsSetIterator();
        };

    public:
        static ReadsSourceIteratorTemplate<uint_read_len_max>* createManagedReadsIterator(const string &srcFile,
                                                                                          const string &pairFile = "",
                                                                                          const string &divisionFile = "",
                                                                                          bool divisionComplement = false,
                                                                                          bool revComplPairFile = false);

        static vector<uint_reads_cnt_max> getReadsOriginalIndexes(const string &divisionFile = "",
                                                                  bool divisionComplement = false,
                                                                  uint64_t readsCount = 0);


        static void writeOutputDivision(const vector<uint_reads_cnt_max> &orgIndexesMapping, const vector<uint32_t> &readsFilterResult,
                        const uint32_t readNotMatchedValue, string divisionFile, bool divisionComplement);
    };

}


#endif //PGTOOLS_READSSETPERSISTENCE_H
