#ifndef PGTOOLS_READSSETPERSISTENCE_H
#define PGTOOLS_READSSETPERSISTENCE_H

#include "../DefaultReadsSet.h"
#include "../PackedConstantLengthReadsSet.h"
#include "../iterator/ReadsSetIterator.h"

namespace PgReadsSet {

    class ReadsSetPersistence {
    private:

        class ManagedReadsSetIterator: public ReadsSourceIteratorTemplate<uint_read_len_max> {
        private:
            vector<ReadsSourceIteratorTemplate< uint_read_len_max>*> coreIterators;
            ReadsSourceIteratorTemplate< uint_read_len_max>* readsIterator = nullptr;

            char buf1[1 << 16];
            char buf2[1 << 16];

            ifstream* srcSource = nullptr;
            ifstream* pairSource = nullptr;
            ifstream* divSource = nullptr;
        public:
            ManagedReadsSetIterator(const string &srcFile, const string &pairFile = "", bool revComplPairFile = false,
                    const string &divisionFile = "", bool divisionComplement = false,
                    bool ignoreNReads = false, bool ignoreNoNReads = false);

        public:
            bool moveNext() override;
            string& getRead() override;
            string& getQualityInfo() override;
            uint_read_len_max getReadLength() override;
            void rewind() override;

            ~ManagedReadsSetIterator() override;

            IndexesMapping* retainVisitedIndexesMapping() override;
        };

    public:
        static ReadsSourceIteratorTemplate<uint_read_len_max>* createManagedReadsIterator(const string &srcFile,
                                                                                          const string &pairFile = "",
                                                                                          bool revComplPairFile = false,
                                                                                          const string &divisionFile = "",
                                                                                          bool divisionComplement = false,
                                                                                          bool ignoreNReads = false,
                                                                                          bool ignoreNoNReads = false);


        template<typename filter_res_t>
        static void writeOutputDivision(IndexesMapping* orgIndexesMapping, const vector<filter_res_t> &readsFilterResult,
                        const filter_res_t readSelectedValue, string divisionFile, bool divisionComplement) {
            std::ofstream divDest(divisionFile, std::ios::out | std::ios::binary);
            if (divDest.fail()) {
                fprintf(stderr, "cannot write to division file %s\n", divisionFile.c_str());
                exit(EXIT_FAILURE);
            }
            divDest << (plainTextWriteMode?TEXT_MODE_ID:BINARY_MODE_ID) << endl;
            writeValue(divDest, orgIndexesMapping->getReadsTotalCount());

            int64_t i = -1;
            uint_reads_cnt_max visitedReadsCount = orgIndexesMapping->getMappedReadsCount();

            if (divisionComplement) {
                int64_t counter = -1;
                while (++i < visitedReadsCount) {
                    while (++counter != orgIndexesMapping->getReadOriginalIndex(i))
                        writeValue(divDest, counter);
                    if (readsFilterResult[i] != readSelectedValue)
                        writeValue(divDest, orgIndexesMapping->getReadOriginalIndex(i));
                }
                while (++counter != orgIndexesMapping->getReadsTotalCount())
                    writeValue(divDest, counter);
            } else {
                for(uint64_t i = 0; i < visitedReadsCount; i++)
                    if (readsFilterResult[i] == readSelectedValue)
                        writeValue(divDest, orgIndexesMapping->getReadOriginalIndex(i));
            }

            writeValue(divDest, UINT64_MAX);
            divDest.close();
        }
    };

}


#endif //PGTOOLS_READSSETPERSISTENCE_H
