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
                    bool divisionComplement = false, bool revComplPairFile = false,
                    bool ignoreNReads = false, bool ignoreNoNReads = false);

        public:
            bool moveNextVirtual();
            string getReadVirtual();
            string getQualityInfoVirtual();
            uint_read_len_max getReadLengthVirtual();
            void rewindVirtual();

            virtual ~ManagedReadsSetIterator();

            const vector<uint_reads_cnt_max> getVisitedIndexesMapping() override;
        };

    public:
        static ReadsSourceIteratorTemplate<uint_read_len_max>* createManagedReadsIterator(const string &srcFile,
                                                                                          const string &pairFile = "",
                                                                                          const string &divisionFile = "",
                                                                                          bool divisionComplement = false,
                                                                                          bool revComplPairFile = false,
                                                                                          bool ignoreNReads = false,
                                                                                          bool ignoreNoNReads = false);

        template<typename filter_res_t>
        static void writeOutputDivision(const vector<uint_reads_cnt_max> &orgIndexesMapping, const vector<filter_res_t> &readsFilterResult,
                        const filter_res_t readSelectedValue, string divisionFile, bool divisionComplement) {
            std::ofstream divDest(divisionFile, std::ios::out | std::ios::binary);
            if (divDest.fail()) {
                fprintf(stderr, "cannot write to division file %s\n", divisionFile.c_str());
                exit(EXIT_FAILURE);
            }
            divDest << (plainTextWriteMode?TEXT_MODE_ID:BINARY_MODE_ID) << endl;

            int64_t i = -1;
            uint64_t visitedReadsCount = orgIndexesMapping.size() - 1;

            if (divisionComplement) {
                int64_t counter = -1;
                while (++i < visitedReadsCount) {
                    while (++counter != orgIndexesMapping[i])
                        writeValue(divDest, counter);
                    if (readsFilterResult[i] != readSelectedValue)
                        writeValue(divDest, orgIndexesMapping[i]);
                }
                while (++counter != orgIndexesMapping[i])
                    writeValue(divDest, counter);
            } else {
                for(uint64_t i = 0; i < visitedReadsCount; i++)
                    if (readsFilterResult[i] == readSelectedValue)
                        writeValue(divDest, orgIndexesMapping[i]);
            }

            writeValue(divDest, UINT64_MAX);
            divDest.close();
        }
    };

}


#endif //PGTOOLS_READSSETPERSISTENCE_H
