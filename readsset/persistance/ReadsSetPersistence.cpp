#include "ReadsSetPersistence.h"

#include "../iterator/DivisionReadsSetDecorators.h"

namespace PgSAReadsSet {

    vector<uint_reads_cnt_max>
    ReadsSetPersistence::getReadsOriginalIndexes(const string &divisionFile, bool divisionComplement,
                                                 uint64_t readsCount) {
        vector<uint_reads_cnt_max> mapping(readsCount);
        if (divisionFile == "") {
            for(uint64_t i = 0; i < readsCount; i++)
                mapping[i] = i;
        } else {
            std::ifstream divSource(divisionFile, std::ios::in | std::ios::binary);
            if (divSource.fail()) {
                fprintf(stderr, "cannot open reads division file %s\n", divisionFile.c_str());
                exit(EXIT_FAILURE);
            }
            uint64_t i = 0;
            uint64_t counter = 0;
            uint64_t currentDivIdx = 0;
            bool plainTextReadMode = readReadMode(divSource);
            readValue(divSource, currentDivIdx, plainTextReadMode);
            while (i < readsCount) {
                if (divisionComplement) {
                    while (counter == currentDivIdx) {
                        readValue(divSource, currentDivIdx, plainTextReadMode);
                        counter++;
                    }
                    mapping[i++] = counter++;
                } else {
                    while (counter != currentDivIdx)
                        counter++;
                    mapping[i++] = counter++;
                    readValue(divSource, currentDivIdx, plainTextReadMode);
                }
            }
        }
        return mapping;
    }

    void ReadsSetPersistence::writeOutputDivision(const vector<uint_reads_cnt_max> &orgIndexesMapping,
                                                  const vector<uint32_t> &readsFilterResult,
                                                  const uint32_t readNotMatchedValue, string divisionFile,
                                                  bool divisionComplement) {
        std::ofstream divDest(divisionFile, std::ios::out | std::ios::binary);
        if (divDest.fail()) {
            fprintf(stderr, "cannot write to division file %s\n", divisionFile.c_str());
            exit(EXIT_FAILURE);
        }
        divDest << (plainTextWriteMode?TEXT_MODE_ID:BINARY_MODE_ID) << endl;

        int64_t i = -1;
        uint64_t readsCount = orgIndexesMapping.size();

        if (divisionComplement) {
            uint64_t counter = 0;
            while (++i < readsCount) {
                while (counter != orgIndexesMapping[i]) {
                    writeValue(divDest, counter);
                    counter++;
                }
                if (readsFilterResult[i] != readNotMatchedValue)
                    writeValue(divDest, orgIndexesMapping[i]);
            }
        } else {
            for(uint64_t i = 0; i < readsCount; i++)
                if (readsFilterResult[i] == readNotMatchedValue)
                    writeValue(divDest, orgIndexesMapping[i]);
        }

        writeValue(divDest, UINT64_MAX);
        divDest.close();
    }

    ReadsSourceIteratorTemplate<uint_read_len_max> *ReadsSetPersistence::createManagedReadsIterator(const string &srcFile,
                                                                                                    const string &pairFile,
                                                                                                    bool ignoreNReads,
                                                                                                    const string &divisionFile,
                                                                                                    bool divisionComplement,
                                                                                                    bool revComplPairFile) {
        return new ManagedReadsSetIterator(srcFile, pairFile, ignoreNReads, divisionFile, divisionComplement, revComplPairFile);
    }

    using namespace PgTools;

    ReadsSetPersistence::ManagedReadsSetIterator::ManagedReadsSetIterator(const string &srcFile, const string &pairFile,
            bool ignoreNReads, const string &divisionFile, bool divisionComplement, bool revComplPairFile) {
        srcSource = new ifstream(srcFile, ios_base::in | ios_base::binary);
        if (srcSource->fail()) {
            fprintf(stderr, "cannot open reads file %s\n", srcFile.c_str());
            exit(EXIT_FAILURE);
        }
        if (pairFile != "") {
            pairSource = new ifstream(pairFile, ios_base::in | ios_base::binary);
            if (pairSource->fail()) {
                fprintf(stderr, "cannot open reads pair file %s\n", pairFile.c_str());
                exit(EXIT_FAILURE);
            }
        }

        if (srcFile.substr(srcFile.length() - 6) == ".fasta")
            readsIterator = new FASTAReadsSourceIterator<uint_read_len_max>(srcSource, pairSource);
        else if (srcFile.substr(srcFile.length() - 6) == ".fastq")
            readsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(srcSource, pairSource);
        else
            readsIterator = new ConcatenatedReadsSourceIterator<uint_read_len_max>(srcSource);

        if ((pairFile != "") && revComplPairFile) {
            coreIterators.push_back(readsIterator);
            readsIterator = new RevComplPairReadsSetIterator<uint_read_len_max>(readsIterator);
        }

        if (divisionFile != "") {
            coreIterators.push_back(readsIterator);
            divSource = new ifstream(divisionFile, std::ios::in | std::ios::binary);
            if (divSource->fail()) {
                fprintf(stderr, "cannot open reads division file %s\n", divisionFile.c_str());
                exit(EXIT_FAILURE);
            }
            readsIterator = new DividedReadsSetIterator<uint_read_len_max>(readsIterator, divSource, divisionComplement);
        }

        if (ignoreNReads) {
            coreIterators.push_back(readsIterator);
            readsIterator = new IgnoreNReadsSetIterator<uint_read_len_max>(readsIterator);
        }
    }

    ReadsSetPersistence::ManagedReadsSetIterator::~ManagedReadsSetIterator() {
        delete(readsIterator);
        for(ReadsSourceIteratorTemplate<uint_read_len_max>* coreIterator: coreIterators)
            delete(coreIterator);
        srcSource->close();
        if (pairSource)
            pairSource->close();
        if (divSource)
            divSource->close();
        delete(srcSource);
        delete(pairSource);
        delete(divSource);
    }

    bool ReadsSetPersistence::ManagedReadsSetIterator::moveNextVirtual() {
        return readsIterator->moveNextVirtual();
    }

    string ReadsSetPersistence::ManagedReadsSetIterator::getReadVirtual() {
        return readsIterator->getReadVirtual();
    }

    string ReadsSetPersistence::ManagedReadsSetIterator::getQualityInfoVirtual() {
        return readsIterator->getQualityInfoVirtual();
    }

    uint_read_len_max ReadsSetPersistence::ManagedReadsSetIterator::getReadLengthVirtual() {
        return readsIterator->getReadLengthVirtual();
    }

    void ReadsSetPersistence::ManagedReadsSetIterator::rewindVirtual() {
        readsIterator->rewindVirtual();
    }

}