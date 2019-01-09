#include "ReadsSetPersistence.h"

#include "../iterator/DivisionReadsSetDecorators.h"

namespace PgSAReadsSet {

    ReadsSourceIteratorTemplate<uint_read_len_max> *ReadsSetPersistence::createManagedReadsIterator(const string &srcFile,
                                                                                                    const string &pairFile,
                                                                                                    const string &divisionFile,
                                                                                                    bool divisionComplement,
                                                                                                    bool revComplPairFile,
                                                                                                    bool ignoreNReads,
                                                                                                    bool ignoreNoNReads) {
        return new ManagedReadsSetIterator(srcFile, pairFile, divisionFile, divisionComplement, revComplPairFile,
                                           ignoreNReads, ignoreNoNReads);
    }

    using namespace PgTools;

    ReadsSetPersistence::ManagedReadsSetIterator::ManagedReadsSetIterator(const string &srcFile, const string &pairFile,
            const string &divisionFile, bool divisionComplement, bool revComplPairFile, bool ignoreNReads, bool ignoreNoNReads) {
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
            readsIterator = new DividedReadsSetIterator<uint_read_len_max>(readsIterator, divSource, divisionComplement,
                    ignoreNReads, ignoreNoNReads);
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

    IndexesMapping* ReadsSetPersistence::ManagedReadsSetIterator::retainVisitedIndexesMapping() {
        return readsIterator->retainVisitedIndexesMapping();
    }
}