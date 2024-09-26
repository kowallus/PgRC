#include "ReadsSetPersistence.h"

#include "../iterator/DivisionReadsSetDecorators.h"

namespace PgReadsSet {

    ReadsSourceIteratorTemplate<uint_read_len_max> *ReadsSetPersistence::createManagedReadsIterator(const string &srcFile,
                                                                                                    const string &pairFile,
                                                                                                    bool revComplPairFile,
                                                                                                    const string &divisionFile,
                                                                                                    bool divisionComplement,
                                                                                                    bool ignoreNReads,
                                                                                                    bool ignoreNoNReads) {
        return new ManagedReadsSetIterator(srcFile, pairFile, revComplPairFile, divisionFile, divisionComplement,
                                           ignoreNReads, ignoreNoNReads);
    }

    using namespace PgTools;

    ReadsSetPersistence::ManagedReadsSetIterator::ManagedReadsSetIterator(const string &srcFile, const string &pairFile,
            bool revComplPairFile, const string &divisionFile, bool divisionComplement, bool ignoreNReads, bool ignoreNoNReads) {
        srcSource = new ifstream(srcFile, ios_base::in | ios_base::binary);
        if (srcSource->fail()) {
            fprintf(stderr, "cannot open reads file %s\n", srcFile.c_str());
            exit(EXIT_FAILURE);
        }
        srcSource->rdbuf()->pubsetbuf(buf1, 1 << 16);
        if (pairFile != "") {
            pairSource = new ifstream(pairFile, ios_base::in | ios_base::binary);
            if (pairSource->fail()) {
                fprintf(stderr, "cannot open reads pair file %s\n", pairFile.c_str());
                exit(EXIT_FAILURE);
            }
            pairSource->rdbuf()->pubsetbuf(buf2, 1 << 16);
        }
        char firstSymbol = srcSource->get();
        srcSource->clear();
        srcSource->seekg(0);
        switch(firstSymbol) {case '@': readsIterator = new FASTQReadsSourceIterator<uint_read_len_max>(srcSource, pairSource);
                break;
            case ';':
            case '>':
                readsIterator = new FASTAReadsSourceIterator<uint_read_len_max>(srcSource, pairSource);
                break;
            default:
                readsIterator = new ConcatenatedReadsSourceIterator<uint_read_len_max>(srcSource);
        }
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

    bool ReadsSetPersistence::ManagedReadsSetIterator::moveNext() {
        return readsIterator->moveNext();
    }

    string& ReadsSetPersistence::ManagedReadsSetIterator::getRead() {
        return readsIterator->getRead();
    }

    string& ReadsSetPersistence::ManagedReadsSetIterator::getQualityInfo() {
        return readsIterator->getQualityInfo();
    }

    uint_read_len_max ReadsSetPersistence::ManagedReadsSetIterator::getReadLength() {
        return readsIterator->getReadLength();
    }

    void ReadsSetPersistence::ManagedReadsSetIterator::rewind() {
        readsIterator->rewind();
    }

    IndexesMapping* ReadsSetPersistence::ManagedReadsSetIterator::retainVisitedIndexesMapping() {
        return readsIterator->retainVisitedIndexesMapping();
    }
}