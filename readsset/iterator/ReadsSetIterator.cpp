#include "ReadsSetIterator.h"

namespace PgSAReadsSet {

    using namespace PgSAHelpers;

    void VectorMapping::saveMapping(string mappingFile) {
        std::ofstream mappingDest(mappingFile, std::ios::out | std::ios::binary | std::ios::trunc);
        if (mappingDest.fail()) {
            fprintf(stderr, "cannot write to indexes mapping file %s\n", mappingFile.c_str());
            exit(EXIT_FAILURE);
        }
        writeReadMode(mappingDest, PgSAHelpers::plainTextWriteMode);
        writeValue(mappingDest, readsCount);
        for(uint_reads_cnt_max orgIdx: this->mappingWithGuard) {
            writeValue(mappingDest, orgIdx);
        }
        writeValue(mappingDest, readsCount);
        mappingDest.close();
    }

    VectorMapping *VectorMapping::loadMapping(string mappingFile) {
        ifstream* divSource = new ifstream(mappingFile, std::ios::in | std::ios::binary);
        if (divSource->fail()) {
            fprintf(stderr, "cannot open reads indexes mapping file %s\n", mappingFile.c_str());
            exit(EXIT_FAILURE);
        }
        bool plainTextReadMode = confirmTextReadMode(*divSource);
        uint_reads_cnt_max readsCount;
        readValue(*divSource, readsCount, plainTextReadMode);
        VectorMapping* mapping = new VectorMapping({}, readsCount);
        uint_reads_cnt_max orgIdx = 0;
        do {
            readValue(*divSource, orgIdx, plainTextReadMode);
            mapping->mappingWithGuard.push_back(orgIdx);
        } while (orgIdx != readsCount);
        delete(divSource);
        return mapping;
    }

    vector<uint_reads_cnt_max> &VectorMapping::getMappingVector() {
        return mappingWithGuard;
    }

    template<typename uint_read_len>
    ConcatenatedReadsSourceIterator<uint_read_len>::ConcatenatedReadsSourceIterator(std::istream* source) {
        this->source = source;
    }

    template<typename uint_read_len>
    ConcatenatedReadsSourceIterator<uint_read_len>::~ConcatenatedReadsSourceIterator() {
    }

    template<typename uint_read_len>
    string ConcatenatedReadsSourceIterator<uint_read_len>::getRead() {
        return line.substr(0, length);
    }

    template<typename uint_read_len>
    uint_read_len ConcatenatedReadsSourceIterator<uint_read_len>::getReadLength() {
        return length;
    }

    template<typename uint_read_len>
    bool ConcatenatedReadsSourceIterator<uint_read_len>::moveNext() {
        if (!std::getline(*source, line)) 
            return false;

        for (length = 0; length < line.length(); length++)
            if (!isalpha(line[length]))
                break;

        counter++;
        return true;
    }
    
    template<typename uint_read_len>
    void ConcatenatedReadsSourceIterator<uint_read_len>::rewind() {
        counter = -1;
        source->clear();
        source->seekg(0);
    }

    template<typename uint_read_len>
    IndexesMapping* ConcatenatedReadsSourceIterator<uint_read_len>::retainVisitedIndexesMapping() {
        return new DirectMapping(counter + 1);
    }

    template<typename uint_read_len>
    FASTAReadsSourceIterator<uint_read_len>::FASTAReadsSourceIterator(std::istream* source) {
        this->source = source;
    }

    template<typename uint_read_len>
    FASTAReadsSourceIterator<uint_read_len>::~FASTAReadsSourceIterator() {

    }

    template<typename uint_read_len>
    FASTAReadsSourceIterator<uint_read_len>::FASTAReadsSourceIterator(std::istream* source, std::istream* pairSource) {
        this->source = source;
        this->pairSource = pairSource;
    }

    template<typename uint_read_len>
    string FASTAReadsSourceIterator<uint_read_len>::getRead() {
        return line.substr(0, length);
    }

    template<typename uint_read_len>
    uint_read_len FASTAReadsSourceIterator<uint_read_len>::getReadLength() {
        return length;
    }

    template<typename uint_read_len>
    bool FASTAReadsSourceIterator<uint_read_len>::moveNext() {
        std::istream* src = source;
        if (pair && pairSource)
            src = pairSource;
        pair = !pair;
        do {
            if (!std::getline(*src, line))
                return false;
        } while (line.find('>') == 0);
            
        for (length = 0; length < line.length(); length++)
            if (!isalpha(line[length]))
                break;

        counter++;
        return true;
    }

    template<typename uint_read_len>
    void FASTAReadsSourceIterator<uint_read_len>::rewind() {
        counter = -1;
        source->clear();
        source->seekg(0);
        if (pairSource) { 
            pairSource->clear();
            pairSource->seekg(0);
        }
        pair = false;
    }

    template<typename uint_read_len>
    IndexesMapping* FASTAReadsSourceIterator<uint_read_len>::retainVisitedIndexesMapping() {
        return new DirectMapping(counter + 1);
    }

    template<typename uint_read_len>
    FASTQReadsSourceIterator<uint_read_len>::FASTQReadsSourceIterator(const string &srcFile, const string &pairFile) {
        ownStreams = true;
        source = new ifstream(srcFile, ios_base::in | ios_base::binary);
        if (source->fail()) {
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
    }

    template<typename uint_read_len>
    FASTQReadsSourceIterator<uint_read_len>::FASTQReadsSourceIterator(std::ifstream* source, std::ifstream* pairSource) {
        ownStreams = false;
        this->source = source;
        this->pairSource = pairSource;
    }

    template<typename uint_read_len>
    FASTQReadsSourceIterator<uint_read_len>::~FASTQReadsSourceIterator() {
        if (ownStreams) {
            source->close();
            if (pairSource)
                pairSource->close();
            delete(source);
            delete(pairSource);
        }
    }

    template<typename uint_read_len>
    string FASTQReadsSourceIterator<uint_read_len>::getRead() {
        return line.substr(0, length);
    }

    template<typename uint_read_len>
    string FASTQReadsSourceIterator<uint_read_len>::getQualityInfo() {
        return quality.substr(0, length);
    }

    template<typename uint_read_len>
    uint_read_len FASTQReadsSourceIterator<uint_read_len>::getReadLength() {
        return length;
    }

    template<typename uint_read_len>
    bool FASTQReadsSourceIterator<uint_read_len>::moveNext() {
        std::istream* src = source;
        if (pair && pairSource)
            src = pairSource;
        pair = !pair;

        if (!std::getline(*src, id))
            return false;
        std::getline(*src, line);
        std::getline(*src, opt_id);
        std::getline(*src, quality);

        for (length = 0; length < line.length(); length++)
            if (!isalpha(line[length]))
                break;

        counter++;
        return true;
    }

    template<typename uint_read_len>
    void FASTQReadsSourceIterator<uint_read_len>::rewind() {
        counter = -1;
        source->clear();
        source->seekg(0);
        if (pairSource) { 
            pairSource->clear();
            pairSource->seekg(0);
        }
        pair = false;
    }


    template<typename uint_read_len>
    IndexesMapping* FASTQReadsSourceIterator<uint_read_len>::retainVisitedIndexesMapping() {
        return new DirectMapping(counter + 1);
    }


    template<typename uint_read_len>
    ReadsSourceIteratorTemplate<uint_read_len>::~ReadsSourceIteratorTemplate() {
    }

    template<typename uint_read_len>
    RevComplPairReadsSetIterator<uint_read_len>::RevComplPairReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator): coreIterator(coreIterator) {
    }

    template<typename uint_read_len>
    bool RevComplPairReadsSetIterator<uint_read_len>::moveNext() {
        if (coreIterator->moveNext()) {
            counter++;
            return true;
        }
        return false;
    }

    template<typename uint_read_len>
    string RevComplPairReadsSetIterator<uint_read_len>::getRead() {
        string read = coreIterator->getRead();
        if (counter % 2)
            PgSAHelpers::reverseComplementInPlace(read);
        return read;
    }

    template<typename uint_read_len>
    string RevComplPairReadsSetIterator<uint_read_len>::getQualityInfo() {
        string q = coreIterator->getQualityInfo();
        if (counter % 2)
            std::reverse(q.begin(), q.end());
        return q;
    }

    template<typename uint_read_len>
    uint_read_len RevComplPairReadsSetIterator<uint_read_len>::getReadLength() {
        return coreIterator->getReadLength();
    }

    template<typename uint_read_len>
    void RevComplPairReadsSetIterator<uint_read_len>::rewind() {
        counter = -1;
        coreIterator->rewind();
    }

    template<typename uint_read_len>
    IndexesMapping *RevComplPairReadsSetIterator<uint_read_len>::retainVisitedIndexesMapping() {
        return coreIterator->retainVisitedIndexesMapping();
    }


    template<typename uint_read_len>
    IgnoreNReadsSetIterator<uint_read_len>::IgnoreNReadsSetIterator(
            ReadsSourceIteratorTemplate<uint_read_len> *coreIterator): coreIterator(coreIterator) {
    }

    template<typename uint_read_len>
    bool IgnoreNReadsSetIterator<uint_read_len>::moveNext() {
        while (coreIterator->moveNext()) {
            if (isFreeOfN()) {
                indexesMapping.push_back(++counter);
                return true;
            }
        }
        return false;
    }

    template<typename uint_read_len>
    string IgnoreNReadsSetIterator<uint_read_len>::getRead() {
        return coreIterator->getRead();
    }

    template<typename uint_read_len>
    string IgnoreNReadsSetIterator<uint_read_len>::getQualityInfo() {
        return coreIterator->getQualityInfo();
    }

    template<typename uint_read_len>
    uint_read_len IgnoreNReadsSetIterator<uint_read_len>::getReadLength() {
        return coreIterator->getReadLength();
    }

    template<typename uint_read_len>
    void IgnoreNReadsSetIterator<uint_read_len>::rewind() {
        counter = -1;
        indexesMapping.clear();
        coreIterator->rewind();
    }

    template<typename uint_read_len>
    bool IgnoreNReadsSetIterator<uint_read_len>::isFreeOfN() {
        return coreIterator->getRead().find('N') == string::npos;
    }

    template<typename uint_read_len>
    IndexesMapping* IgnoreNReadsSetIterator<uint_read_len>::retainVisitedIndexesMapping() {
        IndexesMapping *const coreMapping = coreIterator->retainVisitedIndexesMapping();
        uint_reads_cnt_max readsTotalCount = coreMapping->getReadsTotalCount();
        delete(coreMapping);
        return new VectorMapping(std::move(indexesMapping), readsTotalCount);
    }


    template class ReadsSourceIteratorTemplate<uint_read_len_min>;
    template class ReadsSourceIteratorTemplate<uint_read_len_std>;
    template class ConcatenatedReadsSourceIterator<uint_read_len_min>;
    template class ConcatenatedReadsSourceIterator<uint_read_len_std>;
    template class FASTAReadsSourceIterator<uint_read_len_min>;
    template class FASTAReadsSourceIterator<uint_read_len_std>;
    template class FASTQReadsSourceIterator<uint_read_len_min>;
    template class FASTQReadsSourceIterator<uint_read_len_std>;
    template class RevComplPairReadsSetIterator<uint_read_len_min>;
    template class RevComplPairReadsSetIterator<uint_read_len_std>;
    template class IgnoreNReadsSetIterator<uint_read_len_min>;
    template class IgnoreNReadsSetIterator<uint_read_len_std>;
}
