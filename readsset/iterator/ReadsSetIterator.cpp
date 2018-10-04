#include "ReadsSetIterator.h"

namespace PgSAReadsSet {

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
    uint_read_len ConcatenatedReadsSourceIterator<uint_read_len>::getReadLengthVirtual() {
        return getReadLength();
    }

    template<typename uint_read_len>
    string ConcatenatedReadsSourceIterator<uint_read_len>::getReadVirtual() {
        return getRead();
    }

    template<typename uint_read_len>
    bool ConcatenatedReadsSourceIterator<uint_read_len>::moveNext() {
        if (!std::getline(*source, line)) 
            return false;

        for (length = 0; length < line.length(); length++)
            if (!isalpha(line[length]))
                break;

        return true;
    }
    
    template<typename uint_read_len>
    void ConcatenatedReadsSourceIterator<uint_read_len>::rewind() {
        source->clear();
        source->seekg(0);
    }

    template<typename uint_read_len>
    bool ConcatenatedReadsSourceIterator<uint_read_len>::moveNextVirtual() {
        return moveNext();
    }

    template<typename uint_read_len>
    void ConcatenatedReadsSourceIterator<uint_read_len>::rewindVirtual() {
        rewind();
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
    uint_read_len FASTAReadsSourceIterator<uint_read_len>::getReadLengthVirtual() {
        return getReadLength();
    }

    template<typename uint_read_len>
    string FASTAReadsSourceIterator<uint_read_len>::getReadVirtual() {
        return getRead();
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

        return true;
    }

    template<typename uint_read_len>
    void FASTAReadsSourceIterator<uint_read_len>::rewind() {
        source->clear();
        source->seekg(0);
        if (pairSource) { 
            pairSource->clear();
            pairSource->seekg(0);
        }
    }

    template<typename uint_read_len>
    bool FASTAReadsSourceIterator<uint_read_len>::moveNextVirtual() {
        return moveNext();
    }

    template<typename uint_read_len>
    void FASTAReadsSourceIterator<uint_read_len>::rewindVirtual() {
        rewind();
    }

    template<typename uint_read_len>
    FASTQReadsSourceIterator<uint_read_len>::FASTQReadsSourceIterator(std::istream* source, std::istream* pairSource) {
        this->source = source;
        this->pairSource = pairSource;
    }

    template<typename uint_read_len>
    FASTQReadsSourceIterator<uint_read_len>::~FASTQReadsSourceIterator() {
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
    uint_read_len FASTQReadsSourceIterator<uint_read_len>::getReadLengthVirtual() {
        return getReadLength();
    }

    template<typename uint_read_len>
    string FASTQReadsSourceIterator<uint_read_len>::getReadVirtual() {
        return getRead();
    }

    template<typename uint_read_len>
    string FASTQReadsSourceIterator<uint_read_len>::getQualityInfoVirtual() {
        return getQualityInfo();
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

        return true;
    }

    template<typename uint_read_len>
    void FASTQReadsSourceIterator<uint_read_len>::rewind() {
        source->clear();
        source->seekg(0);
        if (pairSource) { 
            pairSource->clear();
            pairSource->seekg(0);
        }
    }
    
    template<typename uint_read_len>
    bool FASTQReadsSourceIterator<uint_read_len>::moveNextVirtual() {
        return moveNext();
    }

    template<typename uint_read_len>
    void FASTQReadsSourceIterator<uint_read_len>::rewindVirtual() {
        rewind();
    }

    template<typename uint_read_len>
    ReadsSourceIteratorTemplate<uint_read_len>::~ReadsSourceIteratorTemplate() {
    }

    template class ReadsSourceIteratorTemplate<uint_read_len_min>;
    template class ReadsSourceIteratorTemplate<uint_read_len_std>;
    template class ConcatenatedReadsSourceIterator<uint_read_len_min>;
    template class ConcatenatedReadsSourceIterator<uint_read_len_std>;
    template class FASTAReadsSourceIterator<uint_read_len_min>;
    template class FASTAReadsSourceIterator<uint_read_len_std>;
    template class FASTQReadsSourceIterator<uint_read_len_min>;
    template class FASTQReadsSourceIterator<uint_read_len_std>;
    
}
