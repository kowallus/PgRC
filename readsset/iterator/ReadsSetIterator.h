#ifndef READSSETITERATOR_H_INCLUDED
#define READSSETITERATOR_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include <vector>
#include "../../utils/helper.h"
#include "../../pgsaconfig.h"

namespace PgSAReadsSet {

    typedef uint_read_len_std uint_read_len_max;

    class IndexesMapping {
    public:
        virtual uint_reads_cnt_max getReadOriginalIndex(uint_reads_cnt_max idx) = 0;
        virtual uint_reads_cnt_max getMappedReadsCount() = 0;
        virtual uint_reads_cnt_max getReadsTotalCount() = 0;
    };

    class DirectMapping : public IndexesMapping {
    private:
        uint_reads_cnt_max readsCount;
    public:
        DirectMapping(uint_reads_cnt_max readsCount) : readsCount(readsCount) {}

        uint_reads_cnt_max getReadOriginalIndex(uint_reads_cnt_max idx) override { return idx; }
        uint_reads_cnt_max getMappedReadsCount() override { return readsCount; }
        uint_reads_cnt_max getReadsTotalCount() override { return readsCount; }
    };

    class VectorMapping : public IndexesMapping {
    private:
        std::vector<uint_reads_cnt_max> mapping;
        uint_reads_cnt_max readsCount;
    public:
        VectorMapping(vector<uint_reads_cnt_max> &&mapping, uint_reads_cnt_max readsCount) :
        mapping(std::move(mapping)), readsCount(readsCount) {}

        uint_reads_cnt_max getReadOriginalIndex(uint_reads_cnt_max idx) override { return mapping[idx]; }
        uint_reads_cnt_max getMappedReadsCount() override { return mapping.size(); }
        uint_reads_cnt_max getReadsTotalCount() override { return readsCount; }
    };

    template < typename uint_read_len >
    class ReadsSourceIteratorTemplate
    {
        public:

            virtual ~ReadsSourceIteratorTemplate();

            virtual bool moveNext() = 0;
            virtual string getRead() = 0;
            virtual string getQualityInfo() { return std::string(); };
            virtual uint_read_len getReadLength() = 0;
            virtual void rewind() = 0;

            virtual IndexesMapping* retainVisitedIndexesMapping() = 0;
    };

    template < typename uint_read_len >
    class ConcatenatedReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = 0;
            int64_t counter = -1;

        public:

            ConcatenatedReadsSourceIterator(std::istream* source);

            ~ConcatenatedReadsSourceIterator();

            bool moveNext();
            string getRead();
            uint_read_len getReadLength();
            void rewind();

            IndexesMapping* retainVisitedIndexesMapping() override;
    };
    
    template < typename uint_read_len >
    class FASTAReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = 0;
            std::istream* pairSource = 0;
            bool pair = false;
            int64_t counter = -1;
            
        public:

            FASTAReadsSourceIterator(std::istream* source);
            
            FASTAReadsSourceIterator(std::istream* source, std::istream* pairSource);

            ~FASTAReadsSourceIterator();

            bool moveNext();
            string getRead();
            uint_read_len getReadLength();
            void rewind();

            IndexesMapping* retainVisitedIndexesMapping() override;
    };
    
    template < typename uint_read_len >
    class FASTQReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string id, line, opt_id, quality;
            uint_read_len length;
            std::ifstream* source = 0;
            std::ifstream* pairSource = 0;
            bool ownStreams = false;
            bool pair = false;
            int64_t counter = -1;
            
        public:

            FASTQReadsSourceIterator(const string &srcFile, const string &pairFile = std::string());
            FASTQReadsSourceIterator(std::ifstream* source, std::ifstream* pairSource);

            ~FASTQReadsSourceIterator();

            bool moveNext();
            string getRead();
            string getQualityInfo();
            uint_read_len getReadLength();
            void rewind();

            IndexesMapping* retainVisitedIndexesMapping() override;
    };

    template < typename uint_read_len >
    class RevComplPairReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t counter = -1;

    public:
        RevComplPairReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator);

        bool moveNext();
        string getRead();
        string getQualityInfo();
        uint_read_len getReadLength();
        void rewind();
        IndexesMapping* retainVisitedIndexesMapping() override;
    };

    template < typename uint_read_len >
    class IgnoreNReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t counter = -1;
        vector<uint_reads_cnt_max> indexesMapping;

        bool isFreeOfN();
    public:
        IgnoreNReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator);

        bool moveNext();
        string getRead();
        string getQualityInfo();
        uint_read_len getReadLength();
        void rewind();
        IndexesMapping* retainVisitedIndexesMapping() override;
    };
}

#endif // READSSETITERATOR_H_INCLUDED
