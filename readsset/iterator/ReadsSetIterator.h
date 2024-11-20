#ifndef READSSETITERATOR_H_INCLUDED
#define READSSETITERATOR_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include <vector>
#include "../../utils/helper.h"
#include "../../pgrc/pg-config.h"

namespace PgReadsSet {

    typedef uint_read_len_std uint_read_len_max;

    class IndexesMapping {
    public:
        virtual ~IndexesMapping() = default;

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
        vector<uint_reads_cnt_max> mappingWithGuard;
        uint_reads_cnt_max readsCount;
    public:
        VectorMapping(vector<uint_reads_cnt_max> &&mapping, uint_reads_cnt_max readsCount) :
        mappingWithGuard(std::move(mapping)), readsCount(readsCount) {
            if (mappingWithGuard.empty() || mappingWithGuard.back() != readsCount)
                mappingWithGuard.push_back(readsCount);
        }

        uint_reads_cnt_max getReadOriginalIndex(uint_reads_cnt_max idx) override { return mappingWithGuard[idx]; }
        uint_reads_cnt_max getMappedReadsCount() override { return mappingWithGuard.size() - 1; }
        uint_reads_cnt_max getReadsTotalCount() override { return readsCount; }

        vector<uint_reads_cnt_max> &getMappingVector();

        void saveMapping(string mappingFile);
        static VectorMapping* loadMapping(string mappingFile);
    };

    class SumOfMappings : public IndexesMapping {
    private:
        IndexesMapping* im1 = nullptr;
        IndexesMapping* im2 = nullptr;
        uint_reads_cnt_max idxBeg2 = 0;
        inline IndexesMapping *getIm(uint_reads_cnt_max i) const { return (i < idxBeg2 ? im1 : im2); }
        inline uint_reads_cnt_max getImIdx(uint_reads_cnt_max i) const { return i < idxBeg2 ? i : i - idxBeg2; }
    public:
        SumOfMappings(IndexesMapping *im1, IndexesMapping *im2) : im1(im1), im2(im2),
            idxBeg2(im1->getMappedReadsCount()){}

        uint_reads_cnt_max getReadOriginalIndex(uint_reads_cnt_max idx) override {
            return getIm(idx)->getReadOriginalIndex(getImIdx(idx));
        }

        uint_reads_cnt_max getMappedReadsCount() override {
            return im1->getMappedReadsCount() + im2->getMappedReadsCount();
        }

        uint_reads_cnt_max getReadsTotalCount() override {
            return im1->getReadsTotalCount();
        }
    };

    template < typename uint_read_len >
    class ReadsSourceIteratorTemplate
    {
        private:
            string empty;

        public:

            virtual ~ReadsSourceIteratorTemplate();

            virtual bool moveNext() = 0;
            virtual string& getRead() = 0;
            virtual string& getQualityInfo() { return empty; };
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
            std::istream* source = nullptr;
            int64_t counter = -1;

        public:

            ConcatenatedReadsSourceIterator(std::istream* source);

            ~ConcatenatedReadsSourceIterator() override;

            bool moveNext() override;
            string& getRead() override;
            uint_read_len getReadLength() override;
            void rewind() override;

            IndexesMapping* retainVisitedIndexesMapping() override;
    };
    
    template < typename uint_read_len >
    class FASTAReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = nullptr;
            std::istream* pairSource = nullptr;
            bool pair = false;
            int64_t counter = -1;
            
        public:

            FASTAReadsSourceIterator(std::istream* source);
            
            FASTAReadsSourceIterator(std::istream* source, std::istream* pairSource);

            ~FASTAReadsSourceIterator() override;

            bool moveNext() override;
            string& getRead() override;
            uint_read_len getReadLength() override;
            void rewind() override;

            IndexesMapping* retainVisitedIndexesMapping() override;
    };
    
    template < typename uint_read_len >
    class FASTQReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string id, line, opt_id, quality;
            uint_read_len length = 0;
            std::ifstream* source = nullptr;
            std::ifstream* pairSource = nullptr;
            bool ownStreams = false;
            bool pair = false;
            int64_t counter = -1;
            
        public:

            FASTQReadsSourceIterator(const string &srcFile, const string &pairFile = std::string());
            FASTQReadsSourceIterator(std::ifstream* source, std::ifstream* pairSource);

            ~FASTQReadsSourceIterator() override;

            bool moveNext() override;
            string& getRead() override;
            string& getQualityInfo() override;
            uint_read_len getReadLength() override;
            void rewind() override;

            IndexesMapping* retainVisitedIndexesMapping() override;
    };

    template < typename uint_read_len >
    class RevComplPairReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t counter = -1;
        bool read2Reverse, qual2Reverse, reverseQualityStream;

    public:
        RevComplPairReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator,
                                     bool reverseQualityStream = false);

        bool moveNext() override;
        string& getRead() override;
        string& getQualityInfo() override;
        uint_read_len getReadLength() override;
        void rewind() override;
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

        bool moveNext() override;
        string& getRead() override;
        string& getQualityInfo() override;
        uint_read_len getReadLength() override;
        void rewind() override;
        IndexesMapping* retainVisitedIndexesMapping() override;
    };
}

#endif // READSSETITERATOR_H_INCLUDED
