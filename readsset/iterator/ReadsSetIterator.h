#ifndef READSSETITERATOR_H_INCLUDED
#define READSSETITERATOR_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include "../../utils/helper.h"
#include "../../pgsaconfig.h"

namespace PgSAReadsSet {

    typedef uint_read_len_std uint_read_len_max;

    template < typename uint_read_len >
    class ReadsSourceIteratorTemplate
    {
        public:

            virtual ~ReadsSourceIteratorTemplate();

            virtual bool moveNextVirtual() = 0;
            virtual string getReadVirtual() = 0;
            virtual string getQualityInfoVirtual() { return std::string(); };
            virtual uint_read_len getReadLengthVirtual() = 0;
            virtual void rewindVirtual() = 0;
    };

    template < typename uint_read_len >
    class ConcatenatedReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string line;
            uint_read_len length;
            std::istream* source = 0;

        public:

            ConcatenatedReadsSourceIterator(std::istream* source);

            ~ConcatenatedReadsSourceIterator();

            bool moveNext();

            string getRead();

            uint_read_len getReadLength();

            void rewind();
            
            bool moveNextVirtual();
            string getReadVirtual();
            uint_read_len getReadLengthVirtual();
            void rewindVirtual();

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
            
        public:

            FASTAReadsSourceIterator(std::istream* source);
            
            FASTAReadsSourceIterator(std::istream* source, std::istream* pairSource);

            ~FASTAReadsSourceIterator();

            bool moveNext();

            string getRead();

            uint_read_len getReadLength();

            void rewind();
            
            bool moveNextVirtual();
            string getReadVirtual();
            uint_read_len getReadLengthVirtual();
            void rewindVirtual();
    };
    
    template < typename uint_read_len >
    class FASTQReadsSourceIterator: public ReadsSourceIteratorTemplate< uint_read_len >
    {
        private:
            std::string id, line, opt_id, quality;
            uint_read_len length;
            std::istream* source = 0;
            std::istream* pairSource = 0;
            bool pair = false;
            
        public:
            
            FASTQReadsSourceIterator(std::istream* source, std::istream* pairSource = std::string());

            ~FASTQReadsSourceIterator();

            bool moveNext();

            string getRead();
            string getQualityInfo();

            uint_read_len getReadLength();

            void rewind();
            
            bool moveNextVirtual();
            string getReadVirtual();
            string getQualityInfoVirtual();
            uint_read_len getReadLengthVirtual();
            void rewindVirtual();
    };

    template < typename uint_read_len >
    class RevComplPairReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t counter = -1;

    public:
        RevComplPairReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator);

        bool moveNextVirtual();
        string getReadVirtual();
        string getQualityInfoVirtual();
        uint_read_len getReadLengthVirtual();
        void rewindVirtual();
    };

    template < typename uint_read_len >
    class IgnoreNReadsSetIterator: public ReadsSourceIteratorTemplate< uint_read_len > {
    private:
        ReadsSourceIteratorTemplate<uint_read_len>* coreIterator;
        int64_t counter = -1;

        bool isFreeOfN();
    public:
        IgnoreNReadsSetIterator(ReadsSourceIteratorTemplate<uint_read_len> *coreIterator);

        bool moveNextVirtual();
        string getReadVirtual();
        string getQualityInfoVirtual();
        uint_read_len getReadLengthVirtual();
        void rewindVirtual();
    };
}

#endif // READSSETITERATOR_H_INCLUDED
