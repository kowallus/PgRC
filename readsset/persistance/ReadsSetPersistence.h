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
            ReadsSourceIteratorTemplate< uint_read_len_max>* coreIterator = 0;
            ReadsSourceIteratorTemplate< uint_read_len_max>* readsIterator = 0;

            ifstream* srcSource = 0;
            ifstream* pairSource = 0;
            ifstream* divSource = 0;
        public:
            ManagedReadsSetIterator(string srcFile, string pairFile = "", string divisionFile = "", bool divisionComplement = false);

        public:
            bool moveNextVirtual();
            string getReadVirtual();
            string getQualityInfoVirtual();
            uint_read_len_max getReadLengthVirtual();
            void rewindVirtual();

            virtual ~ManagedReadsSetIterator();
        };

    public:
        static ReadsSourceIteratorTemplate<uint_read_len_max>* createManagedReadsIterator(string srcFile,
                                                                                          string pairFile = "",
                                                                                          string divisionFile = "",
                                                                                          bool divisionComplement = false);



    };

}


#endif //PGTOOLS_READSSETPERSISTENCE_H
