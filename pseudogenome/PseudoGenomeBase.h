#ifndef PSEUDOGENOMEBASE_H_INCLUDED
#define PSEUDOGENOMEBASE_H_INCLUDED

#include "../readsset/ReadsSetBase.h"
#include "../readsset/DefaultReadsSet.h"

using namespace PgSAReadsSet;

namespace PgSAIndex {

    class PseudoGenomeBase: public ReadsSetBase
    {
        protected:
            const uint_pg_len_max length;

            PseudoGenomeBase(uint_pg_len_max length, ReadsSetProperties* properties)
            :   ReadsSetBase(properties),
                length(length)
            { };

            PseudoGenomeBase(uint_pg_len_max length, std::istream& src)
            :   ReadsSetBase(src),
                length(length)
            { };

        public:
            
            virtual ~PseudoGenomeBase() { };
            
            bool isPGLengthStd() { return PgSAIndex::isPGLengthStd(length); };
            bool isPGLengthMax() { return PgSAIndex::isPGLengthMax(length); };

            const uint_pg_len_max getPseudoGenomeLength() { return this->length; };

            virtual string getTypeID() = 0;
            virtual void write(std::ostream& dest) = 0;

            virtual bool validateUsing(DefaultReadsSet* readsSet) {
                fprintf(stderr, "Error: validation is not implemented.\n");
                return false;
            }
            
            virtual void buildRepetitiveReadsFilter() { throw(errno); };

            virtual const string getPseudoGenomeVirtual() = 0;

            const static string PSEUDOGENOME_FILE_SUFFIX;
    };

    class PseudoGenomeHeader {
        private:
            string type;
            bool constantReadLength = true;
            uint_read_len_max maxReadLength;
            uint_reads_cnt_max readsCount;
            uint_pg_len_max pgLength;

        public:

            static const string PSEUDOGENOME_HEADER;

            PseudoGenomeHeader() {}

            PseudoGenomeHeader(PseudoGenomeBase* base) {
                this->type = base->getTypeID();
                this->constantReadLength = base->getReadsSetProperties()->constantReadLength;
                this->maxReadLength = base->getReadsSetProperties()->maxReadLength;
                this->readsCount = base->getReadsSetProperties()->readsCount;
                this->pgLength = base->getPseudoGenomeLength();
            }

            PseudoGenomeHeader(std::istream& src) {

                string line;
                src >> line;
                if (line != PSEUDOGENOME_HEADER)
                    cout << "WARNING: wrong PSEUDOGENOME_HEADER";

                src >> type;
                src >> constantReadLength;
                src >> maxReadLength;
                src >> readsCount;
                src >> pgLength;
                src.get();
            }

            ~PseudoGenomeHeader() { };
            
            void write(std::ostream& dest) {

                dest << PSEUDOGENOME_HEADER << "\n";
                dest << type << "\n";
                dest << constantReadLength << "\n";
                dest << maxReadLength << "\n";
                dest << readsCount << "\n";
                dest << pgLength << "\n";
            }

            string getType() { return this->type; };

            bool isReadLengthConstant() { return constantReadLength; };

            bool isReadLengthMin() { return PgSAReadsSet::isReadLengthMin(maxReadLength); };
            bool isReadLengthStd() { return PgSAReadsSet::isReadLengthStd(maxReadLength); };

            bool isReadsCountStd() { return PgSAReadsSet::isReadsCountStd(readsCount); };
            bool isReadsCountMax() { return PgSAReadsSet::isReadsCountMax(readsCount); };

            bool isPGLengthStd() { return PgSAIndex::isPGLengthStd(pgLength); };
            bool isPGLengthMax() { return PgSAIndex::isPGLengthMax(pgLength); };


            uint_reads_cnt_max getReadsCount() const {
                return readsCount;
            }

            void setReadsCount(uint_reads_cnt_max readsCount) {
                PseudoGenomeHeader::readsCount = readsCount;
            }

            void setPseudoGenomeLength(uint_pg_len_max pgLength) {
                PseudoGenomeHeader::pgLength = pgLength;
            }

            uint_read_len_max getMaxReadLength() const {
                return maxReadLength;
            };

            uint_pg_len_max getPseudoGenomeLength() { return this->pgLength; };

            ReadsSetProperties* generateReadsSetProperties() {
                ReadsSetProperties* rsProp = new ReadsSetProperties();
                rsProp->constantReadLength = this->constantReadLength;
                rsProp->minReadLength = this->constantReadLength?this->maxReadLength:1;
                rsProp->maxReadLength = this->maxReadLength;
                rsProp->allReadsLength = this->constantReadLength?this->readsCount*this->maxReadLength:-1;
                rsProp->readsCount = this->readsCount;
                return rsProp;
            }
    };

}

#endif // PSEUDOGENOMEBASE_H_INCLUDED
