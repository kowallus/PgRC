#include "DefaultPseudoGenome.h"

namespace PgIndex {

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::DefaultPseudoGenome(uint_pg_len pgLength, ReadsSetProperties* properties)
    : PseudoGenomeBase(pgLength, properties) {
        this->sequence = new char_pg[this->getLengthWithGuard()];
    }
    
    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::DefaultPseudoGenome(uint_pg_len pgLength, std::istream& src)
    : PseudoGenomeBase(pgLength, src) {

        uint_pg_len arraySize;
        src >> arraySize;
        src.get(); // '/n'

        if (arraySize != getLengthWithGuard())
            cout << "WARNING: wrong size of pseudogenome.";

        sequence = (char_pg*) PgHelpers::readArray(src, getLengthWithGuard() * sizeof(char_pg));
        src.get(); // '/n'
        this->readsList = new ReadsListClass(maxReadLength(), src);
    }

    template < typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass >
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::~DefaultPseudoGenome() {
        if (readsList)
            delete(readsList);
        delete[]sequence;      
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::castBase(PseudoGenomeBase* base) {
        // TODO: validate
        return static_cast<DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>*> (base);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getTypeID() {
        return PGTYPE_DEFAULT;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    void DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::write(std::ostream& dest) {
        this->getReadsSetProperties()->write(dest);
        dest << getLengthWithGuard() << "\n";
        PgHelpers::writeArray(dest, sequence, getLengthWithGuard() * sizeof(char_pg));
        dest << "\n";
        this->readsList->write(dest);
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getReadsList() 
    {
        return readsList; 
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_pg_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getLengthWithGuard() {
        return this->length + this->properties->maxReadLength;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline char DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSymbol(uint_pg_len pos) {
        return sequence[pos];
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline const char_pg* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSuffix(const uint_pg_len pos) {
        return sequence + pos;
    }
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline const char_pg* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset) {
        return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline const char_pg* DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSuffixPtrByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos) {
        return sequence + this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::maxReadLength() {
        return this->properties->maxReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline uint_reads_cnt DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readsCount() {
        return this->properties->readsCount;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline bool DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::isReadLengthConstant() {
        return this->properties->constantReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline const string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getRead(const uint_reads_cnt originalIdx) {
        uint_pg_len pos = this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
        return string(sequence + pos, readLength(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readLength(const uint_reads_cnt originalIdx) {
        return this->readsList->getReadLength(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline const char_pg DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getSymbolImpl(const uint_pg_len posIdx) {
        return *(sequence + posIdx);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getPartImpl(
            const uint_pg_len posIdx, const uint_pg_len length) {
        return string(sequence + posIdx, length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    inline const uint_pg_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getLengthImpl() {
        return this->getPseudoGenomeLength();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getPseudoGenomeVirtual() {
        return string(this->getSuffix(0), this->getLength());
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    const string DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::getReadVirtual(uint_reads_cnt i) {
        return getRead(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readLengthVirtual(uint_reads_cnt i) {
        return readLength(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::isReadLengthConstantVirtual() {
        return isReadLengthConstant();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_reads_cnt DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::readsCountVirtual() {
        return readsCount();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    uint_read_len DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::maxReadLengthVirtual() {
        return maxReadLength();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass>
    bool DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>::validateUsing(DefaultReadsSet* readsSet) {
        if (readsSet->getReadsSetProperties()->compareWith(this->getReadsSetProperties())) 
        {
            for (uint_reads_cnt i = 0; i < this->readsCount(); i++)
                if (readsSet->getRead(i).compare(this->getRead(i)) != 0) {
                    fprintf(stderr, "%s", (toString(i) + "-th read is incorrect (below original vs pg).\n\n").c_str());
                    fprintf(stderr, "%s", (readsSet->getRead(i) + "\n").c_str());
                    fprintf(stderr, "%s", (this->getRead(i) + "\n").c_str());
                    return false;
                }
            return true;
        }
        fprintf(stderr, "Reads sets are different (below original vs pg).\n\n");
        readsSet->printout();
        this->getReadsSetProperties()->printout();
        return false;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::GeneratedPseudoGenome(uint_pg_len sequenceLength, ReadsSetProperties* properties)
    : DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass >(sequenceLength, properties) {
        this->pos = 0;

        this->genReadsList = new GeneratedReadsListClass(this->properties->maxReadLength, this->properties->readsCount, this->length);

        this->readsList = genReadsList;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::~GeneratedPseudoGenome() {
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    void GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::append(uint_read_len_max length, uint_read_len_max overlap,
                                                uint_reads_cnt_max orgIdx) {
        genReadsList->add(pos, length, orgIdx);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    char_pg* GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::getSequencePtr() const {
        return this->sequence;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    void GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::append(const string& read, uint_read_len length, uint_read_len overlap, uint_reads_cnt orgIdx) {
        genReadsList->add(pos, length, orgIdx);

        uint_read_len len = length - overlap;
        if (len > 0) {
            strncpy(this->sequence + pos, read.data(), len);
            pos += len;
        }
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    void GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::validate() {

        if (pos != this->length)
            cout << "WARNING: Generated " << pos << " pseudogenome instead of " << this->length << "\n";

        genReadsList->validate();

        // adding guard
        for (uint_pg_len i = pos; i < this->getLengthWithGuard(); i++)
            this->sequence[i] = this->properties->symbolsList[this->properties->symbolsCount - 1];
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    void GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::buildRepetitiveReadsFilter() {
        time_checkpoint();

        uint_max filterCount = 0;

        genReadsList->setDuplicateFilterKmerLength(this->duplicateFilterKmerLength);
        
        uint_pg_len* lookup = new uint_pg_len[getTableLength()]();
        
        uint_max lookupIdx = 0;
        for (uint_read_len i = 0; i < this->duplicateFilterKmerLength; i++)
            lookupIdx = lookupIdx * this->properties->symbolsCount + this->properties->symbolOrder[(unsigned char) this->getSymbol(i)];

        uint_max popSymbol[UCHAR_MAX] = {0};
        for(uchar j = 1; j < this->properties->symbolsCount; j++) {
            popSymbol[j] = j;
            for(uchar k = 1; k < this->duplicateFilterKmerLength; k++)
                popSymbol[j] *= this->properties->symbolsCount;
        }

        uint_reads_cnt readsListIndex = 0;
        
        uint_pg_len i = 0;
        while (i < this->length) {
            while ((int_max) genReadsList->getReadPosition(readsListIndex) < (int_max) i + this->duplicateFilterKmerLength - this->maxReadLength())
                readsListIndex++;            
            
            while (lookup[lookupIdx] >= genReadsList->getReadPosition(readsListIndex) + 1) {
                genReadsList->setDuplicateFilterFlag(readsListIndex);
                filterCount++;
                readsListIndex++;
            }
                
            lookup[lookupIdx] = i + 1;
            lookupIdx -= popSymbol[this->properties->symbolOrder[(unsigned char) this->getSymbol(i)]];
            lookupIdx = lookupIdx * this->properties->symbolsCount + this->properties->symbolOrder[(unsigned char) this->getSymbol(i++ + this->duplicateFilterKmerLength)];
        }

        delete[] lookup;
        
        cout << "Found " << filterCount << " reads containing duplicate " << (int) genReadsList->getDuplicateFilterKmerLength()
                << "-mers in " << time_millis() << " msec!\n";
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    uint_max GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::getTableLength() {
        return powuint(this->properties->symbolsCount, this->duplicateFilterKmerLength) + 1;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class GeneratedReadsListClass>
    uint_max GeneratedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, GeneratedReadsListClass>::pgSuffixToTableIdx(const char_pg* suffix) {
        uint_max idx = 0;
        for (uint_read_len i = 0; i < this->duplicateFilterKmerLength; i++)
            idx = idx * this->properties->symbolsCount + this->properties->symbolOrder[(unsigned char) *suffix++];

        return idx;
    }

    template class DefaultPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class DefaultPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class DefaultPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    
    template class GeneratedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class GeneratedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class GeneratedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class GeneratedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;
}
