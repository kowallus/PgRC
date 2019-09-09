#include "PackedPseudoGenome.h"

namespace PgSAIndex {
    
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::PackedPseudoGenome(DefaultPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* srcPseudoGenome, uchar symbolsPerElement)
    : PackedPseudoGenomeBase(srcPseudoGenome->getLength(), srcPseudoGenome->getReadsSetProperties(), symbolsPerElement, sizeof (uint_pg_element)) {
        this->sequence = new uint_pg_element[this->getElementsCountWithGuard()];
        sPacker = new SymbolsPackingFacility(this->getReadsSetProperties(), symbolsPerElement);

        // first element is shifted one symbol...
        sequence[0] = sPacker->packPrefixSymbols(srcPseudoGenome->getSuffix(0), symbolsPerElement - 1);

        uint_max i = sPacker->packSequence(srcPseudoGenome->getSuffix(symbolsPerElement - 1), srcPseudoGenome->getLengthWithGuard(), sequence + 1);
        while (++i < getElementsCountWithGuard())
            sequence[i] = sPacker->getMaxValue();

        this->readsList = srcPseudoGenome->getReadsList();
        srcPseudoGenome->unmanageReadsList();
        
        delete(srcPseudoGenome);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::PackedPseudoGenome(uint_pg_len pgLength, std::istream& src)
    : PackedPseudoGenomeBase(pgLength, src, sizeof (uint_pg_element)) {

        uint_pg_len arraySize;
        src >> arraySize;
        src.get(); // '/n'

        if (arraySize != getElementsCountWithGuard())
            cout << "WARNING: wrong size of pseudogenome.";

        sequence = (uint_pg_element*) PgSAHelpers::readArray(src, getElementsCountWithGuard() * sizeof (uint_pg_element));
        src.get(); // '/n'
        this->readsList = new ReadsListClass(maxReadLength(), src);

        sPacker = new SymbolsPackingFacility(this->getReadsSetProperties(), symbolsPerElement);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::~PackedPseudoGenome() {
        delete[]sequence;
        delete(readsList);
        delete(sPacker);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    void PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::write(std::ostream& dest) {
        //FIXME: delegate to base type
        this->getReadsSetProperties()->write(dest);
        dest << (int) this->symbolsPerElement << "\n";

        dest << getElementsCountWithGuard() << "\n";
        PgSAHelpers::writeArray(dest, sequence, getElementsCountWithGuard() * sizeof (uint_pg_element));
        dest << "\n";
        this->readsList->write(dest);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const uint_pg_element* PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getRawSuffix(const uint_pg_len rawPos) {
        return sequence + rawPos;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const uint_pg_element* PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getRawSuffix(const uint_reads_cnt readsListIdx, const uint_read_len pos) {
        return getRawSuffix((this->readsList->getReadPosition(readsListIdx) + pos + 1) / symbolsPerElement);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const string PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSuffix(const uint_pg_len pos, const uint_pg_len length) {
        return sPacker->reverseSequence(sequence, pos + 1, length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const string PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSuffix(const uint_reads_cnt readsListIdx, const uint_read_len offset, const uint_pg_len length) {
        return getSuffix(this->readsList->getReadPosition(readsListIdx) + offset, length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const string PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getRead(uint_reads_cnt originalIdx) {
        return getSuffix(this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)), readLength(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>* PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::castBase(PseudoGenomeBase* base) {
        // TODO: validate
        return static_cast<PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>*> (base);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_pg_len PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getElementsCountWithGuard() {
        return 2 + (this->length + this->properties->maxReadLength - 1) / symbolsPerElement;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline uint_reads_cnt PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readsCount() {
        return this->properties->readsCount;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    bool PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::isReadLengthConstant() {
        return this->properties->constantReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline uint_read_len PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::maxReadLength() {
        return this->properties->maxReadLength;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    bool PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::isReadLengthConstantVirtual() {
        return isReadLengthConstant();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline void PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getKmerByPosition(const uint_reads_cnt originalIdx, const uint_read_len pos, const uint_read_len kmerLength, char_pg* kmerPtr) {
        sPacker->reverseSequence(sequence, this->readsList->getReadPosition(this->readsList->getReadsListIndexOfOriginalIndex(originalIdx)) + pos + 1, kmerLength, kmerPtr);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, ReadsListClass>* PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getReadsList() {
        return readsList;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline char PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSymbol(uint_pg_len pos) {
        uint_max division = divideBySmallInteger(pos + 1, symbolsPerElement);
        return sPacker->reverseValue(sequence[division], moduloBySmallInteger((uint_max) pos + 1, symbolsPerElement, division));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    SymbolsPackingFacility* PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSymbolsPacker() {
        return sPacker;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    string PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getTypeID() {
        return PGTYPE_PACKED;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline uint_read_len PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readLength(uint_reads_cnt originalIdx) {
        return this->readsList->getReadLength(
                this->readsList->getReadsListIndexOfOriginalIndex(originalIdx));
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const char_pg PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getSymbolImpl(const uint_pg_len posIdx) {
        char tmp;
        sPacker->reverseSequence(sequence, posIdx + 1, 1, &tmp);
        return tmp;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string
    PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getPartImpl(
            const uint_pg_len posIdx, const uint_pg_len length) {
        return sPacker->reverseSequence(sequence, posIdx + 1, length);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    inline const uint_pg_len PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getLengthImpl() {
        return this->getPseudoGenomeLength();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_read_len PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::maxReadLengthVirtual() {
        return maxReadLength();
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_read_len PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readLengthVirtual(uint_reads_cnt i) {
        return readLength(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    uint_reads_cnt PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::readsCountVirtual() {
        return readsCount();
    }
   
    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getReadVirtual(uint_reads_cnt i) {
        return getRead(i);
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, typename uint_pg_element, class ReadsListClass>
    const string PackedPseudoGenome<uint_read_len, uint_reads_cnt, uint_pg_len, uint_pg_element, ReadsListClass>::getPseudoGenomeVirtual() {
        return getSuffix(0, this->getLength());
    }

    template class PackedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class PackedPseudoGenome<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max>::Type>;
    template class PackedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std>::Type>;
    template class PackedPseudoGenome<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max, uint_ps_element_min, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max>::Type>;

}
