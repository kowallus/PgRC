#ifndef PGTOOLS_DIVIDEDPCLREADSSETS_H
#define PGTOOLS_DIVIDEDPCLREADSSETS_H

#include "PackedConstantLengthReadsSet.h"

using namespace PgReadsSet;

namespace PgTools {

    class DividedPCLReadsSets {
    private:
        PackedConstantLengthReadsSet* hqReadsSet = 0;
        PackedConstantLengthReadsSet* lqReadsSet = 0;
        bool nReadsLQ;
        PackedConstantLengthReadsSet* nReadsSet = 0;
        bool separateNReadsSet;

        VectorMapping* lqMapping = 0;
        VectorMapping* nMapping = 0;

    public:
        DividedPCLReadsSets(uint_read_len_max readLength, bool separateNReadsSet = false, bool nReadsLQ = false);

        virtual ~DividedPCLReadsSets();

        PackedConstantLengthReadsSet *getHqReadsSet() const { return hqReadsSet; };

        PackedConstantLengthReadsSet *getLqReadsSet() const { return lqReadsSet; };
        VectorMapping* getLqReadsIndexesMapping() const { return lqMapping; };
        bool areNReadsLQ() const { return nReadsLQ; };

        PackedConstantLengthReadsSet *getNReadsSet() const { return nReadsSet; };
        VectorMapping* getNReadsIndexesMapping() const { return nMapping; };

        bool isSeparateNReadsSet() const { return separateNReadsSet; };

        static DividedPCLReadsSets* getQualityDivisionBasedReadsSets(
                ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt, uint_read_len_max readLength,
                double error_limit, bool simplified_suffix_mode, bool separateNReadsSet = false, bool nReadsLQ = false);

        static DividedPCLReadsSets *
    getSimpleDividedPCLReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt, uint_read_len_max readLength,
                                 bool separateNReadsSet, bool nReadsLQ);

        static DividedPCLReadsSets *
    loadDivisionReadsSets(ReadsSourceIteratorTemplate<uint_read_len_max> *readsIt, uint_read_len_max readLength,
                          string lqDivisionFile, bool nReadsLQ, string nDivisionFile = "", bool skipHQReadsSet = false);

        void moveLqReadsFromHqReadsSetsToLqReadsSets(const vector<bool> &isReadHqInHqReadsSet);

        IndexesMapping *generateHqReadsIndexesMapping();

        void removeReadsFromLqReadsSet(const vector<bool> &isLqReadMappedIntoHqPg);
        void removeReadsFromNReadsSet(const vector<bool> &isReadMappedIntoHqPg,
                                                           uint_reads_cnt_max nBegIdx);

        void disposeLqReadsSet();
        void disposeNReadsSet();
        void disposeHqReadsSet();
    };

}


#endif //PGTOOLS_DIVIDEDPCLREADSSETS_H
