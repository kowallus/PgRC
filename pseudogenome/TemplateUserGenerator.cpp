#include "TemplateUserGenerator.h"

#include "../strands/DefaultStrandDetector.h"
#include "DefaultPseudoGenome.h"
#include "PackedPseudoGenome.h"

namespace PgSAIndex {

    TemplateUserGenerator::TemplateUserGenerator() {
    }

    TemplateUserGenerator::~TemplateUserGenerator() {
    }

    template<template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass> class PseudoGenomeUserClass, class PseudoGenomeUserBase>
    PseudoGenomeUserBase *TemplateUserGenerator::generatePseudoGenomeUser(PseudoGenomeBase *pgb) {
        if (pgb->isReadLengthMin()) {
            if (pgb->isReadsCountStd()) {
                if (pgb->isPGLengthStd())
                    return generatePseudoGenomeUserTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
                    PseudoGenomeUserClass, PseudoGenomeUserBase>(pgb);
                if (pgb->isPGLengthMax())
                    return generatePseudoGenomeUserTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max,
                            PseudoGenomeUserClass, PseudoGenomeUserBase>(pgb);
            }
        }
        if (pgb->isReadLengthStd()) {
            if (pgb->isReadsCountStd()) {
                if (pgb->isPGLengthStd())
                    return generatePseudoGenomeUserTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std,
                            PseudoGenomeUserClass, PseudoGenomeUserBase>(pgb);
                if (pgb->isPGLengthMax())
                    return generatePseudoGenomeUserTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max,
                            PseudoGenomeUserClass, PseudoGenomeUserBase>(pgb);
            }
        }

        cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
        return 0;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, template<typename uint_read_lenX, typename uint_reads_cntX, typename uint_pg_lenX, class PseudoGenomeClass> class PseudoGenomeUserClass, class PseudoGenomeUserBase>
    PseudoGenomeUserBase *TemplateUserGenerator::generatePseudoGenomeUserTemplate(PseudoGenomeBase *pgb) {
        if (pgb->getTypeID() == PGTYPE_DEFAULT) {
            if (pgb->isReadLengthConstant()) {
                auto pseudoGenome = (DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>*) pgb;
                return new PseudoGenomeUserClass<uint_read_len, uint_reads_cnt, uint_pg_len, DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>>(pseudoGenome);
            }
        }

        if (pgb->getTypeID() == PGTYPE_PACKED) {
            PackedPseudoGenomeBase* ppgb = (PackedPseudoGenomeBase*) pgb;
            if (pgb->isReadLengthConstant()) {
                if (ppgb->isPgElementMinimal()) {
                    auto pseudoGenome = (PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min> *) ppgb;
                    return new PseudoGenomeUserClass<uint_read_len, uint_reads_cnt, uint_pg_len, PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min>>(pseudoGenome);
                }
                if (ppgb->isPgElementStandard()) {
                    auto pseudoGenome = (PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std> *) ppgb;
                    return new PseudoGenomeUserClass<uint_read_len, uint_reads_cnt, uint_pg_len, PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std>>(pseudoGenome);
                }
            }
        }
        cout << "ERROR: unsupported PGTYPE " << pgb->getTypeID() << " with " << (pgb->isReadLengthConstant()?"constant":"variable") << " read length\n";
        return 0;
    }

    template<template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass> class ReadsListUserClass, class ReadsListUserBase>
    ReadsListUserBase *TemplateUserGenerator::generateReadsListUser(PseudoGenomeBase *pgb) {
        if (pgb->isReadLengthMin()) {
            if (pgb->isReadsCountStd()) {
                if (pgb->isPGLengthStd())
                    return generateReadsListUserTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_std,
                            ReadsListUserClass, ReadsListUserBase>(pgb);
                if (pgb->isPGLengthMax())
                    return generateReadsListUserTemplate<uint_read_len_min, uint_reads_cnt_std, uint_pg_len_max,
                            ReadsListUserClass, ReadsListUserBase>(pgb);
            }
        }
        if (pgb->isReadLengthStd()) {
            if (pgb->isReadsCountStd()) {
                if (pgb->isPGLengthStd())
                    return generateReadsListUserTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_std,
                            ReadsListUserClass, ReadsListUserBase>(pgb);
                if (pgb->isPGLengthMax())
                    return generateReadsListUserTemplate<uint_read_len_std, uint_reads_cnt_std, uint_pg_len_max,
                            ReadsListUserClass, ReadsListUserBase>(pgb);
            }
        }

        cout << "ERROR: CANNOT DETERMINE TEMPLATE TYPES";
        return 0;
    }

    template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, template<typename uint_read_lenX, typename uint_reads_cntX, typename uint_pg_lenX, class ReadsListClass> class ReadsListUserClass, class ReadsListUserBase>
    ReadsListUserBase *TemplateUserGenerator::generateReadsListUserTemplate(PseudoGenomeBase *pgb) {
        ReadsListInterface<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type> * readsList = 0;

        if (pgb->getTypeID() == PGTYPE_DEFAULT) {
            if (pgb->isReadLengthConstant())
                readsList = (((DefaultPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len>*) pgb)->getReadsList());
        }

        if (pgb->getTypeID() == PGTYPE_PACKED) {
            PackedPseudoGenomeBase* ppgb = (PackedPseudoGenomeBase*) pgb;
            if (pgb->isReadLengthConstant()) {
                if (ppgb->isPgElementMinimal())
                    readsList = ((PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_min> *) ppgb)->getReadsList();
                if (ppgb->isPgElementStandard())
                    readsList = ((PackedPseudoGenomeOfConstantLengthReadsType<uint_read_len, uint_reads_cnt, uint_pg_len, uint_ps_element_std> *) ppgb)->getReadsList();
            }
        }

        if (readsList)
            return new ReadsListUserClass<uint_read_len, uint_reads_cnt, uint_pg_len, typename ListOfConstantLengthReadsTypeTemplate<uint_read_len, uint_reads_cnt, uint_pg_len>::Type>(readsList);

        cout << "ERROR: unsupported PGTYPE " << pgb->getTypeID() << " with " << (pgb->isReadLengthConstant()?"constant":"variable") << " read length\n";
        return 0;
    }

    using namespace PgTools;

    template StrandDetectorBase* TemplateUserGenerator::generateReadsListUser<DefaultStrandDetector, StrandDetectorBase>(PseudoGenomeBase* pgb);

}