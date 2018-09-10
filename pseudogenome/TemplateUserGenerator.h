#ifndef PGTOOLS_TEMPLATEUSERGENERATOR_H
#define PGTOOLS_TEMPLATEUSERGENERATOR_H

#include "PseudoGenomeBase.h"
namespace PgSAIndex {

    class TemplateUserGenerator {
    private:
        TemplateUserGenerator();

    public:
        virtual ~TemplateUserGenerator();

        template<template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class PseudoGenomeClass> class PseudoGenomeUserClass, class PseudoGenomeUserBase>
        static PseudoGenomeUserBase *generatePseudoGenomeUser(PseudoGenomeBase *pgb);

        template<template< typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, class ReadsListClass> class ReadsListUserClass, class ReadsListUserBase>
        static ReadsListUserBase *generateReadsListUser(PseudoGenomeBase *pgb);

        template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, template<typename uint_read_lenX, typename uint_reads_cntX, typename uint_pg_lenX, class PseudoGenomeClass> class PseudoGenomeUserClass, class PseudoGenomeUserBase>
        static PseudoGenomeUserBase *generatePseudoGenomeUserTemplate(PseudoGenomeBase *pBase);

        template<typename uint_read_len, typename uint_reads_cnt, typename uint_pg_len, template<typename uint_read_lenX, typename uint_reads_cntX, typename uint_pg_lenX, class ReadsListClass> class ReadsListUserClass, class ReadsListUserBase>
        static ReadsListUserBase *generateReadsListUserTemplate(PseudoGenomeBase *pBase);
    };
}


#endif //PGTOOLS_TEMPLATEUSERGENERATOR_H
