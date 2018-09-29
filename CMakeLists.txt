cmake_minimum_required(VERSION 3.10)
project(PgTools)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

set(HELPER_FILES
        byteswap.h helper.cpp helper.h pghelper.cpp pghelper.h)

set(READSSET_FILES
        ${HELPER_FILES}
        readsset/DefaultReadsSet.cpp
        readsset/DefaultReadsSet.h
        readsset/iterator
        readsset/PackedReadsSet.cpp
        readsset/PackedReadsSet.h
        readsset/ReadsSetBase.h
        readsset/ReadsSetInterface.h
        readsset/iterator/ReadsSetIterator.cpp
        readsset/iterator/ReadsSetIterator.h)

set(PG_FILES
        ${READSSET_FILES}
        pgsaconfig.cpp pgsaconfig.h
        pseudogenome/DefaultPseudoGenome.cpp
        pseudogenome/DefaultPseudoGenome.h
        pseudogenome/PackedPseudoGenome.cpp
        pseudogenome/PackedPseudoGenome.h
        pseudogenome/PackedPseudoGenomeBase.h
        pseudogenome/packing
        pseudogenome/persistence
        pseudogenome/PseudoGenomeBase.h
        pseudogenome/PseudoGenomeInterface.h
        pseudogenome/readslist
        pseudogenome/packing/SymbolsPackingFacility.cpp
        pseudogenome/packing/SymbolsPackingFacility.h
        pseudogenome/persistence/PseudoGenomePersistence.cpp
        pseudogenome/persistence/PseudoGenomePersistence.h
        pseudogenome/readslist/iterator
        pseudogenome/readslist/ListOfConstantLengthReads.cpp
        pseudogenome/readslist/ListOfConstantLengthReads.h
        pseudogenome/readslist/ReadsListInterface.h
        pseudogenome/readslist/ReadsListTypes.h
        pseudogenome/readslist/iterator/ReadsListIteratorInterface.h
        pseudogenome/TemplateUserGenerator.h)

set(PG_GEN_FILES
        ${PG_FILES}
        pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.h
        pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h
        pseudogenome/generator/PackedPseudoGenomeGenerator.h
        pseudogenome/generator/PseudoGenomeGeneratorBase.h
        pseudogenome/persistence/SeparatedPseudoGenomePersistence.cpp pseudogenome/persistence/SeparatedPseudoGenomePersistence.h)

set(STRAND_DETECTOR_FILES
        ${PG_FILES}
        strands/DefaultStrandDetector.cpp strands/DefaultStrandDetector.h strands/StrandDetectorBase.h strands/AbstractStrandDetector.cpp strands/AbstractStrandDetector.h)

set(ROLLING_HASH_FILES
        rollinghash/characterhash.h
        rollinghash/cyclichash.h
        rollinghash/generalhash.h
        rollinghash/mersennetwister.cpp
        rollinghash/mersennetwister.h
        rollinghash/rabinkarphash.h
        rollinghash/threewisehash.h)

set(MATCHER_FILES
        ${PG_FILES}
        ${ROLLING_HASH_FILES}
        matching/matcher.cpp matching/matcher.h
        matching/ConstantLengthPatternsOnTextHashMatcher.cpp
        matching/ConstantLengthPatternsOnTextHashMatcher.h
        matching/PgMatcherBase.h matching/DefaultPgMatcher.cpp matching/DefaultPgMatcher.h)

add_executable(PgMatcher PgMatcher.cpp ${MATCHER_FILES})
add_executable(PgExtractor PgExtractor.cpp ${PG_FILES})
add_executable(PgGenerator PgGenerator.cpp ${PG_GEN_FILES})
add_executable(PgStrandDetector PgStrandDetector.cpp ${STRAND_DETECTOR_FILES})