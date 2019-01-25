cmake_minimum_required(VERSION 3.10)
project(PgTools)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

set(HELPER_FILES
        utils/byteswap.h utils/helper.cpp utils/helper.h)

set(READSSET_FILES
        ${HELPER_FILES}
        readsset/tools
        readsset/tools/SymbolsPackingFacility.cpp
        readsset/tools/SymbolsPackingFacility.h
        readsset/DefaultReadsSet.cpp
        readsset/DefaultReadsSet.h
        readsset/iterator
        readsset/PackedConstantLengthReadsSet.cpp
        readsset/PackedConstantLengthReadsSet.h
        readsset/ReadsSetBase.h
        readsset/ReadsSetInterface.h
        readsset/iterator/ReadsSetIterator.cpp
        readsset/iterator/ReadsSetIterator.h
        readsset/iterator/DivisionReadsSetDecorators.cpp readsset/iterator/DivisionReadsSetDecorators.h
        readsset/DividedPCLReadsSets.cpp readsset/DividedPCLReadsSets.h
        readsset/persistance/ReadsSetPersistence.cpp readsset/persistance/ReadsSetPersistence.h readsset/tools/ReadsSetAnalizer.cpp readsset/tools/ReadsSetAnalizer.h)

set(PG_FILES
        ${READSSET_FILES}
        pgsaconfig.cpp pgsaconfig.h
        pseudogenome/DefaultPseudoGenome.cpp
        pseudogenome/DefaultPseudoGenome.h
        pseudogenome/PackedPseudoGenome.cpp
        pseudogenome/PackedPseudoGenome.h
        pseudogenome/PackedPseudoGenomeBase.h
        pseudogenome/persistence
        pseudogenome/PseudoGenomeBase.h
        pseudogenome/PseudoGenomeInterface.h
        pseudogenome/readslist
        pseudogenome/persistence/PseudoGenomePersistence.cpp
        pseudogenome/persistence/PseudoGenomePersistence.h
        pseudogenome/readslist/iterator
        pseudogenome/readslist/ListOfConstantLengthReads.cpp
        pseudogenome/readslist/ListOfConstantLengthReads.h
        pseudogenome/readslist/ReadsListInterface.h
        pseudogenome/readslist/ReadsListTypes.h
        pseudogenome/readslist/iterator/ReadsListIteratorInterface.h
        pseudogenome/TemplateUserGenerator.h
        pseudogenome/persistence/SeparatedPseudoGenomePersistence.cpp pseudogenome/persistence/SeparatedPseudoGenomePersistence.h
        pseudogenome/readslist/SeparatedExtendedReadsList.cpp pseudogenome/readslist/SeparatedExtendedReadsList.h
        pseudogenome/readslist/iterator/ExtendedReadsListIteratorInterface.h
        pseudogenome/readslist/iterator/ReadsListIteratorExtendedWrapper.h pseudogenome/SeparatedPseudoGenome.cpp pseudogenome/SeparatedPseudoGenome.h)

set(PG_GEN_FILES
        ${PG_FILES}
        pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.h
        pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h
        pseudogenome/generator/PackedPseudoGenomeGenerator.h
        pseudogenome/generator/PseudoGenomeGeneratorBase.h)

set(STRAND_DETECTOR_FILES
        ${PG_FILES}
        strands/DefaultStrandDetector.cpp strands/DefaultStrandDetector.h strands/StrandDetectorBase.h
        strands/AbstractStrandDetector.cpp strands/AbstractStrandDetector.h PgRCManager.cpp PgRCManager.h)

set(HASHES_FILES
        rollinghash/characterhash.h
        rollinghash/cyclichash.h
        rollinghash/generalhash.h
        rollinghash/mersennetwister.cpp
        rollinghash/mersennetwister.h
        rollinghash/rabinkarphash.h
        rollinghash/threewisehash.h
        matching/copmem/Hashes.h
        matching/copmem/metrohash64.cpp matching/copmem/metrohash64.h matching/copmem/platform.h
        matching/copmem/xxhash.c matching/copmem/xxhash.h
        matching/copmem/city.cpp matching/copmem/city.h matching/copmem/citycrc.h)

set(MATCHER_FILES
        ${PG_FILES}
        ${HASHES_FILES}
        matching/ReadsMatchers.cpp matching/ReadsMatchers.h
        matching/ConstantLengthPatternsOnTextHashMatcher.cpp
        matching/ConstantLengthPatternsOnTextHashMatcher.h
        matching/DefaultPgMatcher.cpp matching/DefaultPgMatcher.h
        matching/SimplePgMatcher.cpp matching/SimplePgMatcher.h
        matching/TextMatchers.cpp matching/TextMatchers.h
        matching/copmem/CopMEMMatcher.cpp matching/copmem/CopMEMMatcher.h)

add_executable(PgExtractor PgExtractor.cpp ${PG_FILES})
add_executable(PgGenerator PgGenerator.cpp ${PG_GEN_FILES})
add_executable(PgStrandDetector PgStrandDetector.cpp ${STRAND_DETECTOR_FILES})
add_executable(ReadsDivider ReadsDivider.cpp ${READSSET_FILES})
add_executable(ReadsMatcher ReadsMatcher.cpp ${MATCHER_FILES})
add_executable(PgMatcher PgMatcher.cpp ${MATCHER_FILES})
add_executable(PgRewriter PgRewriter.cpp ${PG_FILES})
add_executable(PgValidator PgValidator.cpp ${PG_FILES})
add_executable(Pgs2Pairs Pgs2Pairs.cpp ${PG_FILES})

set(CMAKE_EXE_LINKER_FLAGS " -static")
add_executable(PgRC PgRC.cpp ${MATCHER_FILES} ${PG_GEN_FILES} PgRCManager.cpp PgRCManager.h)
target_link_libraries(PgRC -static-libgcc -static-libstdc++)