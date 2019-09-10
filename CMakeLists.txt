cmake_minimum_required(VERSION 3.4)
project(PgTools)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2")

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(HELPER_FILES
        utils/byteswap.h utils/helper.cpp utils/helper.h)

set(READSSET_FILES
        ${HELPER_FILES}
        readsset/tools
        utils/SymbolsPackingFacility.cpp
        utils/SymbolsPackingFacility.h
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

set(LZMA_FILES
        utils/LzmaLib.h utils/LzmaLib.cpp
        lzma/LzmaEnc.h lzma/LzmaEnc.c
        lzma/LzmaDec.h lzma/LzmaDec.c
        lzma/Alloc.h lzma/Alloc.c
        lzma/LzFindMt.h lzma/LzFindMt.c
        lzma/LzFind.h lzma/LzFind.c
        lzma/Ppmd7.h lzma/Ppmd7.c
        lzma/Ppmd7Enc.c lzma/Ppmd7Dec.c
        lzma/Threads.h lzma/Threads.c
        utils/VarLenDNACoder.cpp utils/VarLenDNACoder.h)

set(PG_FILES
        ${READSSET_FILES}
        ${LZMA_FILES}
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
        pseudogenome/SeparatedPseudoGenomeBase.cpp pseudogenome/SeparatedPseudoGenomeBase.h
        pseudogenome/SeparatedPseudoGenome.cpp pseudogenome/SeparatedPseudoGenome.h
        pseudogenome/readslist/SeparatedExtendedReadsList.cpp pseudogenome/readslist/SeparatedExtendedReadsList.h
        pseudogenome/persistence/SeparatedPseudoGenomePersistence.cpp pseudogenome/persistence/SeparatedPseudoGenomePersistence.h
        pseudogenome/readslist/iterator/ExtendedReadsListIteratorInterface.h
        pseudogenome/readslist/iterator/ReadsListIteratorExtendedWrapper.h)

set(PG_GEN_FILES
        ${PG_FILES}
        pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/AbstractOverlapPseudoGenomeGenerator.h
        pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/GreedySwipingPackedOverlapPseudoGenomeGenerator.h
        pseudogenome/generator/PackedPseudoGenomeGenerator.h
        pseudogenome/generator/PseudoGenomeGeneratorBase.h)

set(HASHES_FILES
        matching/rollinghash/characterhash.h
        matching/rollinghash/cyclichash.h
        matching/rollinghash/generalhash.h
        matching/rollinghash/mersennetwister.cpp
        matching/rollinghash/mersennetwister.h
        matching/rollinghash/rabinkarphash.h
        matching/rollinghash/threewisehash.h
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

add_executable(PgRC PgRC.cpp PgRCManager.cpp PgRCManager.h ${MATCHER_FILES} ${PG_GEN_FILES})
TARGET_LINK_LIBRARIES(PgRC pthread)
add_executable(PgRC-win PgRC.cpp PgRCManager.cpp PgRCManager.h ${MATCHER_FILES} ${PG_GEN_FILES})
TARGET_LINK_LIBRARIES(PgRC-win pthread)
if((CYGWIN OR MINGW) AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Windows static build: ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_NAME}")
    target_link_libraries(PgRC-win  -static -static-libgcc -static-libstdc++)
endif()
add_executable(PgRC-dev PgRC.cpp PgRCManager.cpp PgRCManager.h ${MATCHER_FILES} ${PG_GEN_FILES})
target_compile_definitions(PgRC-dev PUBLIC "-DDEVELOPER_BUILD")
TARGET_LINK_LIBRARIES(PgRC-dev pthread)