cmake_minimum_required(VERSION 3.5)
project(PgTools)

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    message("Processing Apple build")
    if(DEFINED ENV{CONDA_PREFIX})
        set(LIBOMP_PATHS "$ENV{CONDA_PREFIX}/lib/libomp.dylib")
        if(NOT EXISTS ${LIBOMP_PATHS})
            message("libomp not found in ${LIBOMP_PATHS}")
            file(GLOB LIBOMP_PATHS $ENV{CONDA_PREFIX}/*/libomp.dylib)
        endif()
        set(CLANG_PATH $ENV{CONDA_PREFIX}/bin/clang)
        if(NOT EXISTS ${CLANG_PATH})
            message("clang not found in ${CLANG_PATH}")
            unset(CLANG_PATH)
        endif()
        set(CLANGXX_PATH $ENV{CONDA_PREFIX}/bin/clang++)
        if(NOT EXISTS ${CLANGXX_PATH})
            message("clang++ not found in ${CLANGXX_PATH}")
            unset(CLANGXX_PATH)
        endif()
    else()
        message("Conda not found")
        message("Proceeding with backup configuration")
        # Backup configuration may require setting proper paths in
        #     LIBOMP_PATHS AND CUSTOM_COMPILER_PATH variables.
        set(LIBOMP_PATHS "/usr/local/opt/libomp/lib/libomp.dylib")
        if(NOT EXISTS ${LIBOMP_PATHS})
            message("libomp not found in ${LIBOMP_PATHS}")
            unset(LIBOMP_PATHS)
            file(GLOB LIBOMP_PATHS /usr/local/opt/*/libomp.dylib)
            if(NOT LIBOMP_PATHS)
                file(GLOB LIBOMP_PATHS /opt/*/libomp.dylib)
            endif()
            if(NOT LIBOMP_PATHS)
                file(GLOB LIBOMP_PATHS ~/*/libomp.dylib)
            endif()
        endif()
        if((NOT CLANG_PATH) OR (NOT CLANGXX_PATH))
            set(LLVM_PATH "/usr/local/opt/llvm/bin")
            set(CUSTOM_COMPILER_PATH "/usr/local/opt/llvm@13/bin")
            find_program(CLANG_PATH clang ${CUSTOM_COMPILER_PATH} ${LLVM_PATH} NO_DEFAULT_PATH)
            find_program(CLANGXX_PATH clang++ ${CUSTOM_COMPILER_PATH} ${LLVM_PATH} NO_DEFAULT_PATH)
        endif()
    endif()
    if(LIBOMP_PATHS)
        list(GET LIBOMP_PATHS 0 LIBOMP_PATH)
        get_filename_component(LIBOMP_DIRECTORY "${LIBOMP_PATH}" DIRECTORY)
        message("Found LIBOMP: ${LIBOMP_DIRECTORY}")
    endif()
    if((NOT CLANG_PATH) OR (NOT CLANGXX_PATH))
        message("Could not find clang or clang++ in specified paths: ${LLVM_PATH}; ${LLVM16_PATH}; ${LLVM14_PATH}; ${BIN_PATH};")
    else()
        set(CMAKE_C_COMPILER ${CLANG_PATH})
        set(CMAKE_CXX_COMPILER ${CLANGXX_PATH})
    endif()
    message("compilers: ${CMAKE_C_COMPILER};${CMAKE_CXX_COMPILER}")
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.10" CACHE STRING "Minimum OS X deployment version")
    get_filename_component(CLANG_DIR "${CLANG_PATH}" DIRECTORY)
    get_filename_component(COMPILER_DIR "${CLANG_DIR}" DIRECTORY)
    message("${CLANG_PATH};${CLANG_DIR};${COMPILER_DIR}")
    set(OPENMP_LIBRARIES "${COMPILER_DIR}/lib")
    set(OPENMP_INCLUDES "${COMPILER_DIR}/include")
    include_directories("${OPENMP_INCLUDES}")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
    if(LIBOMP_DIRECTORY)
        link_directories("${LIBOMP_DIRECTORY}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-rpath ${LIBOMP_DIRECTORY} -L${LIBOMP_DIRECTORY}")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-format")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")

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
        coders/SymbolsPackingFacility.cpp
        coders/SymbolsPackingFacility.h
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
        readsset/persistance/ReadsSetPersistence.cpp readsset/persistance/ReadsSetPersistence.h readsset/tools/ReadsSetAnalyzer.cpp readsset/tools/ReadsSetAnalyzer.h)

set(CODERS_FILES
        coders/CodersLib.h coders/CodersLib.cpp
        coders/lzma/LzmaEnc.h coders/lzma/LzmaEnc.c
        coders/lzma/LzmaDec.h coders/lzma/LzmaDec.c
        coders/lzma/Alloc.h coders/lzma/Alloc.c
        coders/lzma/CpuArch.h coders/lzma/CpuArch.c
        coders/lzma/LzFindMt.h coders/lzma/LzFindMt.c coders/lzma/LzFindOpt.c
        coders/lzma/LzFind.h coders/lzma/LzFind.c
        coders/lzma/Ppmd7.h coders/lzma/Ppmd7.c
        coders/lzma/Ppmd7Enc.c coders/lzma/Ppmd7Dec.c
        coders/lzma/Threads.h coders/lzma/Threads.c
        coders/LzmaCoder.cpp coders/LzmaCoder.h
        coders/PpmdCoder.cpp coders/PpmdCoder.h
        coders/RangeCoder.cpp coders/RangeCoder.h
        coders/fse/entropy_common.c coders/fse/hist.c
        coders/fse/fse_compress.c coders/fse/fse_decompress.c
        coders/fse/huf_compress.c coders/fse/huf_decompress.c
        coders/FSECoder.cpp coders/FSECoder.h
        coders/VarLenDNACoder.cpp coders/VarLenDNACoder.h
        coders/PropsLibrary.cpp coders/PropsLibrary.h)

set(PG_FILES
        ${READSSET_FILES}
        ${CODERS_FILES}
        pgrc/pg-config.cpp pgrc/pg-config.h
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
        pseudogenome/generator/ParallelGreedySwipingPackedOverlapPseudoGenomeGenerator.cpp
        pseudogenome/generator/ParallelGreedySwipingPackedOverlapPseudoGenomeGenerator.h
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

set(PGRC_FILES
        ${PG_GEN_FILES}
        ${MATCHER_FILES}
        pgrc/pgrc-params.h pgrc/pgrc-data.h
        pgrc/pgrc-encoder.cpp pgrc/pgrc-encoder.h
        pgrc/pgrc-decoder.cpp pgrc/pgrc-decoder.h)

link_libraries(pthread)

add_executable(PgRC PgRC.cpp ${PGRC_FILES})
add_executable(PgRC-dev EXCLUDE_FROM_ALL PgRC.cpp ${PGRC_FILES})
if((CYGWIN OR MINGW) AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Windows static build: ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_NAME}")
    target_link_libraries(PgRC -static -static-libgcc -static-libstdc++)
    target_link_libraries(PgRC-dev -static -static-libgcc -static-libstdc++)
endif()
target_compile_definitions(PgRC-dev PUBLIC "-DDEVELOPER_BUILD")

if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|aarch64)")
    add_executable(PgRC-noavx PgRC.cpp ${PGRC_FILES})
    if((CYGWIN OR MINGW) AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
        target_link_libraries(PgRC-noavx -static -static-libgcc -static-libstdc++)
    endif()
    target_compile_options(PgRC PUBLIC -mavx2)
    target_compile_options(PgRC-dev PUBLIC -mavx2)
endif()
