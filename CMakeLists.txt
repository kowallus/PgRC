cmake_minimum_required(VERSION 3.10)
project(PgTools)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

set(HELPER_FILES
        byteswap.h helper.cpp helper.h pghelper.cpp pghelper.h)

set(MATCHER_FILES
        ${HELPER_FILES}
        matcher.cpp matcher.h)

add_executable(PgMatcher PgMatcher.cpp ${MATCHER_FILES})
add_executable(PgExtractor PgExtractor.cpp ${HELPER_FILES})
add_executable(PgStrandDetector PgStrandDetector.cpp ${HELPER_FILES})