cmake_minimum_required(VERSION 3.10)
project(PgTools)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

set(MATCHER_FILES
        matcher.cpp
        matcher.h)

add_executable(PgMatcher PgMatcher.cpp ${MATCHER_FILES})