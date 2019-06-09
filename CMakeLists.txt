cmake_minimum_required(VERSION 3.14)
project(OddDegreeCutSep)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(OddDegreeCutSep
        main.cpp
        OddDegreeCutSep.cpp
        OddDegreeCutSep.h)
