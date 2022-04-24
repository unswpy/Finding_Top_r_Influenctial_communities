cmake_minimum_required(VERSION 3.10)
project(Weight_Community)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated") # avoid hash_map warning

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif()
#
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#add_executable(Preprocess preprocess.cpp)
#add_executable(SPreprocess spreprocess.cpp serialize_qt.h file_ctrl.h)
add_executable(Weight_Community main.cpp stdafx.h argument.h TopR.h TopR.cpp
        serialize.h serialize.cpp)