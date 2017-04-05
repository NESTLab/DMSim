cmake_minimum_required(VERSION 2.6.4)

include_directories(${CMAKE_SOURCE_DIR})
add_definitions(-std=c++0x)

add_executable(dmsim main.cpp dmsim.h dmsim.cpp)
