cmake_minimum_required(VERSION 3.13)
project(hpc_coursework)

set(CMAKE_CXX_STANDARD 11)

add_executable(hpc_coursework Burgers.cpp Model.h Burgers.h Model.cpp main.cpp ParseException.h)