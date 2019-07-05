cmake_minimum_required(VERSION 3.13)
project(clas12_analysis)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

include_directories(cpp)

add_executable(clas12_analysis
        cpp/colors.hpp
        cpp/constants.hpp
        cpp/deltat.cpp
        cpp/deltat.hpp
        cpp/filehandeler.hpp
        cpp/histogram.cpp
        cpp/histogram.hpp
        cpp/main.cpp
        cpp/main.hpp
        cpp/physics.cpp
        cpp/physics.hpp
        cpp/Branches.cpp
        cpp/Branches.hpp
        cpp/reaction.cpp
        cpp/reaction.hpp)

target_link_libraries(clas12_analysis ${ROOT_LIBRARIES})
