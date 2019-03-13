cmake_minimum_required(VERSION 3.13)
project(clas12_analysis)

set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

include_directories(cpp)
include_directories(cpp_hipo/hipo)
include_directories(cpp_hipo/libcpp)
include_directories(cpp_hipo/lz4/lib)
include_directories(cpp_hipo/lz4/programs)
include_directories(cpp_reindex)
include_directories(speedTest/cpp)



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
        cpp/physics.hpp)
target_link_libraries(clas12_analysis ${ROOT_LIBRARIES})

add_executable(clas12_hipo
        cpp_hipo/hipo/colors.hpp
        cpp_hipo/hipo/constants.hpp
        cpp_hipo/hipo/datahandeler.cpp
        cpp_hipo/hipo/datahandeler.hpp
        cpp_hipo/hipo/deltat.cpp
        cpp_hipo/hipo/deltat.hpp
        cpp_hipo/hipo/filehandeler.hpp
        cpp_hipo/hipo/glob_files.cpp
        cpp_hipo/hipo/glob_files.hpp
        cpp_hipo/hipo/histogram.cpp
        cpp_hipo/hipo/histogram.hpp
        cpp_hipo/hipo/main.cpp
        cpp_hipo/hipo/main.hpp
        cpp_hipo/hipo/physics.cpp
        cpp_hipo/hipo/physics.hpp)

add_definitions(-D__LZ4__)
add_library(hipo3cpp
        cpp_hipo/libcpp/data.cpp
        cpp_hipo/libcpp/data.h
        cpp_hipo/libcpp/dictionary.cpp
        cpp_hipo/libcpp/dictionary.h
        cpp_hipo/libcpp/event.cpp
        cpp_hipo/libcpp/event.h
        cpp_hipo/libcpp/hipoexceptions.h
        cpp_hipo/libcpp/node.cpp
        cpp_hipo/libcpp/node.h
        cpp_hipo/libcpp/reader.cpp
        cpp_hipo/libcpp/reader.h
        cpp_hipo/libcpp/record.cpp
        cpp_hipo/libcpp/record.h
        cpp_hipo/libcpp/text.cpp
        cpp_hipo/libcpp/text.h
        cpp_hipo/libcpp/utils.cpp
        cpp_hipo/libcpp/utils.h
        cpp_hipo/libcpp/wrapper.cpp
        cpp_hipo/libcpp/writer.cpp
        cpp_hipo/libcpp/writer.h
        cpp_hipo/lz4/contrib/gen_manual/gen_manual.cpp
        cpp_hipo/lz4/examples/blockStreaming_doubleBuffer.c
        cpp_hipo/lz4/examples/blockStreaming_lineByLine.c
        cpp_hipo/lz4/examples/blockStreaming_ringBuffer.c
        cpp_hipo/lz4/examples/compress_functions.c
        cpp_hipo/lz4/examples/dictionaryRandomAccess.c
        cpp_hipo/lz4/examples/frameCompress.c
        cpp_hipo/lz4/examples/HCStreaming_ringBuffer.c
        cpp_hipo/lz4/examples/printVersion.c
        cpp_hipo/lz4/examples/simple_buffer.c
        cpp_hipo/lz4/lib/lz4.c
        cpp_hipo/lz4/lib/lz4.h
        cpp_hipo/lz4/lib/lz4frame.c
        cpp_hipo/lz4/lib/lz4frame.h
        cpp_hipo/lz4/lib/lz4frame_static.h
        cpp_hipo/lz4/lib/lz4hc.c
        cpp_hipo/lz4/lib/lz4hc.h
        cpp_hipo/lz4/lib/xxhash.c
        cpp_hipo/lz4/lib/xxhash.h
        cpp_hipo/lz4/programs/bench.c
        cpp_hipo/lz4/programs/bench.h
        cpp_hipo/lz4/programs/datagen.c
        cpp_hipo/lz4/programs/datagen.h
        cpp_hipo/lz4/programs/lz4cli.c
        cpp_hipo/lz4/programs/lz4io.c
        cpp_hipo/lz4/programs/lz4io.h
        cpp_hipo/lz4/programs/platform.h
        cpp_hipo/lz4/programs/util.h
        cpp_hipo/lz4/tests/checkTag.c
        cpp_hipo/lz4/tests/datagencli.c
        cpp_hipo/lz4/tests/frametest.c
        cpp_hipo/lz4/tests/fullbench.c
        cpp_hipo/lz4/tests/fuzzer.c
        cpp_hipo/lz4/tests/roundTripTest.c)
target_link_libraries(clas12_hipo hipo3cpp ${ROOT_LIBRARIES})

add_executable(clas12_reindex
        cpp_reindex/clas12branches.cpp
        cpp_reindex/clas12branches.hpp
        cpp_reindex/colors.hpp
        cpp_reindex/constants.hpp
        cpp_reindex/deltat.cpp
        cpp_reindex/deltat.hpp
        cpp_reindex/histogram.cpp
        cpp_reindex/histogram.hpp
        cpp_reindex/main.cpp
        cpp_reindex/main.hpp
        cpp_reindex/physics.cpp
        cpp_reindex/physics.hpp
        cpp_reindex/reaction.cpp
        cpp_reindex/reaction.hpp)
target_link_libraries(clas12_reindex ${ROOT_LIBRARIES})

add_executable(speed_test
        speedTest/cpp/main.cpp
        speedTest/cpp/main.hpp)
target_link_libraries(speed_test ${ROOT_LIBRARIES})
