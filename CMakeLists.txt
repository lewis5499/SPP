cmake_minimum_required(VERSION 3.10)
project(SPP-EKF VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_EXE_LINKER_FLAGS -static)

include_directories(${PROJECT_SOURCE_DIR}/src)
link_libraries(ws2_32)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    add_definitions( -D_USE_MATH_DEFINES=True)
    add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
else()
    set(CMAKE_CXX_FLAGS "-Wno-unused-result")
endif()

file(GLOB_RECURSE SOURCE
        src/main.cpp
        src/dataIO.cpp
        src/decode.cpp
        src/matrix.cpp
        src/sockets.cpp
        src/solution.cpp)

if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
    message("Debug mode:${CMAKE_C_FLAGS_DEBUG}")
    add_executable(${PROJECT_NAME} ${SOURCE} ./resource/ico.rc)
    add_compile_options(-Wall -Werror -Wstrict-prototypes -Wmissing-prototypes)
elseif(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
    message("Release mode:${CMAKE_C_FLAGS_RELEASE}")
    add_executable(${PROJECT_NAME} ${SOURCE} ./resource/ico.rc)
else()
    message("MODE:${CMAKE_BUILD_TYPE}")
    message("MODE:${CMAKE_C_FLAGS_RELEASE}")
    add_executable(${PROJECT_NAME} ${SOURCE} ./resource/ico.rc)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/log)

set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/bin)

#set(CPACK_PROJECT_NAME ${PROJECT_NAME})
#set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
#set(CPACK_GENERATOR "ZIP")
#include(CPack)

