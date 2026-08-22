cmake_minimum_required(VERSION 3.15)
project(mxspots_c LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(OpenMP)

add_library(mxspots SHARED
    csrc/mxspots.c
)

target_include_directories(mxspots PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/csrc>
    $<INSTALL_INTERFACE:include>
)

if(OpenMP_C_FOUND)
    target_link_libraries(mxspots PRIVATE OpenMP::OpenMP_C)
endif()

if(NOT MSVC)
    target_compile_options(mxspots PRIVATE -O3 -march=native)
else()
    target_compile_options(mxspots PRIVATE /O2)
endif()

if (WIN32)
    set_target_properties(mxspots PROPERTIES PREFIX "lib")
endif()

install(TARGETS mxspots DESTINATION mxspots)
