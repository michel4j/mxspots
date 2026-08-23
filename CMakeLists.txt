cmake_minimum_required(VERSION 3.15)
project(mxspots_c LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(PORTABLE_BUILD "Build portable binaries without host-specific microarchitecture flags" ON)
option(ENABLE_NATIVE_TUNING "Enable -march=native optimizations for the host CPU" OFF)

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

if(MSVC)
    target_compile_options(mxspots PRIVATE /O2)
else()
    target_compile_options(mxspots PRIVATE -O3)
    if(ENABLE_NATIVE_TUNING OR NOT PORTABLE_BUILD)
        target_compile_options(mxspots PRIVATE -march=native)
    endif()
endif()

if (WIN32)
    set_target_properties(mxspots PROPERTIES PREFIX "lib")
endif()

install(TARGETS mxspots DESTINATION mxspots)
