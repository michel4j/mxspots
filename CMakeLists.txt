cmake_minimum_required(VERSION 3.15)
project(mxspots_c LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_library(mxspots SHARED
    csrc/mxspots.c
)

target_include_directories(mxspots PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/csrc>
    $<INSTALL_INTERFACE:include>
)

if (WIN32)
    set_target_properties(mxspots PROPERTIES PREFIX "lib")
endif()

install(TARGETS mxspots DESTINATION mxspots)
