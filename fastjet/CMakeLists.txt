# CMake file for fastjet codes used for benchmarks and validation
#
# (C) CERN, 2024
# MIT License

cmake_minimum_required(VERSION 3.12)
project(FastJetBenchmarks)

# Add package utilities to CMake path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

find_package(HepMC3 REQUIRED)
# find_package(FastJet REQUIRED)

set(CMAKE_CXX_STANDARD 17)

# Fastjet executable runs reconstruction then outputs either
# exclusive or inclusive jets
add_executable(fastjet-finder
    src/fastjet-finder.cc
    src/fastjet-utils.cc
)

target_include_directories(fastjet-finder PRIVATE
    ../../fastjet-install/include
)

target_link_libraries(fastjet-finder 
    HepMC3::HepMC3
    ../../fastjet-install/lib/libfastjet.so
)
