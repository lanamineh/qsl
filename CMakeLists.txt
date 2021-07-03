cmake_minimum_required(VERSION 3.12.4) # To get std=c++20

# set the project name
project(qsl VERSION 0.1 LANGUAGES CXX ASM)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(OpenMP REQUIRED)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # using Clang
    if (MSVC)
        #set(OpenMP_CXX_FLAGS, "-openmp:llvm") # Experimental version of openmp
        #set(OpenMP_CXX_FLAGS "")
        set(CMAKE_CXX_FLAGS "/EHsc /openmp:llvm")  
    else()
        set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -pg -O2 -fopenmp")  
    endif (MSVC)

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # using GCC
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -pg -O2 -fmax-errors=10 -fopenmp -fconcepts-diagnostics-depth=5 -fprofile-arcs -ftest-coverage")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -pg -O2 -fmax-errors=3 -fopenmp")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
  #message( FATAL_ERROR "Microsoft Visual Studio is not supported yet." )  
  #set(CMAKE_CXX_FLAGS "/EHsc /std:c++17")
endif()

include_directories(${PROJECT_SOURCE_DIR}/include/)

# add source directories
add_subdirectory(src EXCLUDE_FROM_ALL) # The main library directory for all the simulators
add_subdirectory(other/tests EXCLUDE_FROM_ALL) # Small speed-testing directory 

# Example
add_executable(example examples/example.cpp)
target_link_libraries(example PUBLIC qsl)

add_executable(benchmark examples/benchmark.cpp)
target_link_libraries(benchmark PUBLIC qsl)

add_executable(verify examples/verify.cpp)
target_link_libraries(verify PUBLIC qsl)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if(CMAKE_COMPILER_IS_GNUCXX)

  # Output is used by codecov
  include(CodeCoverage)
  setup_target_for_coverage_lcov(
    NAME coverage
    EXECUTABLE basic-tests
    BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/catch2")
    
  include(CodeCoverage)
  setup_target_for_coverage_gcovr_html(
    NAME coverage-html
    EXECUTABLE basic-tests
    BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/catch2")

endif()
 
add_executable(basic-tests src/tests/basic-tests.cpp)
target_include_directories(basic-tests PRIVATE src/)
target_link_libraries(basic-tests PRIVATE qsl)
