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

# Test install to build directory
#set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/qsl-install-test)
set(CMAKE_INSTALL_PREFIX /usr)

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

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
  )
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(run-gtest src/gtest/utils.cpp src/gtest/sim-constructors.cpp src/gtest/gates.cpp)
target_link_libraries(run-gtest gtest_main qsl armadillo)

include(GoogleTest)
gtest_discover_tests(run-gtest)

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

add_executable(resize examples/resize.cpp)
target_link_libraries(resize PUBLIC qsl)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if(CMAKE_COMPILER_IS_GNUCXX)

  # Output is used by codecov
  include(CodeCoverage)
  #set(COVERAGE_EXCLUDES "dir1/*" "dir2/*")
  setup_target_for_coverage_lcov(
    NAME coverage
    EXECUTABLE ./bin/run-gtest
    BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/catch2/*"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/gtest/*"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/test/*"
    EXCLUDE "${CMAKE_BINARY_DIR}/*"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include/qsl/benchmark/*"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/examples/*"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include/qsl/utils/timer.hpp"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/timer.cpp"
    EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/print.cpp"
    EXCLUDE "/usr/*")
endif()

install(TARGETS qsl LIBRARY DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/qsl DESTINATION include)
