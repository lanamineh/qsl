cmake_minimum_required(VERSION 3.12.4) # To get std=c++20

# set the project name
project(qsl-test VERSION 0.1 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# using GCC
# set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -Wextra -g -pg -fmax-errors=10 -fopenmp -fconcepts-diagnostics-depth=5 -fprofile-arcs -ftest-coverage")
# set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall -Wextra -O2 -fopenmp")

add_compile_options(-Wall -Wextra -O2 -fopenmp -g)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
  )
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  qsl-test
  qsl-testing.cpp
  qsl.cpp
)
target_link_libraries(
  qsl-test
  gtest_main
)

add_executable(
  main
  main.cpp
  qsl.cpp
)

add_library(
  qsl
  SHARED
  qsl.cpp
  )

include(GoogleTest)
gtest_discover_tests(qsl-test)
