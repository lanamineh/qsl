

# A program to test whether a flat state vector of doubles
# (rather than a complex type) is better or worse than using
# a complex number type.
add_executable(flat_test flat_test.cpp)
target_link_libraries(flat_test utils)

# A program to test whether prefetching state vector amplitudes
# using the gcc __builtin_prefetch has any effect on speed. 
add_executable(cache_test cache_test.cpp)
target_link_libraries(cache_test utils)

# A program to test whether the use of C-style arrays for the
# state vector has any effect on the speed of the gates.
add_executable(c_test c_test.cpp)
target_link_libraries(c_test utils)

# A program for testing the use of compiler intrinsics in
# gate implementations.
add_executable(intrinsics_test intrinsics_test.cpp)
target_link_libraries(intrinsics_test utils)

# A program to test the phase shift gate.
add_executable(phase_test phase_test.cpp)
target_link_libraries(phase_test utils quest_seq)

# A program to test the speed of Quest.
add_executable(quest_test quest_test.cpp)
target_link_libraries(quest_test utils quest_seq)

# A program to test OpenMP in gate implementations.
find_package(OpenMP REQUIRED)
add_executable(omp_test omp_test.cpp)
target_link_libraries(omp_test PRIVATE utils OpenMP::OpenMP_CXX)

# A program to test whether traversing the state vector
# backwards effects the gate time.
add_executable(reverse_test reverse_test.cpp)
target_link_libraries(reverse_test utils)

# A program to test whether it makes sense to split
# up the real and imaginary parts of the state vector
# (as opposed to interleaving the real and imaginary
# parts).
add_executable(separate_test separate_test.cpp)
target_link_libraries(separate_test utils)

# A program to test the speed of state vectors made from
# std::complex<double> (instead of a custom type).
add_executable(stdcmplx_test stdcmplx_test.cpp)
target_link_libraries(stdcmplx_test utils)

# A program to test the speed of a replacement struct
# for std::complex.
add_executable(struct_test struct_test.cpp)
target_link_libraries(struct_test utils)

# A program for testing the difference between
# for loops and while loops in gate implementations
# (to better understand loop overhead).
add_executable(while_test while_test.cpp)
target_link_libraries(while_test utils)

# A program for testing whether a single for
# loop with masking to obtain the indices is
# better than two nested for loops
add_executable(oneloop_test oneloop_test.cpp)
target_link_libraries(oneloop_test utils)

set_target_properties(flat_test cache_test c_test intrinsics_test phase_test
  quest_test omp_test reverse_test stdcmplx_test struct_test while_test
  separate_test oneloop_test
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests/"
)
