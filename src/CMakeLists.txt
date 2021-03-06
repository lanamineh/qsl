# Utilities source file list
set(UTILS_SRC utils/complex.cpp utils/quantum.cpp
  utils/timer.cpp utils/misc.cpp utils/random.cpp
  utils/print.cpp)

# Qubits object source file list
set(QUBITS_SRC qubits/default/default.cpp qubits/default/gates.cpp qubits/default/measure.cpp qubits/omp/omp.cpp qubits/omp/gates.cpp qubits/omp/measure.cpp qubits/np/np.cpp qubits/np/gates.cpp qubits/np/measure.cpp qubits/resize/resize.cpp qubits/resize/gates.cpp qubits/resize/measure.cpp qubits/omp_np/omp_np.cpp qubits/omp_np/gates.cpp qubits/omp_np/measure.cpp)

# Make the utilities library
add_library(utils ${UTILS_SRC})

# Make the qubits library
find_package(OpenMP REQUIRED)
add_library(qubits ${QUBITS_SRC})
target_link_libraries(qubits PRIVATE OpenMP::OpenMP_CXX)

# Make qsl shared library
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
add_library(qsl SHARED ${QUBITS_SRC} ${UTILS_SRC} ${THREADS_SRC})
target_link_libraries(qsl PRIVATE OpenMP::OpenMP_CXX Threads::Threads)

# Install targets
install(TARGETS qsl LIBRARY DESTINATION lib)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/qsl DESTINATION include)
