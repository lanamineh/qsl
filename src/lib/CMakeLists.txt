# Utilities source file list
set(UTILS_SRC utils/complex.cpp utils/quantum.cpp
  utils/timer.cpp utils/misc.cpp utils/random.cpp)

# Threads source file list
set(THREADS_SRC threads/joiner.cpp threads/pool.cpp)

# Qubits object source file list
set(QUBITS_SRC qubits/default/default.cpp qubits/omp/omp.cpp qubits/np/np.cpp
qubits/qubits_thread.cpp qubits/qubits_asm.cpp qubits/x86_64/pauliX.S qubits/x86_64/rotateX.S)

# Make the utilities library
add_library(utils ${UTILS_SRC})

# Make the threads library
find_package(Threads REQUIRED)
add_library(threads ${THREADS_SRC})
target_link_libraries(threads PRIVATE Threads::Threads)

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
