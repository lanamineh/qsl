mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DNUM_QUBITS=12 -DTEST_LEN=4
cmake --build . -j8
echo -e "\nFinished. Binaries are in build/bin/tests/"
