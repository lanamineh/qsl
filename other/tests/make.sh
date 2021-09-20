rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DNUM_QUBITS=12 -DTEST_LEN=20000
cmake --build . -j8
echo -e "\nFinished. Binaries are in build/bin/tests/"
