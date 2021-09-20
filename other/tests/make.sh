mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DNUM_QUBITS=12 -DTEST_LEN=20000
cmake --build . -j8
echo -e "\nFinished. Binaries are in build/bin/tests/"
echo -e "\nstd::complex test"
./bin/tests/stdcmplx_test
echo -e "\nstruct test"
./bin/tests/struct_test
echo -e "\nflat test"
./bin/tests/flat_test
echo -e "\nseparate test"
./bin/tests/separate_test
