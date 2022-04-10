/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://www.hdfgroup.org/licenses.               *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 *  This example illustrates how to create a dataset that is a 4 x 6
 *  array. It is used in the HDF5 Tutorial.
 */

#include <iostream>
#include <string>
#include <concepts>

#include "H5Cpp.h"

// Removed using namespace H5. It is not clear when something is in this
// namespace and when it is in the global namespace.

const std::string filename{"test.h5"};
const std::string dataset_name{"dset"};
const int x{4}; // dataset dimensions
const int y{6};
const int rank{2};

template<std::floating_point A>
class Thingy{};

int main()
{
    
    try {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        H5::Exception::dontPrint();

        // Create a new file using the default property lists.
	H5::H5File file{filename, H5F_ACC_TRUNC};

        // Create the data space for the dataset.
        hsize_t dims[2]; // dataset dimensions
        dims[0] = x;
        dims[1] = y;
        H5::DataSpace dataspace{rank, dims};

        // Create the dataset.
        H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::STD_I32BE, dataspace);

    } catch (const H5::FileIException & error) {
	// catch failure caused by the H5File operations
        error.printErrorStack();
        return -1;
    } catch (const H5::DataSetIException & error) {
	// catch failure caused by the DataSet operations
        error.printErrorStack();
        return -1;
    } catch (const H5::DataSpaceIException & error) {
	// catch failure caused by the DataSpace operations
        error.printErrorStack();
        return -1;
    }

    return 0; // successfully terminated
}
