This section is all about verifying that the simulator works correctly, and analysing its performance.

Verification
############

A single class, ``Verify``, checks that the simulators behave correctly. There are two types of verification that can be performed:

* Comparison of the results of a simulator with a known true result
* Comparison of the results from two different simulators
    
The first test is a stronger test for correctness, but the second check is often sufficient (iff all the simulators agree, then it is likely they are all correct).

Some of the operations that the simulators perform are deterministic, while others aren't. For the deterministic operations, such as gate applications, it is possible to check that the result is equal to another given result. For probabilistic operations (such as measurement, which involves a random collapse of the state vector), it is necessary to perform a large number of identical tests and check that the average result is close to another given result.


Benchmarking
############
.. doxygenfile:: bench.hpp
