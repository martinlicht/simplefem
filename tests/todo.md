
# TODO-List for Unit tests

There are several issues in the management of unit tests
that need be addressed in the future of this project.



## Return values

All tests should return 0 on successful completion


## structure of tests 

The tests should have a similar form. 
For example, they should state the test topic and when they are done.
This could be handled by macros.

    OPEN "Test for correctness of some feature";
    INFO "Some statement";
    WARNING "something has happened";
    ERROR "error has occured";
    CLOSE "Test finished";
    
There is a difference between output of the program state 
and output of data. The former occurs at sequence points 
and typically does not depend on the parameters of the computation,
it is fixed, so to say. The latter, by contrast, typically 
involves the data and is only written to demonstrate the results 
of the computation. It's length does depend on data.
The logging output is of interest to the general debugging,
whereas the data output is only helpful for visual inspection.


## Global structure

It makes sense to order the tests mimicking the source directory:
./unittests/subtopic/classname.cpp 

This allows to simplify the file names of the tests.

The makefile should remain within the test directory.


## Global testing routine 

How can we organize the unit tests such that 
the testing is automatized. Maybe all the tests are run
and errors are just collected.


## valgrind and other flags 

if tests are automated, then they should include flags 
such as valgrind.


## Virtual files 

In some circumstances, such as testing the input/output routines,
it seems reasonable to avoid actual harddrive files in favor of 
purely virtual files. 

It seems natural that the input/routines 
should be completely oblivious to whether 
the streams are files at all or mere stringstreams.



## Control of output verbosity 

The verbosity of the output should be managed explicitly.
Typically, the test itself should eject little to no output 
depending on the global test format (see below)

It would be nice to have a command line switch to increase verbosity (such as -v).



