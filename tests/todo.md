


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



## Command line 

The handling of command line arguments will be facilitated 
by a set of functions/classes written precisely for that purpose.

This should be a mere extractor class
and be written in the C-conforming subset of C++.



## Logging 

In the long run, it would be nice to use a logging class 
that allows for prefixes, git version, date, time, etc.

Setting this up will require some careful thinking 
and refactoring of the entire code. 
A reasonable approach would be a replacement
of cout and cerr throughout the entire code 
by new derivations of the stream class
which facilitate more behavior.

In a first step, this is just two streams 
with the some functionality as cout and cerr.

In a second step, more functionality may be added.

The logging classes that I have seen use macros to emulate 
different log streams, their usage looks like 
    LOG << "here is a message";
alternatively, I would like to skip the shift operator alltogether 
and perhaps replace by a macro to read 
    LOG "Here is a message";
The nice thing is that the log messages get accummulated in the data structure 
and only on destruction of the temporary object the message gets actually written
in the actual logging object. Thus one can impose various 
prefixes and postfixes. 




