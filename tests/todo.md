


# TODO-List for Unit tests

There are several issues in the management of unit tests
that need be addressed in the future of this project.



## Return values

All tests should return 0 on successful completion


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



