# Check Methods for Debugging in Computational Science

This note explains the usage of `check` methods in C++. I believe they are a helpful tool
in debugging computational science software. We discuss the concept and some pecularities 
specific to class inheritance in C++. 

Objects in C++ conceptually have an internal state that is represented by their member objects 
and the state of their base classes. Typically, not all conceivable internal states of an object 
are supposed to occur during runtime. However, inconsistent internal states might still arise 
during the development phase (and later) due to faulty implementations of member methods. 

We try to catch such inconsistent states by checking for them explicitly and raising an error
on occurance. Once an inconsistent state has been detected, we abort the program.
This note is concerned with practices in scientific computing, 
where software is used by experts, so we are content with simply aborting on detection of an inconistent state.
This is acceptable for codes that do a few but complicated tasks. 
If your needs differ, you might use exceptions instaed or something entirely different. 

It seems reasonable to summarize these sanity checks in a specific method of each class. 
For the purpose of this note, we call such methods `check`.
There are different reasonable signatures.
For example, `void check() const`.
We want `const` because the check method is not supposed to modify the object.

For example, if one member variable `x` is implemented as an `int`
and supposed to count the number of something, then negative values of `x` would represent an inconsistent
state which is never to occur in practice. Declaring `x` is an `unsigned int` would make such an inconsistent
state non-representable in the first place, but often the solution is not as simple as that.

[Example with matrix]



a counter member is implemented as an `int`, 
then negative numbers are not supposed to occurQuite often, we consider only Only some of those states are 'consistent' and 
