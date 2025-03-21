CODING STYLE GUIDE
======================



Naming Conventions
================================

Use `CapitalizedCamelCase` for classes, structs, enums and namespaces.

Use `CapitalizedCamelCase`, `camelCase` or `snake_case` for functions.

Use `snake_case` for variables.

Underscores also be used in the names of classes, structs, enums, and functions, if they indicate variations of a common principle, such as different implementations of the same idea.

As by common and nearly universal convention, macros must always be in `UPPERCASE_WITH_UNDERSCORES`.

In general, do not abbreviate when naming classes and functions. An exception is using common abbreviations (such as 'VTK') specific to the domain of knowledge. 




Indentation and tabs
================================

- No character in the entire source code must be a tab.
- Generally, indentation of blocks uses four spaces though exceptions may apply.









Assertions and checks
======================

Make extensive use of the assert macro. The assertions must have no side effects except possibly output.

At the beginning of every function, call the `check` method of that class except for the `check` method itself and functions that must avoid this for technical reasons. 

At the beginning of the function call, all object arguments need to be checked.

At the end of a method, all non-const object arguments need to be checked, and the `this` pointer of any non-const method as well. Moreover, any return value must be checked. 

There are circumstances when a check is not performed, due to implementation-specific situations. 
In that case, instead, there should be a comment like:

```
// No check possible
```

If there are non-semantic reasons, such as performance, to disable the check,
then comment out the check and add an explanation:

```
// check(); // NOTE: performance
```












 
Enum Policy
===========
 
All enumerations must be `enum class` as opposed to the usual enums.
The program is possibly compiled with `-fstrict-enum`, which effectively prohibits any fancy arithmetics.
All enumerations must have a base type, typically `unsigned char`, large enough to represent all values.

Contrary to the usual C conventions, the enumeration values have simpler names, there is no need to include
the name of the enumeration type, because the values are only accessible through the enumeration identifier
in the first place.
 



Exception policy
================

The entire code adheres to a fail-fast philosophy: any error should lead to an immediate program termination. 
Typically, this is done with program termination once any assertion is triggered. 

For reasons of performance and simplicity, the program is generally compiled with the flag `-fno-exceptions`, disabling all exception handling.
Nevertheless, some methods should always be marked as `noexcept` to accommodate style conventions and silence linters. 
The following methods are always declared `noexcept`:

- move constructor
- move assignment operator
- destructor

Should it be desired to throw exceptions, then is done from within the custom assertion macro. 
A feature test reads:

```cpp
// gcc and clang define __cpp_exceptions = 199711 if exceptions are on
// gcc and clang also define __EXCEPTIONS
#if __cpp_exceptions == 199711
    // with exceptions...
#else
    // without exceptions...
#endif
// see also: https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_exceptions.html
```


Virtual destructor policy
=========================

The destructors of all classes shall be virtual unless the class meets the following conditions:

- The class is declared final
- The class has no parent class

As a rationale, we recall that virtual destructors are mandatory if the class is part of a inheritance hierarchy,
so the conditions describe exactly when a destructor is permitted to be non-virtual. 
A non-virtual destructor will generally have slightly better performance.


C++ Attributes 
=========================

The following attributes are defined via macros and may be used within the code:

- HOTCALL: signifies that a function will be called repeatedly and should be inlined 
- PACKED: marks data structures that should be packed, without regard to alignment 
- LIKELY/UNLIKELY: indicates code paths that are LIKELY/UNLIKELY to be followed
- UNUSED: marks variables and parameters that are not used


`final` keyword
=========================

Classes that are not part of a class hierarchy should be declared `final`.
Rationale: the classes are currently only used within the code,
and any class being part of a hierarchy, should be a conscious decision.
This avoids the vtable in some cases.


Rule of Five
=========================

We adhere to the following rule of Five. If a class features one of the following methods,

- Copy constructor
- Move constructor
- Copy assignment operator
- Move assignment operator
- Destructor

then it must provide all of these five methods. 
This can be via implementation or via `default` or `delete`.
It is preferable to use the `default` implementation whenever possible. 
These methods can be `delete`d until the need for them arises. 


Explicit keyword
=========================

All constructors with arguments, except for the copy and move constructors, are declared `explicit`.
This is irrespective of the number of arguments that these constructors take.

**Rationale:** Marking almost all constructors explicit makes the overload resolution more robust to code change, e.g., when an argument is made optional or an accidental assignment is made. 



General guidelines on resources managed by objects 
==================================================

We want to provide general guidelines for the copy and move operations and destructors. 
The objects in this library manage resources that can be categorized very broadly into the following:

- plain data 
- STL classes
- raw pointers to single objects and arrays

While the handling of plain data and STL classes is canonical, we discuss some policies for raw pointers:

- raw pointers are used to improve performance, even when optimization is disabled 
- after the constructor, raw pointers are non-null
- raw pointers are nullptr if the object is in the moved-from state 
- operations with nullpointers (except object destruction) may lead to program termination at any time 

This last item contrasts with the STL. It allows us to improve performance while assuming that moves happen at the end of an object's practical lifetime.

These are general guidelines that suit our particular situation. 

## Copy constructor guidelines

- The copy constructor will perform copy construction on all plain data members and members with copy constructor, except for pointers to resources
- pointers will be initialized in accordance with the properties of the source
- They will be initialized as nullptr if the above is too complicated or the original pointer was nulltpr
- Any copy operations take place in the body 
- The copy constructor finishes with a check and with `return *this`

## Move constructor guidelines

- The move constructor will perform move construction on all plain data members and members with copy constructor, including pointers to resources
- source pointers are set to zero
- The move constructor finishes with a check and with `return *this`

## Copy assignment operator guidelines 

- The operator finishes with a check and with `return *this`
- The operator checks for self-assignment and performs anything else only then
- Plain data members and members with copy assignment operator experience copy assignment
- Pointers to resources are checked for nullptr; then the data are copied manually

## Move assignment operator guidelines 

- The operator finishes with a check and with `return *this`
- The operator checks for self-assignment and performs anything else only then
- Plain data members and members with move assignment operator move copy assignment
- Pointers to resources are checked for nulltpr; then the pointer is moved to the destination and the source gets a nullptr

## Destructor guideline

- If the object is not a moved-from state, then optionally check at the beginning 
- If pointers to resources are not nullptr, then delete them



Prefixes for code readability
==================================================

Using prefixes such as `get`, `set`, `find`, `count`, `compute`, and `calculate` should improve the readability of the code. Possible convention:

* `is` or `has` for checking boolean conditions: **Example**: `bool is_symmetric()`

* `get` used to retrieve or access a value or object, typically not modifying the state of an object and returning a value or reference (const or non-const) to a member variable. **Example**: `double get_data(int i)`

* `set` used to assign or update a value or object, typically taking one or more parameters and update the state of an object without returning a value. **Example**: `void set_value(int i, double d)`

* `compute` for functions that perform a computation or processing to produce and return a new value, often without side effects. **Example**: `compute_determinant()`

* `calculate` is similar to `compute`, mainly used for performing mathematical or logical calculations and returning the result. **Example**: `calculate_average()`

* `generate` is for functions that involve a lambda, typically for procedural generation of data.

* `create` or `make` for functions that instantiate objects or values and return them.

* `update` should be for functions that modify the internal state of an object, possibly with side effects.

* `reset` or `clear` for resetting an object.

* `add` or `insert` and `remove` or `delete`** for adding or removing objects from a collection.

* `count` counts the number of elements that satisfy a particular condition or the total number of elements in a collection. Returns an integer representing the count. **Example**: `countMembers()`

* `find` for searching for a particular element or condition within a collection or range and returning its position or a reference/pointer to it. It might return a special value or indicator (like end iterator, null pointer) if the search is unsuccessful. **Example**: `findValue( std::vector<int>, int )`

More prefixes: prefixes: 

- eval
- calc
- comp
- generate
- derive
- check
- complete
- count
- give
- take
- call
- load
- save
- print



NOTES Scott Meyers book
==================================================

-  7: declare destructors virtual in polymorphic base classes OK
-  9: never call virtual functions during construction or destruction OK
- 10: have assignments return reference to *this OK
- 11: handle self-assignment in operator=	 OK

- The prefix and postfix forms of the increment/decrement operators have different semantics
- No overloading of the operators &&, ||, or ,.
- If a class dynamically allocates memory, then define a copy constructor and an assignment operator.
- Use initialization before constructor body over within constructor body.
- Assignment operators return reference to *this.
- Do not return a reference when the return type is an object

Better make operators non-member or hidden friends
==================================================

- https://www.reddit.com/r/cpp_questions/comments/1cjgvfe/overloading_operators_from_outside_a_class_as/
- some mandatory member operators 
- other operators should be non-members
- hidden friend idiom: affects argument-dependent lookup because member friends are not in the global namespace
- DRY Principle: Implement non-mutating operators (e.g., +) in terms of mutating members (e.g., +=).
Tell me more about the custom comparison / hashing and the member functions?



