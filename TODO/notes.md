## Exception policy

The entire code adheres to a fail-fast philosophy: any error should lead to an immediate program termination. 
For that reason, the compilation is generally compiled with the flag `-fno-exceptions`, disabling all exception handling.
Nevertheless, some methods should always be marked as `noexcept` to accommodate style conventions and silence linters. 
The following methods are always declared `noexcept`:

- move constructor
- move assignment operator
- destructor

## Virtual destructor policy

The destructors of all classes shall be virtual unless the class meets the following conditions:

- The class is declared final
- The class has no parent class

As a rationale, we recall that virtual destructors are mandatory if the class is part of a inheritance hierarchy,
so the conditions describe exactly when a destructor is permitted to be non-virtual. 
A non-virtual destructor will generally have slightly better performance.

## Final keyword

Classes that are not part of a class hierarchy should be declared `final`.
Rationale: the classes are currently only used within the code,
and any class being part of a hierarchy, should be a conscious decision.
This avoids the vtable in some cases.


## Rule of Five

We adhere to the following rule of Five. If a class features one of the following methods

- Copy constructor
- Move constructor
- Copy assignment operator
- Move assignment operator
- Destructor

then it must provide all of these five methods. 
This can be via implementation or via `default` or `delete`.


# General guidelines on resources managed by objects 

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

# Copy constructor guidelines

- The copy constructor will perform copy construction on all plain data members and members with copy constructor, except for pointers to resources
- pointers will be initialized in accordance with the properties of the source
- They will be initialized as nullptr if the above is too complicated or the original pointer was nulltpr
- Any copy operations take place in the body 
- The copy constructor finishes with a check and with `return *this`

# Move constructor guidelines

- The move constructor will perform move construction on all plain data members and members with copy constructor, including pointers to resources
- source pointers are set to zero
- The move constructor finishes with a check and with `return *this`

# Copy assignment operator guidelines 

- The operator finishes with a check and with `return *this`
- The operator checks for self-assignment and performs anything else only then
- Plain data members and members with copy assignment operator experience copy assignment
- Pointers to resources are checked for nullptr; then the data are copied manually

# Move assignment operator guidelines 

- The operator finishes with a check and with `return *this`
- The operator checks for self-assignment and performs anything else only then
- Plain data members and members with move assignment operator move copy assignment
- Pointers to resources are checked for nulltpr; then the pointer is moved to the destination and the source gets a nullptr

# Destructor guideline

- If the object is not a moved-from state, then optionally check at the beginning 
- If pointers to resources are not nullptr, then delete them

# Explicit keyword

All constructors with arguments, except for the copy and move constructors, are declared `explicit`.
This is irrespective of the number of arguments that these constructors take.


## Prefixes for code readability

Using prefixes such as `get`, `set`, `find`, `count`, `compute`, and `calculate` should improve the readability of the code. Possible convention:

1. `is` or `has` for checking boolean conditions: **Example**: `bool is_symmetric()`

2. `get` used to retrieve or access a value or object, typically not modifying the state of an object and returning a value or reference (const or non-const) to a member variable. **Example**: `double get_data(int i)`

3. `set` used to assign or update a value or object, typically taking one or more parameters and update the state of an object without returning a value. **Example**: `void set_value(int i, double d)`

4. `compute` for functions that perform a computation or processing to produce and return a new value, often without side effects. **Example**: `compute_determinant()` 

5. `calculate` is similar to `compute`, mainly used for performing mathematical or logical calculations and returning the result. **Example**: `calculate_average()` 

6. `create` or `make` for functions that instantiate objects or values and return them. 

7. `update` should be for functions that modify the internal state of an object, possibly with side effects

8. `reset` or `clear` for resetting an object.

9. `add` or `insert` and `remove` or `delete`** for adding or removing objects from a collection.

10. `count` counts the number of elements that satisfy a particular condition or the total number of elements in a collection. Returns an integer representing the count. **Example**: `countMembers()` 

11. `find` for searching for a particular element or condition within a collection or range and returning its position or a reference/pointer to it. It might return a special value or indicator (like end iterator, null pointer) if the search is unsuccessful. **Example**: `findValue( std::vector<int>, int )` 




## NOTES Scott Meyers book

- 7: declare destructors virtual in polymorphic base classes OK
- 9: never call virtual functions during construction or destruction OK
- 10: have assignments return reference to *this OK
- 11: handle self-assignment in operator=	 OK


### Better make operators non-member

- https://www.reddit.com/r/cpp_questions/comments/1cjgvfe/overloading_operators_from_outside_a_class_as/
- some mandatory member operators 
- other operators should be non-members
- hidden friend idiom: affects argument-dependent lookup because member friends are not in the global namespace
- DRY Principle: Implement non-mutating operators (e.g., +) in terms of mutating members (e.g., +=).
Tell me more about the custom comparison / hashing and the member functions?



2D pullback:

```
k=1	max		| PF |             = max                          = max
k=2	max * min       | PF |  max        = max * min / min              = max

3D pullback:
k=1	max 		| PF |             = max 			  = max 
k=2	max * mid       | PF |  max 	   = max * mid / min 		  = max/min * mid 
k=3	max * mid * min | PF |  max * mid  = max * mid * min / min * mid  = max 
```

--- PRODUCT NAME:

    Coffeec.org
    Coffeecpp.org
    FEECpp.org
    FEEC++

  http://asp-software.org/www/misv_resources/business-articles/how-to-name-an-app-a-program-a-company-or-a-service/

--- SOFTWARE:

 * Provide software for rapid prototyping of tensor-valued finite element methods over simplicial grids.
 * Gain practical experience in the implementation of FEEC 
 
 * YES: Maintainability, Flexibility, Produce qualitative results
 * NO: Performance, industry

 * Future features:
 * * mesh refinement in several dimensions (10h)
 * * better preconditioners (4h)
 * * parallelism (lots and lots)
 
    -----------------------
    prefixes: eval, get, set, calc, compute, generate, derive, check, complete, count, 




