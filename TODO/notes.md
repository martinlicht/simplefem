2D pullback:
k=1	max		| PF |             = max                          = max 
k=2	max * min       | PF |  max        = max * min / min              = max

3D pullback:
k=1	max 		| PF |             = max 			  = max 
k=2	max * mid       | PF |  max 	   = max * mid / min 		  = max/min * mid 
k=3	max * mid * min | PF |  max * mid  = max * mid * min / min * mid  = max 


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




# Prefixes for code readability 

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


