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
    prefixes: 
        eval, get, set, calc, compute, generate, derive, check, complete, count, 





Using prefixes like `get`, `set`, `find`, `count`, `compute`, and `calculate` in your function and method names can significantly improve the readability and understandability of your code by clearly indicating the intended action or purpose of each function. Here's a breakdown of what functions and methods with those prefixes should do, along with suggestions for additional relevant prefixes:

### 1. `get`
- **Purpose**: Used to retrieve or access a value or object. These functions typically do not modify the state of an object and return a value or reference (const or non-const) to a member variable.
- **Example**: `getColor()` returns the color property of an object.

### 2. `set`
- **Purpose**: Used to assign or update a value or object. These methods typically take one or more parameters and update the state of an object without returning a value.
- **Example**: `setColor(Color newColor)` assigns a new color to an object.

### 3. `find`
- **Purpose**: Used for searching for a particular element or condition within a collection or range and returning its position or a reference/pointer to it. It might return a special value or indicator (like end iterator, null pointer) if the search is unsuccessful.
- **Example**: `findUserById(int id)` searches for a user by their ID and returns a pointer or reference to the user object.

### 4. `count`
- **Purpose**: Used to count the number of elements that satisfy a particular condition or the total number of elements in a collection. Returns an integer representing the count.
- **Example**: `countActiveUsers()` returns the number of users marked as active.

### 5. `compute`
- **Purpose**: Used for functions that perform a computation or processing to produce and return a new value, often without side effects.
- **Example**: `computeTotalCost()` calculates and returns the total cost based on internal factors or parameters.

### 6. `calculate`
- **Purpose**: Similar to `compute`, it is used for performing mathematical or logical calculations and returning the result.
- **Example**: `calculateInterestRate(float principal, float rate, int time)` calculates the interest over a given period.

### Additional Relevant Prefixes

- **`is` or `has`**: Used for functions that return a boolean value based on a condition or state. For example, `isEmpty()` checks if a container is empty, `hasPermission()` checks if the current user has a specific permission.
- **`create` or `make`**: Used for functions that instantiate objects or values and return them. For example, `createUser()` might instantiate and return a new user object.
- **`update`**: Indicates that the function modifies the state of an object or system based on given parameters, often with side effects. For example, `updateRecord(Record record)` might update a record in a database.
- **`remove` or `delete`**: Used for functions that remove elements or data. For example, `removeUser(int userId)` might remove a user from a collection.
- **`add` or `insert`**: Indicates adding or inserting elements into a collection. For example, `addProduct(Product product)` would add a product to a catalog.
- **`reset` or `clear`**: Used for functions that reset an object's state or clear contents of a collection. For example, `resetSettings()` resets settings to default values.

Each prefix suggests a specific action or outcome, making it easier for developers to understand what a function does at a glance. Consistently applying these naming conventions can make your code more intuitive and maintainable.