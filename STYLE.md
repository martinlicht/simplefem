
STYLE GUIDE
===========


Naming Conventions
------------------

Use `CapitalizedCamelCase` for classes, structs, enums and namespaces.

Use `CapitalizedCamelCase`, `camelCase` or `snake_case` for functions.

Use `snake_case` for variables.

Underscores also be used in the names of classes, structs, enums, and functions, if they indicate variations of a common principle, such as different implementations of the same idea.

As by common and nearly universal convention, macros must always be in `UPPERCASE_WITH_UNDERSCORES`.

In general, do not abbreviate when naming classes and functions. An exception is using common abbreviations (such as 'VTK') specific to the domain of knowledge. 



No Tabs
-------

No character in the entire source code must be a tab.





Assertions and checks
---------------------

Make extensive use of the assert macro. 
All asserts must be `const` to their arguments, no side effects should occur.

At the beginning of every function,
call the check() method of that class except for
- the check method itself
- methods that avoid check for good reason, see below.
Likewise, all object arguments need to be checked at the beginning of the function.

At the end of a method, all non-const object arguments need to be checked,
and the this pointer of any non-const method as well.


There are circumstances when a check is not performed.


If a check is not supposed to be called for semantic reasons,
instead of check there should be a comment like 

```
// No check possible
```

If there are non-semantic reasons, such as performance, to disable the check,
then comment out the check and add an explanation:

```
// check(); // NOTE: performance
```





Indentation
-----------

Generally, indentation of blocks uses four spaces though exceptions may apply.

