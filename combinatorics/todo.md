

# MODULE Combinatorics

This module provides utilities for the handling of mappings between integer ranges. 
This includes multiindices and alternator indices.

The class IndexRange merely abstracts the notion of intervals of integers. 
Any IndexRange is defined uniquely by the integers that are officially contained in it.

IndexMaps abstract mappings from one interval of integers into another one. 
Any mapping is uniquely defined by its domain and its codomain.

Multiindices are just IndexMaps with enhanced functionality.
It is understood that Multiindices map into the non-negative integers.

Finally, this module also provides utility functions
that generate specific (ordered) sets of indexmaps (and multiindices).
    
    
## TODO: equip print method with verbosity signal



The print method should receive an additional (optional) parameter 
to print pure data or embellish the output with some text. 
    
    Indexrange
    * check
    * check im constructor 
    * unittest komplett 
    
    IndexMap 
    * check 
    * check im constructor 
    * unittest komplett 
    
    MultiIndex
    * check 
    * check im constructor 
    * unittest komplett 
    
    generateMultiIndices 
    * robusten
    * unittest komplett 
    
    GenerateIndexMaps
    * robusten
    * unittest komplett 
    
    ---------------------------
    
    IndexMap
    - interne daten besser konzentrieren auf wenige function calls 
    - C++ vector ersetzen durch C style array 
    
