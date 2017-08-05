
    
    
    o Rewrite core float vector class 
      
      Write it up in a manner that is close to the STL vector class.
      Perhaps even make it a descendant of std::vector<Float>
      and wrap it only thinly.
      https://stackoverflow.com/questions/2034916/is-it-okay-to-inherit-implementation-from-stl-containers-rather-than-delegate
      There seem to be complications, so it should be delayed until further notice.
      There is rather a speed-up if we replace it by generic C++ memory allocation.
      In particular, it does not really mesh with later efforts of parallelization. 
      Furthermore, it is better to entirely hide the implementation from the user.
      
    
    o Parallelization of Float Vector class
      
      Many of the methods in the float vector class 
      can be parallelized with openMP.
      - Constructors
      - zero, scale
      - NOT random
      - scalarproductwith
      - norm, maxnorm, lpnorm
      - add vectors 
      
      
    o Implement vector slices 
      
      A vector slice refers to a part of a vector.
      The slice knows the original vector and 
      some data determine how to access the original members.
      
      Best approach would be to introduce an abstract class
      for vectors that captures the interface. 
      Then fork off the original class of vectors 
      and the new slice implementation. 
      SEE ALSO Implement lambda-based vectors
      
      
    o Implement lambda-based vectors 
      
      The get/set methods can then be given in terms 
      of lambdas that produce the required terms/references 
      on the spot. This gives the most general functionality.
      
      Note that read-only vectors can be implemented 
      by having the set operation cause an error.
      
    
    
    
    
    
    
    
    
    
    
