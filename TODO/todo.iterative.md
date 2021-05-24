  
  
# (LOW) Iterative Methods to implement
  
The following iterative solvers can be implemented.
  - [x] Residual Minimizing Descent
  - [x] Conjugate Residual Method 
  - [x] Conjugate Residual Method on Normal Equations
  - [ ] Richardson iteration 
  - [ ] Gradient energy descent 
  - [ ] Gradient residual descent 
  - [ ] Symmetric Lanczos minimum residual method 

# (LOW) GMRES with Restart 

Implement the generalized minimal residual method
where the search directions are rebuilt from scratch
after a fixed number of iteration vectors have 
been constructed.

# (LOW) Rewrite algorithms to be complex number stable 
  
All algorithms should be written in a manner 
that is also correct when using complex numbers. 
This should be accompanied by a written exposition
of Krylov subspace methods.


  
# Preconditioners to implement 

  - [ ] Jacobi preconditioner 
  - [ ] different scaling preconditioners
  - [ ] Gauss-Seidel preconditioner
  - [x] SOR + SSOR preconditioner 
  - [ ] block diagonal preconditioner 
  - [ ] block gauss-seidel preconditioner 
  - [ ] adjustable gauss-seidel preconditioner 
  - [ ] Polynomial preconditioners 
  

# (LOW) Provide Preconditioned variants for all iterative methods
  
For each iterative method there should be a preconditioned 
method available. New iterative methods should only be added
if the preconditioned variant is added too.
  

  
    
    
  
  
  
  
  
  
