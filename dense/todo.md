

# Dense Matrix class 

## TODO List 

* check
* check im constructor 
- unittest komplett 

## Future features 
- interne daten umstellen auf C-style arrays 
- high speed product 
- einige algorithmen intern spezialisieren:
- determinante 
- cofactor matrix / subdeterminant
- inverse
- iterative matrix inversion: 
- newton-based methods 
- fixpoint iteration (Richardson) with preconditioner 
- others?


## Dense Matrix algorithms 

matrix tensor product 
* unit test 

Matrix functions: det, inv, transpose 
* unit test : transpose
* unit test : inverse and determinant 
- unit test : subdeterminant matrix 

Scalar functions
* unit test
- infty as special value 

Simplesolver
- unit test

LU Factorization
- algorithmen 
- unit test 

Cholesky 
* unit test
- in place inversion 
- different access patterns 

QR Factorization (Gram-Schmidt)
- unit test




## Dense Matrix algorithms to implement:

- Gerschgorin (row/column) circles, and maximal estimate 
perhaps more estimators 

- Diagonal solve 
- Left triangular solve 
- Right triangular solve 
- Unit Left triangular solve 
- Unit Right triangular solve 
- (Averages between left and right triangular solves)

- check whether diagonal 
- check whether (unit) left/right triangular 

- LU decomposition 
- LU decomposition, row pivot
- LU decomposition, column pivot 
- LU decomposition, full pivot 

- Cholesky 
- Cholesky, pivot 

- stabilized Gram-Schmidt 

- QR decomposition 


## MatrixAlgorithms
- intern niederdimensionale spezialisierungen verwenden  

## MatrixAlgorithms
* check
* check im constructor 
- unittest komplett 












