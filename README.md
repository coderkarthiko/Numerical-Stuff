## Numerical Stuff

C++ implementations of various numerical algorithms (mostly from Linear Algebra). 

## determinant.cpp
Finds determinant of an N * N matrix "A" in O(N^3) using Gaussian elimination. det(A)=tr(U) where U is upper triangular form of A.

## inverse.cpp
Finds inverse of an N * N matrix using Gauss-Jordan elimination.

## linear_solver.cpp
Solves Ax = y by computing (A^-1)(yT). Basic stuff.

## simplex.cpp
Finds value of x that maximizes/minimizes (c^T)x subject to Ax <= b and x >= 0 where c is a vector of dimension N, A is a P * N matrix and b is a vector of dimension P. 
