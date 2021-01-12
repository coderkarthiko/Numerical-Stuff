#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> vd;
typedef vector<vd> mat;

// function
mat inv(mat A) {
  // update det as elimination happens
  double det = 1.; 
  // initialize identity matrix
  mat I(A.size(), vd(A.size(), 0));
  for(int i = 0; i < A.size(); i ++) {
    I[i][i] = 1;
  }
  for(int i = 0; i < A.size(); i ++) {
    // pivot element
    double f = A[i][i];
    // if A[i][i] is 0, determinant is 0 => no inverse
    if(f == 0.) {
      cout << "NO INVERSE, |A| = 0";
    }
    // divide row by pivot element
    for(int j = 0; j < A.size(); j ++) {
      A[i][j] /= f;
      I[i][j] /= f;
    }
    // update det
    det *= f;
    // row operations till RREF
    for(int j = 0; j < A.size(); j ++) {
      if(j != i) {
        f = A[j][i];        
        for(int k = 0; k < A.size(); k ++) {
          A[j][k] -= f * A[i][k];
          I[j][k] -= f * I[i][k];
        }
      }
    }
  }
  // print determinant and return inverse
  cout << "INVERSE EXISTS, |A| = " << det << "\n\n";
  return I;
}


int main() {
  mat A = {{1.3, 5., 11., 1., 11., 7.},
           {6., 4., 2.1, 14, 8.},
           {1.7, 8., 9., 3.2, 0.},
           {7., 1.5, 3., 2., 1.},
           {5.6, 1.6, 5., 4.2, 4.}};
  for(auto r: inv(A)) {
    for(auto c: r) {
      cout << c << "  ";
    }
    cout << "\n\n";
  }
}