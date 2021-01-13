// C++ program to find determinant of a matrix 
// Convert matrix M into upper triangular form using Gaussian elimination in O(N^3) time
// Then, det(M) is simply prod(trace(M))
// Laplace formula takes O(N!) time! Way slower...

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> vd;
typedef vector<vd> mat;

double det(mat M) {
  double det = 1; // set initial det val = 1
  for(int i = 0; i < M.size(); i ++) {
    int r = i; // partial pivoting to prevent decimal errors
    for(int j = i + 1; j < M.size(); j ++) {
      if(abs(M[j][i]) > abs(M[r][i])) {
        r = j;
      }
    }
    if(M[i][i] == 0) { // if pivot element is 0, det(M) = 0
      return 0; 
    }
    if(i != r) { 
      // if i =/= r, swap row elements except pivot element and multiply det by -1
      for(int j = i; j < M.size(); j ++) {
        swap(M[i][j], M[r][j]);
      }
      det *= -1.;
    }
    det *= M[i][i]; // multiply det by pivot element
    for(int j = i + 1; j < M.size(); j ++) { // row elimination
      for(int k = i + 1; k < M.size(); k ++) {
        M[j][k] -= M[i][k] * (M[j][i] / M[i][i]);
      }
    }
  }
  return det; // return determinant
}

mat get() { // input matrix
  int n;
  cin >> n;
  mat M(n, vd(n));
  for(int i = 0; i < M.size(); i ++) {
    for(int j = 0; j < M.size(); j ++) {
      cin >> M[i][j];
    }
  }
  return M;
}

int main() {
  cout << fixed << det(get());
}
