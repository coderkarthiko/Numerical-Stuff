#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

using namespace std;

typedef vector<double> vd;
typedef vector<vd> mat;

mat inv(mat M) {
  // make identity matrix of same dim as M
  mat I(M.size(), vd(M.size(), 0));
  for(int i = 0; i < M.size(); i ++) {
    I[i][i] = 1;
  }
  for(int i = 0; i < M.size(); i ++) {
    // partial pivoting to prevent decimal errors
    int r = i;
    for(int j = i + 1; j < M.size(); j ++) { 
      if(abs(M[j][i]) > abs(M[r][i])) {
        r = j;
      }
    }
    swap(M[i], M[r]); // swap row M[i] with M[r]
    swap(I[i], I[r]); // swap row I[i] with I[r]
    double f = M[i][i]; // store pivot element in f
    if(f == 0) {
      // if pivot element becomes 0, inverse doesn't exist because det(M) = 0
      throw invalid_argument("NO INVERSE"); 
    }
    for(int j = 0; j < M.size(); j ++) {
      M[i][j] /= f; // divide row M[i] by f 
      I[i][j] /= f; // divide row I[i] by f
    }
    for(int j = 0; j < M.size(); j ++) { // row elimination 
      if(j != i) {
        f = M[j][i]; // store scale element
        for(int k = 0; k < M.size(); k ++) {
          M[j][k] -= f * M[i][k]; // subtract f * M[i] from M[j] 
          I[j][k] -= f * I[i][k]; // subtract f * M[i] from I[j]
        }
      }
    }
  }
  return I; // return inverse
}

mat mul(mat A, mat B) {
  mat AB(A.size(), vd(B[0].size()));
  for(int i = 0; i < AB.size(); i ++) {
    for(int j = 0; j < AB[0].size(); j ++) {
      for(int k = 0; k < AB.size(); k ++) {
        AB[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return AB;
}

mat transpose(mat A) { // return transpose of matrix A
  mat AT(A[0].size(), vd(A.size()));
  for(int i = 0; i < AT.size(); i ++) {
    for(int j = 0; j < AT[0].size(); j ++) {
      AT[i][j] = A[j][i];
    }
  }
  return AT;
}

vd solve(mat A, vd y) { // solve Ax = y where x = inv(A)(yT)
  vd x;
  for(auto r: mul(inv(A), transpose({y}))) {
    x.push_back(r[0]);
  }
  return x;
}

mat get_mat() { // input matrix
  int n;
  cin >> n;
  mat M(n, vd(n));
  for(int i = 0; i < n; i ++) {
    for(int j = 0; j < n; j ++) {
      cin >> M[i][j];
    }
  }
  return M;
}

vd get_vd() { // input vector
  int n;
  cin >> n;
  vd v(n);
  for(int i = 0; i < n; i ++) {
    cin >> v[i];
  }
  return v;
}

void print_mat(mat M) { // print matrix
  for(auto row: M) {
    for(auto col: row) {
      cout << col << "  ";
    }
    cout << "\n";
  }
}

int main() {
  print_mat(transpose({solve(get_mat(), get_vd())}));
}