#include <iostream>
#include <vector>

using namespace std;

typedef vector<vector<double>> mat;
typedef vector<double> vd;

void solve(mat A, vd y) {
  for(int i = 0; i < A.size(); i ++) {
    for(int j = i; j < A.size(); j ++) {
      double f = A[j][i]; 
      // if A[i][i] is 0, then print A and y
      if(f == 0) {
        for(int i = 0; i < A.size(); i ++) {
          for(int j = 0; j < A.size(); j ++) {
            cout << A[i][j] << "  ";
          }
          cout << "    " << y[i] << "\n\n";
        }
        return;
      }
      // if A[i][i] not zero, divide row by first element
      for(int k = i; k < A.size(); k ++) {
        A[j][k] /= f;
      }
      y[j] /= f;
    }
    // subtract rows till echelon form achieved
    for(int j = i + 1; j < A.size(); j ++) {
      for(int k = i; k < A.size(); k ++) {
        A[j][k] -= A[i][k];
      }
      y[j] -= y[i];
    }
  }
  // back substitution
  for(int i = A.size() - 1; i > -1; i --) {
    for(int j = A.size() - 1; j > i; j --) {
      y[i] -= A[i][j] * y[j];
    }
  }
  // print answer
  for(auto x: y) {
    cout << x << " ";
  }
}

int main() {
  solve({{1., 8., 2.}, 
         {14., 12., 21.},
         {31., 3., 6.}}, {15., 141., 100.});
}
