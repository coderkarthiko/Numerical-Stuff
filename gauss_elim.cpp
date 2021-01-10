// Gaussian Elimination - 10-01-2021

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
    cout << x << "\n\n";
  }
}
int main() {
  solve({{1.00, 0.00, 0.00, 0.00, 0.00, 0.00},
         {1.00, 0.63, 0.39, 0.25, 0.16, 0.10},
         {1.00, 1.26, 1.58, 1.98, 2.49, 3.13},
         {1.00, 1.88, 3.55, 6.70, 12.62, 23.80},
         {1.00, 2.51, 6.32, 15.88, 39.90, 100.28},
         {1.00, 3.14, 9.87, 31.01, 97.41, 306.02}}, 
         {-0.01, 0.61, 0.91, 0.99, 0.60, 0.02});
}
