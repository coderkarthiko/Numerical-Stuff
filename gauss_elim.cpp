#include <iostream>
#include <vector>

using namespace std;

typedef vector<vector<double>> mat;
typedef vector<double> vd;

void solve(mat M, vd y) {
  for(int i = 0; i < M.size(); i ++) {
    for(int j = i; j < M.size(); j ++) {
      double f = M[j][i]; 
      if(f == 0) {
        for(int i = 0; i < M.size(); i ++) {
          for(int j = 0; j < M.size(); j ++) {
            cout << M[i][j] << "  ";
          }
          cout << "    " << y[i] << "\n\n";
        }
        return;
      }
      for(int k = i; k < M.size(); k ++) {
        M[j][k] /= f;
      }
      y[j] /= f;
    }
    for(int j = i + 1; j < M.size(); j ++) {
      for(int k = i; k < M.size(); k ++) {
        M[j][k] -= M[i][k];
      }
      y[j] -= y[i];
    }
  }
  for(int i = M.size() - 1; i > -1; i --) {
    for(int j = M.size() - 1; j > i; j --) {
      y[i] -= M[i][j] * y[j];
    }
  }
  for(auto x: y) {
    cout << x << " ";
  }
}

int main() {
  solve({{1., 8., 2.}, 
         {14., 12., 21.},
         {31., 3., 6.}}, {15., 141., 100.});
}