#include <iostream>
#include <vector>
#include <stdexcept>

using namespace std;

typedef vector<double> vd;
typedef vector<vd> mat;

mat inv(mat M) {
  mat I(M.size(), vd(M.size(), 0));
  for(int i = 0; i < M.size(); i ++) {
    I[i][i] = 1;
  }
  for(int i = 0; i < M.size(); i ++) {
    double f = M[i][i];
    if(f == 0) {
      throw invalid_argument("matrix is singular");
    }
    for(int j = 0; j < M.size(); j ++) {
      M[i][j] /= f;
      I[i][j] /= f;
    }
    for(int j = 0; j < M.size(); j ++) {
      if(j != i) {
        f = M[j][i];        
        for(int k = 0; k < M.size(); k ++) {
          M[j][k] -= f * M[i][k];
          I[j][k] -= f * I[i][k];
        }
      }
    }
  }
  return I;
}

mat get() {
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

void print_mat(mat M) {
  for(auto row: M) {
    for(auto col: row) {
      cout << col << "  ";
    }
    cout << "\n\n";
  }
}

int main() {
  print_mat(inv(get()));
}
