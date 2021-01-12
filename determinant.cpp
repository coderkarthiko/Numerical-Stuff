#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> vd;
typedef vector<vd> mat;

double det(mat M) {
  double det = 1;
  for(int i = 0; i < M.size(); i ++) {
    int r = i; 
    for(int j = i + 1; j < M.size(); j ++) {
      if(abs(M[j][i]) > abs(M[r][i])) {
        r = j;
      }
    }
    if(M[i][i] == 0) {
      return 0;
    }
    if(i != r) {
      for(int j = i; j < M.size(); j ++) {
        swap(M[i][j], M[r][j]);
      }
      det *= -1.;
    }
    det *= M[i][i];
    for(int j = i + 1; j < M.size(); j ++) {
      for(int k = i + 1; k < M.size(); k ++) {
        M[j][k] -= M[i][k] * (M[j][i] / M[i][i]);
      }
    }
  }
  return det;
}

mat get() {
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
