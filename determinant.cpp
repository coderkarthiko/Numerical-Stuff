// determinant using Gaussian elimination - O(N^3)

#include <iostream>
#include <vector>

using namespace std;

typedef vector<vector<int>> mat;

double det(mat M) {
  double det = 1;
  for(int i = 0; i < M.size(); i ++) {
    det *= M[i][i];
    if(det == 0) {
      return 0;
    }
    for(int j = i + 1; j < M.size(); j ++) {
      double f = M[j][i] / M[i][i];
      for(int k = i; k < M.size(); k ++) {
        M[j][k] -= M[i][k] * f; 
      }
    }
  }
  return det;
}

int main() {
  mat M = {{5, 2, 9, 4}, 
           {1, 5, 0, 4}, 
           {7, 8, 1, 3},
           {3, 3, 3, 1}};
  cout << det(M);
}
