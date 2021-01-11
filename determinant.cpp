#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

typedef vector<vector<int>> matrix;

matrix MINOR(matrix M, int r, int c) {
  matrix minor;
  for(int i = 0; i < M.size(); i ++) {
    if(i != r) {
      minor.push_back({});
      for(int j = 0; j < M.size(); j ++) {
        if(j != c) {
          minor.back().push_back(M[i][j]);
        }
      }
    }
  }
  return minor;
}

int DET(matrix M) {
  if(M.size() == 1) {
    return M[0][0];
  }
  else {
    int det = 0;
    for(int i = 0; i < M.size(); i ++) {
      det += pow(-1, i) * M[0][i] * DET(MINOR(M, 0, i));
    }
    return det;
  }
}

int main() {
  matrix mat = {{5, 2, 9, 4}, 
                {1, 5, 0, 4}, 
                {7, 8, 1, 3},
                {3, 3, 3, 1}};
  cout << DET(mat);
}