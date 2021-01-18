#include "types.h"

bool strassen = true;

// print vector
void print_vd(vd v) {
    for (auto vx : v) {
        if (vx > threshold) {
            cout << vx << " ";
        }
        else {
            cout << 0 << " ";
        }
    }
}

// print matrix
void print_mat(mat M) {
    for (auto r : M) {
        for (auto c : r) {
            if (c > threshold) {
                cout << c << " ";
            }
            else {
                cout << 0 << " ";
            }
        }
        cout << "\n";
    }
}

// return identity matrix of dim N
mat eye(int n) {
    mat I(n, vd(n));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1.;
    }
    return I;
}

// set a block of a matrix equal to another matrix
void eq(mat& A, mat B, pii rc) {
    for (int i = 0; i < B.size(); i++) {
        for (int j = 0; j < B.front().size(); j++) {
            A[rc.first + i][rc.second + j] = B[i][j];
        }
    }
}

// return a block of a matrix
mat block(mat M, pii r, pii c) {
    mat B(r.second - r.first, vd(c.second - c.first));
    for (int i = 0; i < r.second - r.first; i++) {
        for (int j = 0; j < c.second - c.first; j++) {
            B[i][j] = M[r.first + i][c.first + j];
        }
    }
    return B;
}

// add two matrices
mat add(mat A, mat B) {
    assert(A.size() == B.size());
    assert(A.front().size() == B.front().size());
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.front().size(); j++) {
            A[i][j] += B[i][j];
        }
    }
    return A;
}

// subtract two matrices
mat subtract(mat A, mat B) {
    assert(A.size() == B.size());
    assert(A.front().size() == B.front().size());
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.front().size(); j++) {
            A[i][j] -= B[i][j];
        }
    }
    return A;
}

// matrix transpose
mat transpose(mat A) { // return transpose of matrix 
    mat AT(A[0].size(), vd(A.size()));
    for (int i = 0; i < AT.size(); i++) {
        for (int j = 0; j < AT[0].size(); j++) {
            AT[i][j] = A[j][i];
        }
    }
    return AT;
}

// mat mul
mat mul(mat A, mat B) {
    mat C(A.size(), vd(B.front().size()));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < B.front().size(); j++) {
            for (int k = 0; k < A.size(); k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// determinant of a matrix
ld det(mat M) {
    ld det = 1; // set initial det val = 1
    for (int i = 0; i < M.size(); i++) {
        int r = i; // partial pivoting to prevent decimal errors
        for (int j = i + 1; j < M.size(); j++) {
            if (abs(M[j][i]) > abs(M[r][i])) {
                r = j;
            }
        }
        if (i != r) {
            // if i =/= r, swap row elements except pivot element and multiply det by -1
            for (int j = i; j < M.size(); j++) {
                swap(M[i][j], M[r][j]);
            }
            det *= -1.;
        }
        if (abs(M[i][i]) <= threshold) { // if pivot element is 0, det(M) = 0
            return 0;
        }
        det *= M[i][i]; // multiply det by pivot element
        for (int j = i + 1; j < M.size(); j++) { // row elimination
            for (int k = i + 1; k < M.size(); k++) {
                M[j][k] -= M[i][k] * (M[j][i] / M[i][i]);
            }
        }
    }
    return det; // return determinant
}

// inverse of a matrix
mat inv(mat M) {
    // make identity matrix of same dim as M
    mat I = eye(M.size());
    for (int i = 0; i < M.size(); i++) {
        I[i][i] = 1.;
    }
    for (int i = 0; i < M.size(); i++) {
        int r = i; // partial pivoting 
        for (int j = i + 1; j < M.size(); j++) {
            if (abs(M[j][i]) > abs(M[r][i])) {
                r = j;
            }
        }
        swap(M[i], M[r]);
        swap(I[i], I[r]);
        double f = M[i][i];
        assert(abs(f) >= threshold);
        for (int j = 0; j < M.size(); j++) { // row / pivot
            M[i][j] /= f;
            I[i][j] /= f;
        }
        for (int j = 0; j < M.size(); j++) { // row elimination 
            if (j != i) {
                f = M[j][i];
                for (int k = 0; k < M.size(); k++) {
                    M[j][k] -= f * M[i][k];
                    I[j][k] -= f * I[i][k];
                }
            }
        }
    }
    return I; // return inverse
}

// input a vector
vd get_vd() {
    int n;
    cin >> n;
    vd v(n);
    for (int i = 0; i < n; i++) {
        std::cin >> v[i];
    }
    return v;
}

// input a matrix
mat get_mat() {
    int n, m;
    cin >> n >> m;
    mat M(n, vd(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cin >> M[i][j];
        }
    }
    return M;
}
