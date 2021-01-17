#include "types.h"

// for matrix multiplication, standard mat mul algorithm - O(N^3)
// Strassen's algorithm is asymptotically faster - ~O(N^2.8704)
// but it has a large constant upfront and it accumulates decimal errors
// therefore, standard mat mul is better for practical purposes
// to use Strassen's algorithm, set variable "strassen" to true

bool strassen = false;

// print vector
void print_vd(vd v) {
    for (auto vx : v) {
        cout << vx << " ";
    }
}

// print matrix
void print_mat(mat M) {
    for (auto r : M) {
        for (auto c : r) {
            cout << c << " ";
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
void eq(mat &A, mat B, pii rc) {
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
    int mdim = max(A.size(), A.front().size()) * max(B.size(), B.front().size());
    if (!strassen || mdim <= 1024) { 
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
    else {
        int dim = pow(2, ceil(log(mdim) / log(4)));
        // Let R = PQ and P, Q and R are all of same dim
        mat P(dim, vd(dim, 0)), Q(dim, vd(dim, 0)), R(dim, vd(dim, 0));
        eq(P, A, { 0, 0 });
        eq(Q, B, { 0, 0 });
        mat P11 = block(P, { 0, dim / 2 }, { 0, dim / 2 });
        mat P12 = block(P, { 0, dim / 2 }, { dim / 2, dim });
        mat P21 = block(P, { dim / 2, dim }, { 0, dim / 2 });
        mat P22 = block(P, { dim / 2, dim }, { dim / 2, dim });
        mat Q11 = block(Q, { 0, dim / 2 }, { 0, dim / 2 });
        mat Q12 = block(Q, { 0, dim / 2 }, { dim / 2, dim });
        mat Q21 = block(Q, { dim / 2, dim }, { 0, dim / 2 });
        mat Q22 = block(Q, { dim / 2, dim }, { dim / 2, dim });
        // M1 = (P11 + P22)(Q11 + Q22)
        mat M1 = mul(add(P11, P22), add(Q11, Q22));
        // M2 = (P21 + P22)Q11
        mat M2 = mul(add(P21, P22), Q11);
        // M3 = P11(Q12 - Q22)
        mat M3 = mul(P11, subtract(Q12, Q22));
        // M4 = P22(Q21 - Q11)
        mat M4 = mul(P22, subtract(Q21, Q11));
        // M5 = (P11 + P12)Q22
        mat M5 = mul(add(P11, P12), Q22);
        // M6 = (P21 - P11)(Q11 + Q12)
        mat M6 = mul(subtract(P21, P11), add(Q11, Q12));
        // M7 = (P12 - P22)(Q21 + Q22)
        mat M7 = mul(subtract(P12, P22), add(Q21, Q22));
        // C11 = M1 + M4 + M5 - M7
        mat C11 = subtract(add(add(M1, M4), M5), M7);
        // C12 = M3 + M5
        mat C12 = add(M3, M5);
        // C21 = M2 + M4
        mat C21 = add(M2, M4);
        // C22 = M1 + M3 + M6 - M2
        mat C22 = subtract(add(add(M1, M3), M6), M2);
        // Combine and make R using Cs
        eq(R, C11, { 0, 0 });
        eq(R, C12, { 0, dim / 2 });
        eq(R, C21, { dim / 2, 0 });
        eq(R, C22, { dim / 2, dim / 2 });
        // get the non-zero part
        R = block(R, { 0, A.size() }, { 0, B.front().size() });
        // return product
        return R;
    }
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
