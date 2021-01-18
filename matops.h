#include "types.h"

bool strassen = true;

// print vector
void vdprint(vd v) {
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
void matprint(mat M) {
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
        for (int j = 0; j < AT.front().size(); j++) {
            AT[i][j] = A[j][i];
        }
    }
    return AT;
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

// concatenate A and B row-wise
mat rconcat(mat A, mat B) {
    assert(A.front().size() == B.front().size());
    for (auto r : B) {
        A.push_back(r);
    }
    return A;
}

// concatenate A and B column-wise
mat cconcat(mat A, mat B) {
    assert(A.size() == B.size());
    for (int i = 0; i < A.size(); i++) {
        for (auto c : B[i]) {
            A[i].push_back(c);
        }
    }
    return A;
}

// vector dot product
ld dot(vd A, vd B) {
    ld prod = 0;
    for (int i = 0; i < A.size(); i++) {
        prod += A[i] * B[i];
    }
    return prod;
}

vd matrow(mat A, int r) {
    return A[r];
}

vd matcol(mat A, int c) {
    vd col = {};
    for (auto r : A) {
        col.push_back(r[c]);
    }
    return col;
}

// for matrices of dims AxB and CxD, if max(A, B, C, D) < 512 - O(N^3) mat mul algorithm
// if max(A, B, C, D) > 512 - Strassen's ~O(2^2.8074) mat mul algorithm
// for disproportionate mat sizes - split, mul then concatenate
mat mul(mat A, mat B) {
    int mdim = max(A.size(), B.front().size());
    if (!strassen || mdim < 512) {
        mat C(A.size(), vd(B.front().size()));
        for (int i = 0; i < A.size(); i++) {
            for (int j = 0; j < B.front().size(); j++) {
                C[i][j] = dot(matrow(A, i), matcol(B, j));
            }
        }
        return C;
    } 
    else {
        mdim += mdim % 2;
        int dim = mdim / 2;
        mat P(mdim, vd(mdim, 0)), Q(mdim, vd(mdim, 0)), R(mdim, vd(mdim, 0));
        eq(P, A, { 0, 0 });
        eq(Q, B, { 0, 0 });
        mat P11 = block(P, { 0, dim }, { 0, dim });
        mat P12 = block(P, { 0, dim }, { dim, mdim });
        mat P21 = block(P, { dim, mdim }, { 0, dim });
        mat P22 = block(P, { dim, mdim }, { dim, mdim });
        mat Q11 = block(Q, { 0, dim }, { 0, dim });
        mat Q12 = block(Q, { 0, dim }, { dim, mdim });
        mat Q21 = block(Q, { dim, mdim }, { 0, dim });
        mat Q22 = block(Q, { dim, mdim }, { dim, mdim });
        mat M1 = mul(add(P11, P22), add(Q11, Q22));
        mat M2 = mul(add(P21, P22), Q11);
        mat M3 = mul(P11, subtract(Q12, Q22));
        mat M4 = mul(P22, subtract(Q21, Q11));
        mat M5 = mul(add(P11, P12), Q22);
        mat M6 = mul(subtract(P21, P11), add(Q11, Q12));
        mat M7 = mul(subtract(P12, P22), add(Q21, Q22));
        eq(R, subtract(add(add(M1, M4), M7), M5), { 0, 0 });
        eq(R, add(M3, M5), { 0, dim });
        eq(R, add(M2, M4), { dim, 0 });
        eq(R, subtract(add(add(M1, M3), M6), M2), { dim, dim });
        return block(R, { 0, A.size() }, { 0, B[0].size() });
    }
}

// input a vector
vd getvd() {
    int n;
    cin >> n;
    vd v(n);
    for (int i = 0; i < n; i++) {
        std::cin >> v[i];
    }
    return v;
}

// input a matrix
mat getmat() {
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
