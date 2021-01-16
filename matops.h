#include "types.h"

// mat multiplication
// for N*M >= 2^12, uses ~O(N^2.8074) Strassen's algorithm 
// otherwise, use O(N^3) version
mat mul(mat A, mat B) { 
    mat AB(A.size(), vd(B[0].size()));
    for (int i = 0; i < AB.size(); i++) {
        for (int j = 0; j < AB[0].size(); j++) {
            for (int k = 0; k < AB.size(); k++) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return AB;
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

// determinant of NxN matrix
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
                std::swap(M[i][j], M[r][j]);
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

// inverse of NxN matrix
mat inv(mat M) {
    // make identity matrix of same dim as M
    mat I(M.size(), vd(M.size()));
    for (int i = 0; i < M.size(); i++) {
        I[i][i] = 1.;
    }
    for (int i = 0; i < M.size(); i++) {
        int r = i; // partial pivoting 
        for (int j = i + 1; j < M.size(); j++) {
            if (std::abs(M[j][i]) > std::abs(M[r][i])) {
                r = j;
            }
        }
        swap(M[i], M[r]);
        swap(I[i], I[r]);
        double f = M[i][i];
        if (f <= threshold) {
            std::cout << "NO INVERSE";
            exit(EXIT_FAILURE);
        }
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
    std::cin >> n;
    vd v(n);
    for (int i = 0; i < n; i++) {
        std::cin >> v[i];
    }
    return v;
}

// input a matrix
mat get_mat() {
    int n, m;
    std::cin >> n >> m;
    mat M(n, vd(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cin >> M[i][j];
        }
    }
    return M;
}

// print vector
void print_vd(vd v) {
    for (auto vx : v) {
        std::cout << vx << " ";
    }
}

// print matrix
void print_mat(mat M) {
    for (auto r : M) {
        for (auto c : r) {
            std::cout << c << " ";
        }
        std::cout << "\n";
    }
}
