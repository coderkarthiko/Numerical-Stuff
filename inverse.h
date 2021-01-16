#pragma once

#include "types.h"

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
