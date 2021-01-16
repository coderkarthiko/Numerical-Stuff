#pragma once

#include "types.h"

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
        if (M[i][i] <= threshold) { // if pivot element is 0, det(M) = 0
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
