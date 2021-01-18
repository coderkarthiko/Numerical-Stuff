#include "matops.h"

vd linsolve(mat A, vd y) { 
    vd x;
    for (auto r : mul(inv(A), transpose({ y }))) {
        x.push_back(r[0]);
    }
    return x;
}
