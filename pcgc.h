
#include "utils.h"

int conjugateGradient(real_t *A, real_t *b, real_t *x,
                      int n, int k, double tol, int maxit,
                      real_t *residuo, real_t *M);