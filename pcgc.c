#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
#include "utils.h"

/**
 * Multiplicação de matriz k-diagonal compacta por vetor: y = A*x
 */
static void matVecKDiagonal(real_t *A, real_t *x, real_t *y, int n, int k) {
    int d = k / 2;
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        int start = (i - d) > 0 ? (i - d) : 0;
        int end   = (i + d) < (n - 1) ? (i + d) : (n - 1);
        for (int j = start; j <= end; j++) {
            y[i] += A[i * k + (j - i + d)] * x[j];
        }
    }
}

/**
 * Aplica o pré-condicionador M: z = M^-1 * r
 * Para M denso (n x n), fazemos multiplicação direta
 */
static void applyPreCond(real_t *M, real_t *r, real_t *z, int n) {
    for (int i = 0; i < n; i++) {
        z[i] = 0.0;
        for (int j = 0; j < n; j++) {
            z[i] += M[i * n + j] * r[j];
        }
    }
}

/**
 * Método de Gradientes Conjugados
 *
 * A : matriz k-diagonal compacta
 * b : vetor de termos independentes
 * x : vetor solução inicial (chute)
 * n : dimensão
 * k : número de diagonais
 * tol : tolerância
 * maxit : número máximo de iterações
 * residuo : saída do resíduo final
 * M : pré-condicionador (pode ser NULL para sem pré-condicionador)
 *
 * Retorna número de iterações executadas
 */
int conjugateGradient(real_t *A, real_t *b, real_t *x,
                      int n, int k, double tol, int maxit,
                      real_t *residuo, real_t *M) {
    real_t *r = (real_t *) malloc(n * sizeof(real_t));
    real_t *p = (real_t *) malloc(n * sizeof(real_t));
    real_t *Ap = (real_t *) malloc(n * sizeof(real_t));
    real_t *z = (real_t *) malloc(n * sizeof(real_t));
    if (!r || !p || !Ap || !z) {
        fprintf(stderr, "Erro: falha ao alocar memória em CG.\n");
        exit(1);
    }

    // r0 = b - A*x0
    matVecKDiagonal(A, x, Ap, n, k);
    for (int i = 0; i < n; i++)
        r[i] = b[i] - Ap[i];

    // z0 = M^-1 * r0
    if (M)
        applyPreCond(M, r, z, n);
    else
        for (int i = 0; i < n; i++)
            z[i] = r[i];

    // p0 = z0
    for (int i = 0; i < n; i++)
        p[i] = z[i];

    double rz_old = 0.0;
    for (int i = 0; i < n; i++)
        rz_old += r[i] * z[i];

    int iter;
    for (iter = 0; iter < maxit; iter++) {
        matVecKDiagonal(A, p, Ap, n, k);

        double alpha = rz_old;
        double pAp = 0.0;
        for (int i = 0; i < n; i++)
            pAp += p[i] * Ap[i];
        alpha /= pAp;

        // x = x + alpha * p
        for (int i = 0; i < n; i++)
            x[i] += alpha * p[i];

        // r = r - alpha * A*p
        for (int i = 0; i < n; i++)
            r[i] -= alpha * Ap[i];

        // verifica tolerância
        double max_err = 0.0;
        for (int i = 0; i < n; i++)
            if (fabs(alpha * p[i]) > max_err)
                max_err = fabs(alpha * p[i]);
        if (max_err < tol)
            break;

        // z = M^-1 * r
        if (M)
            applyPreCond(M, r, z, n);
        else
            for (int i = 0; i < n; i++)
                z[i] = r[i];

        double rz_new = 0.0;
        for (int i = 0; i < n; i++)
            rz_new += r[i] * z[i];

        double beta = rz_new / rz_old;

        // p = z + beta * p
        for (int i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i];

        rz_old = rz_new;
    }

    // calcula resíduo final
    *residuo = 0.0;
    for (int i = 0; i < n; i++)
        *residuo += r[i] * r[i];
    *residuo = sqrt(*residuo);

    free(r);
    free(p);
    free(Ap);
    free(z);

    return iter + 1;
}
