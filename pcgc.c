#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

//Função Gradientes Conjugados Pré-condicionados (Jacobi)

int gradienteConjugado(real_t *A, real_t *b, real_t *x, int n, int maxit, double eps, real_t *M, real_t *normaFinal)
{
    real_t *r = malloc(n * sizeof(real_t));  // resíduo
    real_t *z = malloc(n * sizeof(real_t));  // resíduo pré-condicionado
    real_t *p = malloc(n * sizeof(real_t));  // direção de busca
    real_t *Ap = malloc(n * sizeof(real_t)); // A * p
    if (!r || !z || !p || !Ap) {
        fprintf(stderr, "Erro: falha de alocação de memória.\n");
        return -1;
    }

    // --- Passo 1: calcular r = b - A*x ---
    for (int i = 0; i < n; i++) {
        real_t soma = 0.0;
        for (int j = 0; j < n; j++)
            soma += A[i*n + j] * x[j];
        r[i] = b[i] - soma;
    }

    // --- Passo 2: aplicar pré-condicionador M (Jacobi) ---
    for (int i = 0; i < n; i++) {
        if (M != NULL && ABS(M[i]) > 1e-12)
            z[i] = r[i] / M[i];
        else
            z[i] = r[i]; // sem pré-condicionador
    }

    // --- Passo 3: inicializar direção p = z ---
    for (int i = 0; i < n; i++)
        p[i] = z[i];

    // Produto escalar inicial (rᵗz)
    real_t rz_old = 0.0;
    for (int i = 0; i < n; i++)
        rz_old += r[i] * z[i];

    // --- Loop principal ---
    for (int iter = 1; iter <= maxit; iter++) {
        // --- Ap = A * p ---
        for (int i = 0; i < n; i++) {
            real_t soma = 0.0;
            for (int j = 0; j < n; j++)
                soma += A[i*n + j] * p[j];
            Ap[i] = soma;
        }

        // --- calcular alpha ---
        real_t pAp = 0.0;
        for (int i = 0; i < n; i++)
            pAp += p[i] * Ap[i];

        if (ABS(pAp) < 1e-12) {
            fprintf(stderr, "Erro numérico: p^T A p = 0.\n");
            break;
        }

        real_t alpha = rz_old / pAp;

        // --- atualizar x = x + alpha*p ---
        for (int i = 0; i < n; i++)
            x[i] += alpha * p[i];

        // --- atualizar r = r - alpha*Ap ---
        for (int i = 0; i < n; i++)
            r[i] -= alpha * Ap[i];

        // --- verificar convergência: ||r|| < eps ---
        real_t norma_r = 0.0;
        for (int i = 0; i < n; i++)
            norma_r += r[i] * r[i];
        norma_r = sqrt(norma_r);

        if (norma_r < eps) {
            printf("Convergiu em %d iterações. ||r|| = %.6e\n", iter, norma_r);
            free(r); free(z); free(p); free(Ap);
            return iter;
        }

        // --- aplicar pré-condicionador: z = M^{-1} * r ---
        for (int i = 0; i < n; i++) {
            if (M != NULL && abs(M[i]) > 1e-12)
                z[i] = r[i] / M[i];
            else
                z[i] = r[i];
        }

        // --- beta = (rᵗz novo) / (rᵗz antigo) ---
        real_t rz_new = 0.0;
        for (int i = 0; i < n; i++)
            rz_new += r[i] * z[i];

        real_t beta = rz_new / rz_old;
        rz_old = rz_new;

        // --- atualizar direção p = z + beta*p ---
        for (int i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i];
    }

    printf("Aviso: não convergiu após %d iterações.\n", maxit);
    free(r); free(z); free(p); free(Ap);
    return maxit;
}
