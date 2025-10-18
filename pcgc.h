#ifndef PCGC_H
#define PCGC_H

#include "utils.h"

/**
 * Conjugate Gradient pré-condicionado (suporta M=NULL ou M=diagonal)
 *
 * A: matriz ASC (n*n) armazenada linearmente (A[i*n + j])
 * b: vetor RHS (n)
 * x: vetor solução (entrada: inicial guess; saída: solução). Aqui assumimos x inicial = 0.
 * n: dimensão
 * k: número de diagonais (ímpar)
 * eps: tolerância (critério ||x - x_prev||_inf < eps)
 * maxit: número máximo de iterações
 * residuo_out: saída (norma L2 do resíduo final)
 * M: ponteiro para pré-condicionador. Se NULL -> M = I. Caso Jacobi, M aponta para vetor de dimensão n contendo os elementos da diagonal (D_i).
 * norma_inf_out: saída (norma infinita entre últimas iterações)
 *
 * Retorna: número de iterações realizadas (>0) em caso de sucesso,
 *          0 se solução inicial já é correta (r==0),
 *         -1 em caso de quebra numérica (p^T A p == 0 ou divisão por zero no pré-condicionador).
 */
int gradienteConjugado(real_t *A, real_t *b, real_t *x, int n, int maxit, double eps, real_t *M, real_t *normaFinal, rtime_t *tempoIter);

#endif
