#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"

//gradiente Conjugado Pré condicionado
int gradienteConjugado(real_t *A, real_t *b, real_t *x, int n, int maxit, double eps, real_t *M, real_t *normaFinal, rtime_t *tempoIter)
{
    //constantes do formato dia
    const int k_ASP = N_DIAG;       // 7
    const int offset_center = OFFSET_CENTER; // 3

    //alocação dos vetores auxiliares
    real_t *r = malloc(n * sizeof(real_t));
    real_t *z = malloc(n * sizeof(real_t));
    real_t *p = malloc(n * sizeof(real_t));
    real_t *Ap = malloc(n * sizeof(real_t));

    if (!r || !z || !p || !Ap) return -1;

    //calcular resíduo inicial,
    
    // Inicializa r com b
    for(int i=0; i<n; ++i) r[i] = b[i];

    // Subtrai A*x (Usando a lógica otimizada de diagonais)
    for (int diag_idx = 0; diag_idx < k_ASP; diag_idx++) {
        int offset = diag_idx - offset_center;
        
        // Calcula limites fora do loop para evitar IFs internos
        int inicio = (offset < 0) ? -offset : 0;
        int fim    = (offset > 0) ? n - offset : n;
        
        real_t *diagonal = &A[diag_idx * n]; // Ponteiro para o início da diagonal

        // Loop "limpo" para o compilador vetorizar (Auto-vectorization)
        for (int i = inicio; i < fim; i++) {
            r[i] -= diagonal[i] * x[i + offset];
        }
    }

    //Pré-condicionador -
    for (int i = 0; i < n; i++) {
        z[i] = (M != NULL) ? r[i] / M[i] : r[i];
        p[i] = z[i];
    }

    //calculo do produto escalar inicial
    real_t rz_old = 0.0;
    for (int i = 0; i < n; i++) rz_old += r[i] * z[i];

    int iter;
    *tempoIter = timestamp();

    // Marcador LIKWID)
    LIKWID_MARKER_START("op1");

    for (iter = 1; iter <= maxit; iter++) {
        
        //Loop externo nas diagonais
        
        //zera o vetor Ap
        for(int i=0; i<n; ++i) Ap[i] = 0.0;

        // acumula contribuição de cada diagonal
        for (int diag_idx = 0; diag_idx < k_ASP; diag_idx++) {
            int offset = diag_idx - offset_center;
            
            //Definição dos limites
            int inicio = (offset < 0) ? -offset : 0;
            int fim    = (offset > 0) ? n - offset : n;
            
            real_t *diagonal = &A[diag_idx * n];

            //
            for (int i = inicio; i < fim; i++) {
                Ap[i] += diagonal[i] * p[i + offset];
            }
        }

        // produto escalar
        real_t pAp = 0.0;
        for (int i = 0; i < n; i++) pAp += p[i] * Ap[i];

        // Verificação de divide por zero
        if (fabs(pAp) < 1e-15) break; 

        real_t alpha = rz_old / pAp;

        // Atualizade X e R
        real_t norma_r_sq = 0.0;
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            norma_r_sq += r[i] * r[i];
        }

        real_t norma_r = sqrt(norma_r_sq);
        *normaFinal = norma_r;

        // critério de parada
        if (norma_r < eps) break;

        // Aplica Precondicionador 
        if (M) {
            for (int i = 0; i < n; i++) z[i] = r[i] / M[i];
        } else {
            for (int i = 0; i < n; i++) z[i] = r[i];
        }

        //calculo de Beta
        real_t rz_new = 0.0;
        for (int i = 0; i < n; i++) rz_new += r[i] * z[i];

        real_t beta = rz_new / rz_old;
        rz_old = rz_new;

        //atualiza direção p
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }
    }

    LIKWID_MARKER_STOP("op1");

    *tempoIter = timestamp() - *tempoIter;
    if (iter > 0) *tempoIter = *tempoIter / iter;

    free(r); free(z); free(p); free(Ap);
    return iter;
}