#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "sislin.h" // Deve incluir N_DIAG (7) e OFFSET_CENTER (3)

//Função Gradientes Conjugados Pré-condicionados (Jacobi)
//A: Matriz do sistema (ASP, no formato DIA n x k)
//b: Vetor do lado direito (bsp)
// x: Vetor solução (entrada: palpite inicial; saída: solução aproximada)
// n: Dimensão do sistema
// maxit: Número máximo de iterações (retorno da função)
// eps: Tolerância de convergência (erro absoluto máximo ||r||)
// M: é a Diagonal de A
// normaFinal: Ponteiro para armazenar a norma do resíduo final
int gradienteConjugado(real_t *A, real_t *b, real_t *x, int n, int maxit, double eps, real_t *M, real_t *normaFinal, rtime_t *tempoIter)
{
    // Constantes do formato DIA
    const int k_ASP = N_DIAG;       // 7
    const int offset_center = OFFSET_CENTER; // 3

    // Alocação dos vetores auxiliares
    real_t *r = malloc(n * sizeof(real_t));  // resíduo
    real_t *z = malloc(n * sizeof(real_t));  // resíduo pré-condicionado
    real_t *p = malloc(n * sizeof(real_t));  // direção de busca
    real_t *Ap = malloc(n * sizeof(real_t)); // A * p

    // verifica a alocação de memória
    if (!r || !z || !p || !Ap) {
        printf("Erro: falha de alocação de memória.\n");
        return -1;
    }
    
    // ==========================================================
    // Passo 1: calcular o resíduo inicial r = b - A*x (CORRIGIDO PARA DIA)
    // ==========================================================
    for (int i = 0; i < n; i++) {
        real_t soma = 0.0;
        
        // Multiplicação otimizada da linha i de A pelo vetor x (Formato DIA)
        for (int diag_idx = 0; diag_idx < k_ASP; diag_idx++) {
            int offset = diag_idx - offset_center; // -3, -2, ..., 3
            int j = i + offset; 
            
            // Checa limites de j: 0 <= j < n
            if (j >= 0 && j < n) {
                // Acesso corrigido: A[linha * k_ASP + indice_diag]
                soma += A[i * k_ASP + diag_idx] * x[j];
            }
        }
        r[i] = b[i] - soma;
    }

    // Passo 2: aplicar pré-condicionador M (Jacobi)
    // o pré-condicionador é aplicado para obter o resídulo pré-condicionado 'z'
    for (int i = 0; i < n; i++) {
        if (M != NULL && ABS(M[i]) > __DBL_EPSILON__)
            z[i] = r[i] / M[i];
        else
            z[i] = r[i]; // sem pré-condicionador
    }

    // Passo 3: inicializar direção de busca p = z
    for (int i = 0; i < n; i++)
        p[i] = z[i];

    // Produto escalar inicial (rᵗz)
    real_t rz_old = 0.0;
    for (int i = 0; i < n; i++)
        rz_old += r[i] * z[i];

    // ============== Loop principal ==============
    int iter;
    *tempoIter = timestamp();
    for (iter = 1; iter <= maxit; iter++) {
        
        // ==========================================================
        // Ap = A * p (CORRIGIDO PARA DIA)
        // ==========================================================
        for (int i = 0; i < n; i++) {
            real_t soma = 0.0;
            
            // Loop otimizado sobre as diagonais de A (ASP)
            for (int diag_idx = 0; diag_idx < k_ASP; diag_idx++) {
                int offset = diag_idx - offset_center; // -3, -2, ..., 3
                int j = i + offset; // Coluna j correspondente ao elemento na linha i
                
                // Checa limites de j: 0 <= j < n
                if (j >= 0 && j < n) {
                    // Acesso corrigido: A[linha * k_ASP + indice_diag]
                    soma += A[i * k_ASP + diag_idx] * p[j];
                }
            }
            Ap[i] = soma;
        }

        // calcular alpha (tamanho do passo)
        real_t pAp = 0.0;
        for (int i = 0; i < n; i++)
            pAp += p[i] * Ap[i];

        // tratamento se for zero
        if (ABS(pAp) < __DBL_EPSILON__) {
            break;
        }

        real_t alpha = rz_old / pAp;

        // atualizar x = x + alpha*p
        for (int i = 0; i < n; i++)
            x[i] += alpha * p[i];

        // atualizar r = r - alpha*Ap
        for (int i = 0; i < n; i++)
            r[i] -= alpha * Ap[i];

        // verificar convergência: ||r|| < eps
        real_t norma_r = 0.0;
        for (int i = 0; i < n; i++)
            norma_r += r[i] * r[i]; 
        norma_r = sqrt(norma_r); 
        *normaFinal = norma_r;

        if (norma_r < eps) {
            *tempoIter = timestamp() - *tempoIter;
            *tempoIter = *tempoIter / iter;
            free(r); 
            free(z); 
            free(p); 
            free(Ap);
            return iter;
        }

        // aplicar pré-condicionador: z = M^{-1} * r
        for (int i = 0; i < n; i++) {
            if (M != NULL && ABS(M[i]) > __DBL_EPSILON__)
                z[i] = r[i] / M[i];
            else
                z[i] = r[i];
        }

        // Calcular o fator beta (Gram-Schmidt)
        real_t rz_new = 0.0;
        for (int i = 0; i < n; i++)
            rz_new += r[i] * z[i];

        real_t beta = rz_new / rz_old;
        rz_old = rz_new;

        // atualizar direção p = z + beta*p
        for (int i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i];
    }
    
    // Se o loop terminar sem atingir a tolerância
    *tempoIter = timestamp() - *tempoIter;
    *tempoIter = *tempoIter / iter;
    free(r); 
    free(z); 
    free(p); 
    free(Ap);
    return iter;
}