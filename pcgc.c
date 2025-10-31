#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

//Função Gradientes Conjugados Pré-condicionados (Jacobi)
//A: Matriz do sistema
//b: Vetor do lado direito
// x: Vetor solução (entrada: palpite inicial; saída: solução aproximada)
// n: Dimensão do sistema
// maxit: Número máximo de iterações (retorno da função)
// eps: Tolerância de convergência (erro absoluto máximo ||r||)
//M: é a Diagonal de A
// normaFinal: Ponteiro para armazenar a norma do resíduo final
int gradienteConjugado(real_t *A, real_t *b, real_t *x, int n, int maxit, double eps, real_t *M, real_t *normaFinal, rtime_t *tempoIter)
{
    // Alocação dos vetores auxiliares
    real_t *r = malloc(n * sizeof(real_t));  // resíduo
    real_t *z = malloc(n * sizeof(real_t));  // resíduo pré-condicionado
    real_t *p = malloc(n * sizeof(real_t));  // direção de busca
    real_t *Ap = malloc(n * sizeof(real_t)); // A * p

    //verifica a alocação de memória
    if (!r || !z || !p || !Ap) {
        printf("Erro: falha de alocação de memória.\n");
        return -1; //retorna o erro
    }
    
    // Passo 1: calcular o resíduo inicial r = b - A*x
    for (int i = 0; i < n; i++) {
        real_t soma = 0.0;
        //multiplicação da linha i de A pelo vetor x
        for (int j = 0; j < n; j++)
            soma += A[i*n + j] * x[j]; //matriz A
        r[i] = b[i] - soma;
    }

    // Passo 2: aplicar pré-condicionador M (Jacobi)
    //o pré-condicionador é aplicado para obter o resídulo pré-condicionado 'z'
    //M é a diagonal de A
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
        // Ap = A * p
        //Multiplicação da matriz A pelo vetor de direção de busca p
        for (int i = 0; i < n; i++) {
            real_t soma = 0.0;
            for (int j = 0; j < n; j++)
                soma += A[i*n + j] * p[j];
            Ap[i] = soma;
        }

        // calcular alpha (tamanho do passo)
        real_t pAp = 0.0;
        for (int i = 0; i < n; i++)
            pAp += p[i] * Ap[i];

        //tratamento se for zero
        if (ABS(pAp) < __DBL_EPSILON__) {
            // printf("Erro numérico: p^T A p = 0.\n");
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
            norma_r += r[i] * r[i]; //soma dos quadrados
        norma_r = sqrt(norma_r); //raiz quadrada
        // atribui valor do erro 
        *normaFinal = norma_r;

        if (norma_r < eps) {
            // atribui valor do erro 
            // calcula tempo medio das iteracoes ate agora
            *tempoIter = timestamp() - *tempoIter;
            *tempoIter = *tempoIter / iter;
            free(r); 
            free(z); 
            free(p); 
            free(Ap);
            return iter; //retorna número de iterações
        }

        // aplicar pré-condicionador: z = M^{-1} * r
        //calcula o novo resíduo pré-condicionado
        for (int i = 0; i < n; i++) {
            if (M != NULL && ABS(M[i]) > __DBL_EPSILON__)
                z[i] = r[i] / M[i];
            else
                z[i] = r[i];
        }

        // Calcular o fator beta (Gram-Schmidt)
        // beta = ((r^T)z novo) / ((r~T)z antigo)
        real_t rz_new = 0.0;
        for (int i = 0; i < n; i++)
            rz_new += r[i] * z[i];

        real_t beta = rz_new / rz_old;
        rz_old = rz_new; //atualiza

        // atualizar direção p = z + beta*p
        for (int i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i];
    }
    // calcula tempo medio das iteracoes ate agora
    *tempoIter = timestamp() - *tempoIter;
    *tempoIter = *tempoIter / iter;
    // Se o loop terminar sem atingir a tolerância
    
    //printf("Aviso: não convergiu após %d iterações.\n", maxit);
    free(r); 
    free(z); 
    free(p); 
    free(Ap);
    return iter;
}
