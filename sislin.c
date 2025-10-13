#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
 //prosso alterar isso para gerar uma matriz com diagonal principal forte?
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

//função de alocação da matriz principal A
int alocaKDiagonal(int n, real_t ***matriz) {
    //aloca o ponteiro para as linhas
    *matriz = (real_t **)malloc(n * sizeof(real_t *));
    if (*matriz == NULL) {
        fprintf(stderr, "ERRO: Falha ao alocar memória para as linhas da matriz.\n");
        return 0; //retorna falha
    }

    //aloca as colunas
    for (int i = 0; i < n; ++i) {
        //calloc inicializado elementos com 0.0
        (*matriz)[i] = (real_t *)calloc(n, sizeof(real_t));
        if ((*matriz)[i] == NULL) {
            fprintf(stderr, "ERRO: Falha ao alocar memória para as colunas na linha %d.\n", i);
            //se falhar libera tudo que já foi alocado
            for (int j = 0; j < i; ++j) {
                free((*matriz)[j]);
            }
            free(*matriz);
            *matriz = NULL;
            return 0; //retorna falha
        }
    }
    return 1; //
}

int alocaVetorB(int n, real_t **vetor) {
    *vetor = (real_t *)malloc(n * sizeof(real_t));
    if (*vetor == NULL) {
        fprintf(stderr, "ERRO: Falha ao alocar memória para o vetor.\n");
        return 0; //falha
    }
    return 1; //sucesso    
}

//libera memória da matriz kdiagonal
void liberaKDiagonal(int n, real_t **matriz) {
    if (matriz == NULL) return;
    for (int i = 0; i < n; ++i) {
        free(matriz[i]);
    }
    free(matriz);
}

//libera memoria do vetor B 
void liberaVetorB(real_t *vetor) {
    if (vetor == NULL) return;
    free(vetor);
}

//imprimei o destema linear
void imprimeSistema(int n, real_t **A, real_t *B) {
    printf("--- Matriz A ---\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%5.1f ", A[i][j]);
        }
        printf("\n");
    }

    printf("\n--- Vetor B ---\n");
    for (int i = 0; i < n; ++i) {
        printf("%5.1f\n", B[i]);
    }
}



void criaKDiagonal(int n, int k, real_t *A, real_t *B) {
    //preenche a matriz
    int d = (k - 1) / 2;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (abs(i - j) <= d) {
                // Correct way to access element (i, j) in a 1D array
                A[i * n + j] = generateRandomA(i, j, k);
            }
            // Optional: Initialize non-diagonal elements to 0
            else {
                A[i * n + j] = 0.0;
            }
        }
    }
    
    //preenche o vetor
    for (int i = 0; i < n; ++i) {
        // Correct way to access element i in a 1D array
        B[i] = generateRandomB(k);
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t *ASP, real_t *bsp, real_t *tempo)
{
 *tempo = timestamp();

    // Calcula A' = Aᵗ * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            real_t soma = 0.0;
            // Produto escalar da coluna 'i' de A pela coluna 'j' de A
            for (int k2 = 0; k2 < n; ++k2) {
                soma += A[k2 * n + i] * A[k2 * n + j];
            }
            // Acesso direto, sem o ponteiro extra
            ASP[i * n + j] = soma;
        }
    }

    // Calcula b' = Aᵗ * b
    for (int i = 0; i < n; ++i) {
        real_t soma = 0.0;
        // Produto escalar da coluna 'i' de A pelo vetor b
        for (int k2 = 0; k2 < n; ++k2) {
            soma += A[k2 * n + i] * b[k2];
        }
        // Acesso direto
        bsp[i] = soma;
    }

    *tempo = timestamp() - *tempo; 
}


void geraDLU (real_t *A, int n, int k,
              real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
    *tempo = timestamp();

    /* aloca */
    *D = (real_t *) malloc(n * sizeof(real_t));
    *L = (real_t *) malloc((n > 1 ? (n - 1) : 0) * sizeof(real_t));
    *U = (real_t *) malloc((n > 1 ? (n - 1) : 0) * sizeof(real_t));

    if (!(*D) || (n>1 && (!(*L) || !(*U)))) {
        fprintf(stderr, "Erro: falha ao alocar D/L/U em geraDLU().\n");
        if (*D) free(*D);
        if (*L) free(*L);
        if (*U) free(*U);
        *D = *L = *U = NULL;
        *tempo = timestamp() - *tempo;
        return;
    }

    /* preenche */
    for (int i = 0; i < n; ++i) {
        (*D)[i] = A[i * n + i];
    }
    for (int i = 0; i < n - 1; ++i) {
        /* subdiagonal L: A[(i+1), i] */
        (*L)[i] = A[(i + 1) * n + i];
        /* superdiagonal U: A[i, i+1] */
        (*U)[i] = A[i * n + (i + 1)];
    }

    /* proteção: se alguma diagonal for zero ou muito pequena, aumente para epsilon */
    const real_t eps_diag = 1e-12;
    for (int i = 0; i < n; ++i) {
        if (fabs((*D)[i]) < eps_diag) {
            /* aviso e correção */
            // fprintf(stderr, "Warning: D[%d] muito pequeno (%.3e). Corrigindo para %.3e\n", i, (*D)[i], eps_diag);
            (*D)[i] = eps_diag;
        }
    }

    *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k,
                 real_t **M, rtime_t *tempo)
{
    *tempo = timestamp();

    if (w == -1.0) {
        /* sem pré-condicionador */
        *M = NULL;
        *tempo = timestamp() - *tempo;
        return;
    }

    if (w == 0.0) {
        /* Jacobi: devolve vetor diagonal D (com proteção numérica) */
        *M = (real_t *) malloc(n * sizeof(real_t));
        if (!(*M)) {
            fprintf(stderr, "ERRO: falha ao alocar M em geraPreCond()\n");
            *M = NULL;
            *tempo = timestamp() - *tempo;
            return;
        }
        const real_t eps_diag = 1e-12;
        for (int i = 0; i < n; ++i) {
            real_t val = D[i];
            if (fabs(val) < eps_diag) {
                // evitar divisão por zero no precond
                // fprintf(stderr, "Warning: D[%d] muito pequeno (%.3e). Corrigindo para %.3e\n", i, val, eps_diag);
                val = eps_diag;
            }
            (*M)[i] = val;
        }
        *tempo = timestamp() - *tempo;
        return;
    }

    /* Se chegou aqui: Gauss-Seidel / SSOR não implementados (opção BÔNUS) */
    fprintf(stderr, "Aviso: pré-condicionador com w=%.6g não implementado (use -1 ou 0.0).\n", (double)w);
    *M = NULL;
    *tempo = timestamp() - *tempo;
    return;
}


real_t calcResiduoSL (real_t *A, real_t *b, real_t *X,
		      int n, int k, rtime_t *tempo)
{
    rtime_t t0 = timestamp();

    real_t *r = malloc(n * sizeof(real_t));
    if (!r) {
        fprintf(stderr, "Erro de alocação no cálculo do resíduo.\n");
        return -1;
    }

    // Calcula r = b - A*x
    for (int i = 0; i < n; ++i) {
        real_t Ax_i = 0.0;
        for (int j = 0; j < n; ++j) {
            Ax_i += A[i * n + j] * X[j];
        }
        r[i] = b[i] - Ax_i;
    }

    // Calcula norma euclidiana ||r||
    real_t soma = 0.0;
    for (int i = 0; i < n; ++i) {
        soma += r[i] * r[i];
    }
    real_t norma = sqrt(soma);

    free(r);
    *tempo = timestamp() - t0;
    return norma;
}



