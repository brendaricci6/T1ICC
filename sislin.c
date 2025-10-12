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

    // --- Aloca vetores ---
    *D = (real_t *) malloc(n * sizeof(real_t));
    *L = (real_t *) malloc((n - 1) * sizeof(real_t));
    *U = (real_t *) malloc((n - 1) * sizeof(real_t));

    if (!(*D) || !(*L) || !(*U)) {
        fprintf(stderr, "Erro: falha ao alocar D, L ou U em geraDLU().\n");
        return;
    }

    // --- Inicializa com zeros ---
    for (int i = 0; i < n; i++) {
        (*D)[i] = 0.0;
    }
    
    for (int i = 0; i < n - 1; i++) {
        (*L)[i] = 0.0;
        (*U)[i] = 0.0;
    }

    // --- Preenche D, L e U ---
    // Considerando o formato compactado: cada linha i tem k elementos
    // com a diagonal principal no meio: offset = k/2
    int offset = k / 2;

    for (int i = 0; i < n; i++) {
        int base = i * k;

        // Diagonal principal
        (*D)[i] = A[base + offset];

        // Subdiagonal (elemento logo abaixo da diagonal principal)
        if (i > 0)
            (*L)[i - 1] = A[base + offset - 1];

        // Superdiagonal (elemento logo acima da diagonal principal)
        if (i < n - 1)
            (*U)[i] = A[base + offset + 1];
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


    // Aloca M (matriz densa n x n)
    *M = (real_t *) calloc(n * n, sizeof(real_t));
    if (*M == NULL) {
        fprintf(stderr, "ERRO: Falha ao alocar memória para M.\n");
        exit(1);
    }

    // Preenche M usando fórmula do pré-condicionador
    for (int i = 0; i < n; i++) {
        // diagonal principal
        (*M)[i*n + i] = (D[i] + w * L[i]) / D[i] * (D[i] + w * U[i]);

        // subdiagonal
        if (i > 0)
            (*M)[i*n + (i-1)] = (D[i] + w * L[i]) / D[i] * (w * U[i-1]);

        // superdiagonal
        if (i < n - 1)
            (*M)[i*n + (i+1)] = (w * L[i+1]) * (D[i] + w * U[i]) / D[i];
    }

  *tempo = timestamp() - *tempo;
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



