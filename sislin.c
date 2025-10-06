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



void criaKDiagonal(int n, int k, real_t ***A, real_t **B) {
    //preenche a matriz
    int d = (k - 1) / 2;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (abs(i - j) <= d) {
                (*A)[i][j] = generateRandomA(i, j, k);
            }
        }
    }
    
    //preenche o vetor
    for (int i = 0; i < n; ++i) {
        (*B)[i] = generateRandomB(k);
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t **ASP, real_t **bsp, real_t *tempo)
{
  *tempo = timestamp();

  *tempo = timestamp() - *tempo;
 
}


void geraDLU (real_t *A, int n, int k,
	      real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();


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


  *tempo = timestamp() - *tempo;
}


real_t calcResiduoSL (real_t *A, real_t *b, real_t *X,
		      int n, int k, rtime_t *tempo)
{
  *tempo = timestamp();

  real_t *r = calloc(n, sizeof(real_t));

  

  *tempo = timestamp() - *tempo;
}



