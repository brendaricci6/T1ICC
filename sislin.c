#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

//Funções Auxiliares de geração de coeficiente aleatórios

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

//funções de alocação e liberação de memória

//função de alocação da matriz principal A
int alocaKDiagonal(int n, real_t ***matriz) {
    //aloca o ponteiro para as linhas
    *matriz = (real_t **)malloc(n * sizeof(real_t *));
    if (*matriz == NULL) {
        printf("ERRO: Falha ao alocar memória para as linhas da matriz.\n");
        return 0; //retorna falha
    }

    //aloca as colunas
    for (int i = 0; i < n; ++i) {
        //calloc inicializado elementos com 0.0
        (*matriz)[i] = (real_t *)calloc(n, sizeof(real_t));
        if ((*matriz)[i] == NULL) {
            printf("ERRO: Falha ao alocar memória para as colunas na linha %d.\n", i);
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

//aloca vetor B
int alocaVetorB(int n, real_t **vetor) {
    *vetor = (real_t *)malloc(n * sizeof(real_t));
    if (*vetor == NULL) {
        printf("ERRO: Falha ao alocar memória para o vetor.\n");
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
void imprimeSistema(int n, real_t *A, real_t *B) {
    printf("--- Matriz A ---\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%5.1f ", A[i*n + j]);
        }
        printf("\n");
    }

    // printf("\n--- Vetor B ---\n");
    // for (int i = 0; i < n; ++i) {
    //     printf("%5.1f\n", B[i]);
    // }
}


//função principal para gerar a Matriz A KDiagonal e o Vetor B.
void criaKDiagonal(int n, int k, real_t *A, real_t *B) {
    // d é o raio de diagonais
    // (k é o número total de diagonais, ímpar: k = 2d + 1)
    int d = (k - 1) / 2;

    //preenche A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            //verifica se o elemento está dentro da banda diagonal: |i - j| <= d
            if (abs(i - j) <= d) {
                // Preenche com valor aleatório
                A[i * n + j] = generateRandomA(i, j, k);
            }
            // Se estiver fora da banda, o elemento é zero
            else {
                A[i * n + j] = 0.0;
            }
        }
    }
    
    //preenche o vetor
    for (int i = 0; i < n; ++i) {
        B[i] = generateRandomB(k);
    }
}

//função que transforma um sistema Ax=b em um sistema equivalente simétrico e positivo-definido (SPD)
//método padrão para garantir que a matriz A' seja SPD, permitindo o uso do PCG.
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t *ASP, real_t *bsp, real_t *tempo)
{
 *tempo = timestamp(); //inicia medição de tempo 

    // Calcula A' = A~T * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            real_t soma = 0.0;
            //cálculo do elemento (i, j) de AT A: (AT)_i * A_j
            // A coluna 'i' de AT é a linha 'i' de A. A coluna 'j' de A é A_j.
            // O código implementa o produto escalar da coluna 'i' de A pela coluna 'j' de A:
            // (AT A)i,j = Σ_k A[k, i] * A[k, j]
            for (int k2 = 0; k2 < n; ++k2) {
                soma += A[k2 * n + i] * A[k2 * n + j];
            }
            // Acesso direto, sem o ponteiro extra
            ASP[i * n + j] = soma;
        }
    }

    // Calcula b' = AT * b
    for (int i = 0; i < n; ++i) {
        real_t soma = 0.0;
        // Cálculo do elemento 'i' de AT b: (AT)_i * b
        // Produto escalar da coluna 'i' de A pelo vetor b: Σ_k A[k, i] * b[k]
        for (int k2 = 0; k2 < n; ++k2) {
            soma += A[k2 * n + i] * b[k2];
        }
        // Acesso direto
        bsp[i] = soma;
    }

    *tempo = timestamp() - *tempo; //finaliza medição de tempo 
}

//funç~es de pré-condicionamento

//Gera os vetores de Decomposição DLU para uma matriz A.
//D: Diagonal principal
//L: Sub-diagonal
//U: Super-diagonal
void geraDLU (real_t *A, int n, int k,
              real_t *D, real_t *L, real_t *U, rtime_t *tempo, double eps)
{
    *tempo = timestamp();

    /* preenche o vetor D com a diagonal principal*/
    for (int i = 0; i < n; ++i) {
        D[i] = A[i * n + i];
    }
    //preenche os vetores L e U com as diagonais vizinhas (sub e super)
    for (int i = 0; i < n - 1; ++i) {
        /* subdiagonal L: A[(i+1), i] */
        L[i] = A[(i + 1) * n + i];
        /* superdiagonal U: A[i, i+1] */
        U[i] = A[i * n + (i + 1)];
    }

    /* proteção: se alguma diagonal for zero ou muito pequena, aumente para epsilon */
    for (int i = 0; i < n; ++i) {
        if (ABS(D[i]) < eps) {
            /* aviso e correção */
            // printf("Warning: D[%d] muito pequeno (%.3e). Corrigindo para %.3e\n", i, (*D)[i], eps);
            D[i] = eps; //define um valor mínimo (eps) para evitar instabilidade
        }
    }

    *tempo = timestamp() - *tempo;
}

/**
 * Devolve o pré-condicionador M.
 * w (omega) é o parâmetro de relaxamento.
 * w = 0.0 -> Pré-condicionador Jacobi (M = D)
 * w = -1.0 -> Sem pré-condicionador (M = I)
 * Outros valores de w podem indicar SOR ou SSOR (não implementados aqui)
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t *M, rtime_t *tempo, double eps)
{
    *tempo = timestamp();

    //caso 1: w = -1.0 -> Sem pré-condicionador (M = I)
    // função retorna sem preencher M.
    if (w == -1.0) {
        *tempo = timestamp() - *tempo;
        return;
    }

    //caso 2: w = 0.0 -> Pré-condicionador Jacobi (M = D)
    if (w == 0.0) {
        /* Jacobi: devolve vetor diagonal D (com proteção numérica) */
        if (!M) {
            *tempo = timestamp() - *tempo;
            printf("ERRO: falha ao alocar M em geraPreCond()\n");
            return;
        }
        
        // Define M como a diagonal D
        for (int i = 0; i < n; ++i) {
            real_t val = D[i];
            if (ABS(val) < eps) {
                // evitar divisão por zero no precond
                // printf("Warning: D[%d] muito pequeno (%.3e). Corrigindo para %.3e\n", i, val, eps);
                val = eps;
            }
            M[i] = val; // M armazena a diagonal D
        }
        *tempo = timestamp() - *tempo;
        return;
    }

    //
    printf("pré-condicionador com w=%lf não implementado (usar -1 ou 0.0).\n", w);
    *tempo = timestamp() - *tempo;
    return;
}

//Funções de pós-Processamento

//Calcula a Norma do Resíduo Euclidiano (||r||2 = ||b - A*X||2) para avaliar a precisão
real_t calcResiduoSL (real_t *A, real_t *b, real_t *X,
		      int n, int k, rtime_t *tempo)
{
    rtime_t t0 = timestamp(); // Inicia a medição de tempo

    real_t *r = malloc(n * sizeof(real_t));
    if (!r) {
        printf("Erro de alocação no cálculo do resíduo.\n");
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



