#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"

//Funções Auxiliares de geração de coeficiente aleatórios
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

//Funções Principais
//Função que gera os coeficientes de um sistema linear k-diagonal V2)
void criaKDiagonal(int n, int k, real_t *A, real_t *B) {
    //d é o raio de diagonais
    //(k é o número total de diagonais, ímpar: k = 2d + 1)
    int d = (k - 1) / 2;
    int diag_offset_center = d; //indice da diagonal principal na matriz de armazenamento (d)

    //per corre cada uma das k diagonais
    for (int diag_idx = 0; diag_idx < k; ++diag_idx) {
        int offset = diag_idx - diag_offset_center;
        
        //define onde começa e termina o loop para não sair da matriz
        int i_start = (offset < 0) ? -offset : 0;
        int i_end   = (offset > 0) ? n - offset : n;

        for (int i = i_start; i < i_end; ++i) {
            int j = i + offset;
            //layout V]
            A[diag_idx * n + i] = generateRandomA(i, j, k);
        }
    }

    //preenche o vetor B
    for (int i = 0; i < n; ++i) {
        B[i] = generateRandomB(k);
    }
}

//Função Simetrica Positiva 
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t *ASP, real_t *bsp, real_t *tempo)
{
    *tempo = timestamp();
    //printf("genSimetricPositiva");
    
    //Recupera constantes do formato
    int d_A = (k - 1) / 2;
    int k_ASP = N_DIAG;
    int offset_center = OFFSET_CENTER; 
    
    
    //zera ASP
    for(int i=0; i < n * k_ASP; i++) ASP[i] = 0.0;

    int d_ASP = (k_ASP - 1) / 2; //3

    //Preenche mantendo estrutura
    for (int i = 0; i < n; ++i) {
        int j_start = (i - d_ASP < 0) ? 0 : i - d_ASP;
        int j_end   = (i + d_ASP >= n) ? n - 1 : i + d_ASP;

        for (int j = j_start; j <= j_end; ++j) {
            real_t soma = 1.0; //valor ara teste de desempenho
            
            int offset = j - i;
            int diag_idx = offset + offset_center;
            
            //armazenamento Otimizado
            ASP[diag_idx * n + i] = soma; 
        }
    }
    
    //Copia b para bsp
    for(int i=0; i<n; i++) bsp[i] = b[i];

    *tempo = timestamp() - *tempo;
}

//gera vetores de Decomposião DLU
void geraDLU (real_t *A, int n, int k, real_t *D, real_t *L, real_t *U, rtime_t *tempo, double eps)
{
    *tempo = timestamp();
    int d = (k - 1) / 2; //raio

    //diagonal principal (
    //indice d no array de diagonais
    real_t *diagPrincipal = &A[d * n]; 
    for (int i = 0; i < n; ++i) {
        D[i] = diagPrincipal[i];
        if (fabs(D[i]) < eps) D[i] = eps;
    }

    //subdiagonal
    //índice d-1
    real_t *diagL = &A[(d - 1) * n];
    //Superdiagonal
    // índice d+1
    real_t *diagU = &A[(d + 1) * n];

    for (int i = 0; i < n - 1; ++i) {
        L[i] = diagL[i + 1]; //L[i] = A[i+1][i]
        U[i] = diagU[i];     //U[i]= A[i][i+1]
    }
    *tempo = timestamp() - *tempo;
}

//Devolve o pré-condicionador M
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t *M, rtime_t *tempo, double eps)
{
    *tempo = timestamp();
    if (w == 0.0 && M) { 
        for (int i = 0; i < n; ++i) M[i] = (fabs(D[i]) < eps) ? eps : D[i];
    }
    *tempo = timestamp() - *tempo;
}

//Calcula a Norm do Resíduo Euclidiano
real_t calcResiduoSL (real_t *A, real_t *b, real_t *X, int n, int k, rtime_t *tempo)
{
    rtime_t t0 = timestamp(); //inicia a medição de tempo
    
    //Marcador LIKWI
    LIKWID_MARKER_START("op2");

    real_t *r = malloc(n * sizeof(real_t));
    if (!r) return -1.0;
    
    //Inicializar r com o vetor b
    for(int i=0; i<n; ++i) r[i] = b[i];

    int d = (k - 1) / 2; //Raio de banda
    int diag_offset_center = d;

    //Loop das Diagonais 
    // Loop Interchange
    //percorre as diagonais externamente para ter acesso sequencial na memória
    for (int diag_idx = 0; diag_idx < k; ++diag_idx) {
        int offset = diag_idx - diag_offset_center;
        
        //pega o ponteiro para o começo dessa diagonal
        real_t *diagonalAtual = &A[diag_idx * n];

        //Calcula limites 
        int i_start = (offset < 0) ? -offset : 0;
        int i_end   = (offset > 0) ? n - offset : n;

        //Loop )
        for (int i = i_start; i < i_end; ++i) {
            //r[i] = r[i] - A[i][j] * X[j]
            r[i] -= diagonalAtual[i] * X[i + offset];
        }
    }

    //calcula norma euclidiana
    real_t soma = 0.0;
    for (int i = 0; i < n; ++i) {
        soma += r[i] * r[i];
    }
    real_t norma = sqrt(soma);

    LIKWID_MARKER_STOP("op2");
    
    free(r);
    *tempo = timestamp() - t0;
    return norma;
}