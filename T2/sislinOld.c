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

// Função auxiliar para acessar o elemento (i, j) da matriz DIA
// A: ponteiro para o início da matriz n x 7 (armazenada em formato de matriz)
// i, j: coordenadas
// n: dimensão
// Retorna: o valor A[i][j] (0.0 se fora da banda)
static inline real_t getMatrizDIA(const real_t *A, int i, int j, int n) {
    int offset = j - i; // Deslocamento da diagonal
    // Verifica se está dentro das 7 diagonais (-3 <= offset <= 3)
    if (offset >= -OFFSET_CENTER && offset <= OFFSET_CENTER) {
        int diag_idx = offset + OFFSET_CENTER; // Índice da coluna (0 a 6)
        // O elemento A[i][diag_idx] no armazenamento de matriz
        return A[i * N_DIAG + diag_idx]; 
    }
    return 0.0;
}

//funções de alocação e liberação de memória

void imprimeSistema(int n, real_t *A, real_t *B) {
    printf("--- Matriz A ---\n");
    // Aqui N_DIAG deveria ser k (número de diagonais da matriz A)
    int k_temp = 7; // Assumindo k=7 para Matriz A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Use uma função para acessar o valor A[i][j] no formato DIA
            printf("%5.1f ", getMatrizDIA(A, i, j, n));
        }
        printf("\n");
    }
}

// Imprime a matriz de coeficientes e o vetor B, mostrando apenas as 7 diagonais armazenadas.
// A é a matriz no formato DIA (n x 7)
void imprimeDiagonais(int n, real_t *A, real_t *B) {
    printf("\n--- Matriz de Coeficientes no Formato DIA (%d Diagonais) ---\n", N_DIAG);
    
    // Imprime um cabeçalho para indicar quais diagonais estão sendo exibidas
    printf("Diagonais (Offset j-i):\n");
    printf(" [ -3 ] [ -2 ] [ -1 ] [ 0  ] [ +1 ] [ +2 ] [ +3 ] | Vetor B\n");
    
    // Itera sobre as linhas da matriz (i de 0 a n-1)
    for (int i = 0; i < n; ++i) {
        // Itera sobre as 7 colunas do formato DIA (diag_idx de 0 a 6)
        for (int diag_idx = 0; diag_idx < N_DIAG; ++diag_idx) {
            // O elemento é A[i][diag_idx], acessado em 1D como A[i * N_DIAG + diag_idx]
            printf("%6.1f ", A[i * N_DIAG + diag_idx]);
        }
        // Imprime o elemento b[i] correspondente
        printf(" | %6.1f\n", B[i]);
    }
}


// Para: void criaKDiagonal(int n, int k, MatrizDIA A, real_t *B) {
// ...
void criaKDiagonal(int n, int k, real_t *A, real_t *B) {
    // d é o raio de diagonais
    // (k é o número total de diagonais, ímpar: k = 2d + 1)
    int d = (k - 1) / 2;
    int diag_offset_center = d; // Índice da diagonal principal na matriz de armazenamento (d)

    // preenche A no formato DIA (n x k)
    // O código original aloca n*n, mas o loop preenche n*n.
    // **Assumindo que A foi realocado para n*k em cgSolver.c**
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int offset = j - i;
            // verifica se o elemento está dentro da banda diagonal: |i - j| <= d
            if (abs(offset) <= d) {
                int diag_idx = offset + diag_offset_center;
                // Preenche com valor aleatório no formato DIA: A[i * k + diag_idx]
                A[i * k + diag_idx] = generateRandomA(i, j, k);
            }
            // Não precisa de 'else' pois a memória deve ser inicializada com zero 
            // no cgSolver.c se n*k for alocado (ou inicializado com 0.0, mas 
            // zeros não são armazenados na representação DIA. O código vai depender
            // de inicialização com calloc e uso do loop acima).
        }
    }
    
    //preenche o vetor
    for (int i = 0; i < n; ++i) {
        B[i] = generateRandomB(k);
    }
}

// **Você precisará criar uma versão da getMatrizDIA que aceite k e d_A.**
static inline real_t MatrizDiagonal(const real_t *A, int i, int j, int n, int k, int d) {
    int offset = j - i;
    if (abs(offset) <= d) {
        int diag_idx = offset + d; 
        return A[i * k + diag_idx]; 
    }
    return 0.0;
}
// sislin.c

// Novo formato de A: MatrizDIA (n x k)
// Novo formato de ASP: MatrizDIA (n x 7)
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t *ASP, real_t *bsp, real_t *tempo)
{
    *tempo = timestamp();

    // d_A é o raio de banda de A
    int d_A = (k - 1) / 2;
    int d_ASP = 3; // Raio de banda de ASP (assumindo 7 diagonais: (7-1)/2)
    int k_ASP = N_DIAG; // 7
    int offset_ASP = OFFSET_CENTER; // 3

    // 1. Calcula A' = A~T * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Apenas calcula e armazena se (i, j) está dentro da banda de ASP (7 diagonais)
            if (abs(i - j) <= d_ASP) {
                real_t soma = 0.0;
                // (AT A)i,j = Σ_k A[k, i] * A[k, j]
                // Os limites de k2 podem ser reduzidos (min/max), mas o loop completo funciona
                for (int k2 = 0; k2 < n; ++k2) {
                    // Substitua o acesso A[k2 * n + i] * A[k2 * n + j] 
                    // por acessos via função DIA:
                    real_t val_ki = MatrizDiagonal(A, k2, i, n, k, d_A); // A[k2][i]
                    real_t val_kj = MatrizDiagonal(A, k2, j, n, k, d_A); // A[k2][j]
                    soma += val_ki * val_kj;
                }
                // Armazenamento no formato DIA (n x 7)
                int offset_diag = j - i;
                int diag_idx = offset_diag + offset_ASP;
                ASP[i * k_ASP + diag_idx] = soma;
            }
        }
    }

    // 2. Calcula b' = AT * b
    for (int i = 0; i < n; ++i) {
        real_t soma = 0.0;
        // (AT b)i = Σ_k A[k, i] * b[k]
        for (int k2 = 0; k2 < n; ++k2) {
            // Substitua o acesso A[k2 * n + i] * b[k2]
            // por acesso DIA para A:
            real_t val_ki = MatrizDiagonal(A, k2, i, n, k, d_A); // A[k2][i]
            soma += val_ki * b[k2];
        }
        bsp[i] = soma;
    }

    *tempo = timestamp() - *tempo;
}


//funç~es de pré-condicionamento

//Gera os vetores de Decomposição DLU para uma matriz A.
//D: Diagonal principal
//L: Sub-diagonal
//U: Super-diagonal
// sislin.c

// sislin.c

//Gera os vetores de Decomposição DLU para uma matriz A (no formato DIA: n x k).
//D: Diagonal principal
//L: Sub-diagonal
//U: Super-diagonal
void geraDLU (real_t *A, int n, int k,
              real_t *D, real_t *L, real_t *U, rtime_t *tempo, double eps)
{
    *tempo = timestamp();

    // Raio de banda: (k-1)/2. Este valor é o índice da diagonal principal.
    int d_A = (k - 1) / 2; 
    
    // Índices de coluna no armazenamento DIA (n x k) para as 3 diagonais centrais:
    int diag_D_idx = d_A;     // Índice da Diagonal Principal (offset 0)
    int diag_L_idx = d_A - 1; // Índice da 1ª Subdiagonal (offset -1)
    int diag_U_idx = d_A + 1; // Índice da 1ª Superdiagonal (offset +1)

    // A matriz A é armazenada em formato DIA (n x k). 
    // O acesso é dado por A[linha * k + indice_diagonal].

    /* 1. Preenche o vetor D com a diagonal principal */
    for (int i = 0; i < n; ++i) {
        // Acesso corrigido no formato DIA: D[i] = A[i][diag_D_idx]
        D[i] = A[i * k + diag_D_idx]; 
    }
    
    /* 2. Preenche os vetores L (Subdiagonal) e U (Superdiagonal) */
    // Note que L e U armazenam n-1 elementos.
    for (int i = 0; i < n - 1; ++i) {
        /* Subdiagonal L: L[i] recebe A[i+1, i] (offset -1) */
        // O elemento A[i+1, i] está na linha i+1, na coluna diag_L_idx do array DIA.
        L[i] = A[(i + 1) * k + diag_L_idx];
        
        /* Superdiagonal U: U[i] recebe A[i, i+1] (offset +1) */
        // O elemento A[i, i+1] está na linha i, na coluna diag_U_idx do array DIA.
        U[i] = A[i * k + diag_U_idx];
    }

    /* 3. Proteção: se alguma diagonal for zero ou muito pequena */
    for (int i = 0; i < n; ++i) {
        if (ABS(D[i]) < eps) {
            D[i] = eps; // Define um valor mínimo (eps) para evitar instabilidade
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
        printf("pré-condicionador com w=%lf não implementado (0.0).\n", w);
        exit(1);
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
    exit(1);
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
    int d_A = (k - 1) / 2; // Raio de banda da matriz original A
    for (int i = 0; i < n; ++i) {
        real_t Ax_i = 0.0;
        
        // Loop otimizado sobre as diagonais de A (k diagonais)
        for (int diag_idx = 0; diag_idx < k; diag_idx++) {
            int offset = diag_idx - d_A; // -d_A, ..., d_A
            int j = i + offset; 
            
            // Checa limites de j: 0 <= j < n
            if (j >= 0 && j < n) {
                Ax_i += A[i * k + diag_idx] * X[j];
            }
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



