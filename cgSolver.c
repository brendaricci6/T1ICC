#include "sislin.h" //funções para sistemas lineares
#include "pcgc.h" //funções específicas do método PCG
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    int n;          // dimensão do SL >10
    int k;          // número de diagonais da matriz >1 e ímpar
    double omega;   // pré-condicionador
    int maxit;      // número máx. de iterações
    double epsilon; // erro aprox. absoluto máximo

    //variáveis para armazenar tempos de execução
    rtime_t tDLU, tPrecond, tPCG, tResiduo;
    //variáveis para armazenar normas
    real_t normaFinal = 0.0, norma_residuo = 0.0;

    // -------------------- Leitura da entrada --------------------

    //lê n, k, omega, maxit, epsilon da entrada padrão (STDIN)
    int items_read = scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon);

    //verifica se a leitura foi bem-sucedida
    if (items_read < 5) {
        fprintf(stderr, "Erro: Não foi possível ler todos os 5 valores de entrada.\n");
        return 1;
    }

    //validação dos parâmetros
    if (n <= 10) {
        fprintf(stderr, "Erro: dimensão deve ser > 10\n");
        return 1;
    }

    if (k <= 1 || k % 2 == 0) {
        fprintf(stderr, "Erro: número de diagonais inválido (deve ser ímpar > 1)\n");
        return 1;
    }

    //exibição dos parâmetros lidos
    printf("\n--- Valores lidos ---\n");
    printf("Dimensão (n): %d\n", n);
    printf("Diagonais (k): %d\n", k);
    printf("Pré-condicionador (ω): %f\n", omega);
    printf("Max. Iterações (maxit): %d\n", maxit);
    printf("Tolerância (ε): %g\n", epsilon);
    printf("---------------------\n");

    // -------------------- Geração do sistema --------------------
    printf("Gerando sistema tridiagonal simétrico positivo...\n");

    real_t *A = calloc(n * n, sizeof(real_t)); //aloca matriz A inicializando com 0
    real_t *b = malloc(n * sizeof(real_t)); //aloca vetor B inicializando com 0
    real_t *x = calloc(n, sizeof(real_t)); //aloca o vetor de solução x inicializando com 0 
    
    //verifica a alocação de memória 
    if (!b || !x || !A) {
        fprintf(stderr, "Erro de alocação de memória em b ou x ou a\n");
        return 1;
    }

    //marca tempo de geração da matriz A e vetor b 
    rtime_t tGen = timestamp();

    //chama função que cria a matriz e o vetor B 
    criaKDiagonal(n, k, A, b);
    // imprimeSistema(n, A, b);

    //calcula o tempo gasto 
    tGen = timestamp() - tGen;
    printf("Sistema gerado em %.6es.\n\n", tGen);

    // -------------------- Decomposição DLU --------------------

    //D: diagonal principal (n elementos)
    real_t *D =  malloc(n * sizeof(real_t));
    //L: triangular inferior (n-1 elementos)
    real_t *L = malloc((n - 1) * sizeof(real_t));
    //U: triangular superior (n-1 elementos)
    real_t *U = malloc((n - 1) * sizeof(real_t));

    //função que gera o DLU 
    //calcula a decomposição DLU de A
    //armazena o tempo em tDLU
    geraDLU(A, n, k, D, L, U, &tDLU, epsilon);
    printf("Decomposição DLU gerada em %.6es.\n", tDLU);

    // printf("U->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", U[i]);
    // printf("\n");

    // printf("D->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", D[i]);
    // printf("\n");

    // printf("L->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", L[i]);
    // printf("\n");

    // -------------------- Geração do pré-condicionador --------------------
    
    //aloca vetor M que armazena o pré-condicionados
    real_t *M = malloc(n * sizeof(real_t));

    //gera o pré-condicionador M usando D, L, U
    //o parâmetro omega e armazena o tempo em tPrecond
    geraPreCond(D, L, U, omega, n, k, M, &tPrecond, epsilon);
    printf("Pré-condicionador gerado em %.6es.\n", tPrecond);

    // printf("M->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", M[i]);
    // printf("\n");
    
    // -------------------- Execução do método PCG --------------------
    printf("\nExecutando método de Gradientes Conjugados Pré-Condicionado...\n");

    //iterações 
    int iter = 0;
    real_t res_estimado = 0.0;

    //mede o tempo de execução do GCG
    tPCG = timestamp();

    //executa o pcg 
    //A: matriz
    //b: lado direito
    //x: solução (saída)
    // n: dimensão,
    //maxit: max. iterações
    //epsilon: tolerância
    //M: pré-condicionador
    iter = gradienteConjugado(A, b, x, n, maxit, epsilon, M, &normaFinal);
    
    //calcula tempo total da execução do pcg
    tPCG = timestamp() - tPCG;
    printf("PCG concluído em %.6es (%d iterações)\n", tPCG, iter);

    // printf("X->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", x[i]);
    // printf("\n");


    // -------------------- Cálculo do resíduo --------------------

    //calcula a norma resíduo com os valores de A e x obtidos 
    norma_residuo = calcResiduoSL(A, b, x, n, k, &tResiduo);
    printf("Norma do resíduo calculada em %.6es: %.6e\n", tResiduo, norma_residuo);

    // -------------------- Impressão dos resultados --------------------
    printf("\n\n############## RESULTADOS FINAIS ##############\n");
    printf("%d\n", n);
    for (int i = 0; i < n; ++i)
        printf("%.16g ", x[i]);
    printf("\n");

    printf("Norma do erro: %.8e\n", normaFinal);
    printf("Norma do resíduo: %.8e\n", norma_residuo);
    printf("Tempo pré-cond: %.8e\n", tPrecond);
    printf("Tempo total PCG: %.8e\n", tPCG);
    printf("Tempo cálculo resíduo: %.8e\n", tResiduo);
    printf("Iterações: %d\n", iter);

    // -------------------- Libera memória --------------------
    free(A);
    free(b);
    free(x);
    free(D);
    free(L);
    free(U);
    if (M) free(M);

    return 0;
}
