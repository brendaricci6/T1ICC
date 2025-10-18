#include "sislin.h"
#include "pcgc.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    int n;          // dimensão do SL >10
    int k;          // número de diagonais da matriz >1 e ímpar
    double omega;   // pré-condicionador
    int maxit;      // número máx. de iterações
    double epsilon; // erro aprox. absoluto máximo

    rtime_t tDLU, tPrecond, tPCG, tResiduo;
    real_t norma_erro = 0.0, norma_residuo = 0.0;

    // -------------------- Leitura da entrada --------------------
    int items_read = scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon);

    if (items_read < 5) {
        fprintf(stderr, "Erro: Não foi possível ler todos os 5 valores de entrada.\n");
        return 1;
    }

    if (n <= 10) {
        fprintf(stderr, "Erro: dimensão deve ser > 10\n");
        return 1;
    }

    if (k <= 1 || k % 2 == 0) {
        fprintf(stderr, "Erro: número de diagonais inválido (deve ser ímpar > 1)\n");
        return 1;
    }

    printf("\n--- Valores lidos ---\n");
    printf("Dimensão (n): %d\n", n);
    printf("Diagonais (k): %d\n", k);
    printf("Pré-condicionador (ω): %f\n", omega);
    printf("Max. Iterações (maxit): %d\n", maxit);
    printf("Tolerância (ε): %g\n", epsilon);
    printf("---------------------\n");

    // -------------------- Geração do sistema --------------------
    printf("Gerando sistema tridiagonal simétrico positivo...\n");

    real_t *A = calloc(n * n, sizeof(real_t));
    real_t *b = malloc(n * sizeof(real_t));
    real_t *x = calloc(n, sizeof(real_t));
    if (!b || !x) {
        fprintf(stderr, "Erro de alocação de memória em b ou x\n");
        return 1;
    }
    rtime_t tGen = timestamp();
    criaKDiagonal(n, k, A, b);
    // imprimeSistema(n, A, b);
    tGen = timestamp() - tGen;
    printf("Sistema gerado em %.6es.\n\n", tGen);

    // -------------------- Decomposição DLU --------------------
    real_t *D =  malloc(n * sizeof(real_t));
    real_t *L = malloc((n - 1) * sizeof(real_t));
    real_t *U = malloc((n - 1) * sizeof(real_t));
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
    real_t *M = malloc(n * sizeof(real_t));
    geraPreCond(D, L, U, omega, n, k, M, &tPrecond, epsilon);
    printf("Pré-condicionador gerado em %.6es.\n", tPrecond);

    // printf("M->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", M[i]);
    // printf("\n");
    
    // -------------------- Execução do método PCG --------------------
    printf("\nExecutando método de Gradientes Conjugados Pré-Condicionado...\n");

    int iter = 0;
    real_t res_estimado = 0.0;
    tPCG = timestamp();
    iter = gradienteConjugado(A, b, x, n, maxit, epsilon, M);
    tPCG = timestamp() - tPCG;

    printf("PCG concluído em %.6es (%d iterações)\n", tPCG, iter);
    // printf("X->");
    // for (int i = 0; i < n; i++) 
    //     printf(" %6e", x[i]);
    // printf("\n");
    // -------------------- Cálculo do resíduo --------------------
    norma_residuo = calcResiduoSL(A, b, x, n, k, &tResiduo);
    printf("Norma do resíduo calculada em %.6es: %.6e\n", tResiduo, norma_residuo);

    // -------------------- Impressão dos resultados --------------------
    printf("\n\n############## RESULTADOS FINAIS ##############\n");
    printf("%d\n", n);
    for (int i = 0; i < n; ++i)
        printf("%.16g ", x[i]);
    printf("\n");

    printf("Norma do erro: %.8e\n", norma_erro);
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
