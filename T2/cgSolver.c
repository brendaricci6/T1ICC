#include "sislin.h" //funções para sistemas lineares
#include "pcgc.h" //funções específicas do método PCG
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    int n;          // dimensão do SL >10
    int k = 7;          // número de diagonais da matriz >1 e ímpar
    double omega;   // pré-condicionador
    int maxit;      // número máx. de iterações
    double epsilon; // erro aprox. absoluto máximo

    //variáveis para armazenar tempos de execução
    rtime_t tDLU = 0.0, tPrecond = 0.0, tempoIter = 0.0, tResiduo = 0.0;
    //variáveis para armazenar normas
    real_t normaFinal = 0.0, norma_residuo = 0.0;

    // ========== Leitura da entrada ==========

    //lê n, k, omega, maxit, epsilon da entrada padrão (STDIN)
    int items_read = scanf("%d %lf %d %lf", &n, &omega, &maxit, &epsilon);

    //verifica se a leitura foi bem-sucedida
    if (items_read < 4) {
        printf("Erro: Não foi possível ler todos os 5 valores de entrada.\n");
        return 1;
    }

    //validação dos parâmetros
    if (n <= 10) {
        printf("Erro: dimensão deve ser > 10\n");
        return 1;
    }

    #ifdef DEBUG
    // exibição dos parâmetros lidos
    printf("\n--- Valores lidos ---\n");
    printf("Dimensão (n): %d\n", n);
    printf("Diagonais (k): %d\n", k);
    printf("Pré-condicionador (ω): %f\n", omega);
    printf("Max. Iterações (maxit): %d\n", maxit);
    printf("Tolerância (ε): %g\n", epsilon);
    printf("==========-\n");
    #endif

    // ========== Geração do sistema ==========
    #ifdef DEBUG
    printf("Gerando sistema tridiagonal simétrico positivo...\n");
    #endif

    real_t *A = calloc(n * k, sizeof(real_t)); //aloca matriz A inicializando com 0
    real_t *b = calloc(n, sizeof(real_t)); //aloca vetor B inicializando com 0
    real_t *x = calloc(n, sizeof(real_t)); //aloca o vetor de solução x inicializando com 0 
    
    //verifica a alocação de memória 
    if (!b || !x || !A) {
        printf("Erro de alocação de memória em b ou x ou a\n");
        return 1;
    }

    //marca tempo de geração da matriz A e vetor b 
    rtime_t tGen = timestamp();

    //chama função que cria a matriz e o vetor B 
    criaKDiagonal(n, k, A, b);

    #ifdef DEBUG
    imprimeSistema(n, A, b);
    #endif

    #ifdef DEBUG
    imprimeDiagonais(n, A, b);
    #endif

    
    real_t *ASP = calloc(n * k, sizeof(real_t));
    real_t *bsp = calloc(n, sizeof(real_t));

    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tGen);

    #ifdef DEBUG
    printf("--- Matriz ASP Gerada ---\n");
    imprimeSistema(n, ASP, bsp); // <--- Nova chamada (recomendada)
    imprimeDiagonais(n, ASP, bsp);
    #endif

    //calcula o tempo gasto 
    tGen = timestamp() - tGen;
    #ifdef DEBUG
    printf("Sistema gerado em %.6es.\n\n", tGen);
    #endif


    // ========== Decomposição DLU ==========

    real_t *D = malloc(n * sizeof(real_t));
    real_t *L = malloc((k-1)/2 * n * sizeof(real_t));
    real_t *U = malloc((k-1)/2 * n * sizeof(real_t));

    //função que gera o DLU 
    //calcula a decomposição DLU de A
    //armazena o tempo em tDLU
    geraDLU(ASP, n, k, D, L, U, &tDLU, epsilon);

    #ifdef DEBUG
    printf("Decomposição DLU gerada em %.6es.\n", tDLU);

    printf("U->");
    for (int i = 0; i < (k-1)/2 * n; i++) 
        printf(" %6e", U[i]);
    printf("\n");

    printf("D->");
    for (int i = 0; i < n; i++) 
        printf(" %6e", D[i]);
    printf("\n");

    printf("L->");
    for (int i = 0; i < (k-1)/2 * n; i++) 
        printf(" %6e", L[i]);
    printf("\n");
    #endif
    // ========== Geração do pré-condicionador ==========
    
    //aloca vetor M que armazena o pré-condicionados
    real_t *M;
    //gera o pré-condicionador M usando D, L, U
    //o parâmetro omega e armazena o tempo em tPrecond
    if (omega != -1.0) {
        M = malloc(n * sizeof(real_t));
        if (!M) {
            printf("Erro de alocacao M\n");
            free(A); free(b); free(x); free(D); free(L); free(U);
            return 1;
        }
        geraPreCond(D, L, U, omega, n, k, M, &tPrecond, epsilon);
    } else {
        M = NULL;
    }

    #ifdef DEBUG
    printf("Pré-condicionador gerado em %.6es.\n", tPrecond);ac

    printf("M->");
    for (int i = 0; i < n; i++) 
        printf(" %6e", M[i]);
    printf("\n");
    
    // ========== Execução do método PCG ==========
    printf("\nExecutando método de Gradientes Conjugados Pré-Condicionado...\n");
    #endif

    //iterações 
    int iter = 0;
    real_t res_estimado = 0.0;

    //mede o tempo de execução do GCG
    //executa o pcg 
    //A: matriz
    //b: lado direito
    //x: solução (saída)
    // n: dimensão,
    //maxit: max. iterações
    //epsilon: tolerância
    //M: pré-condicionador
    iter = gradienteConjugado(ASP, bsp, x, n, maxit, epsilon, M, &normaFinal, &tempoIter);
    
    //calcula tempo total da execução do pcg
    #ifdef DEBUG
    printf("PCG concluído em %.6es (%d iterações)\n", tempoIter, iter);

    printf("X->");
    for (int i = 0; i < n; i++) 
        printf(" %6e", x[i]);
    printf("\n");
    #endif


    // ========== Cálculo do resíduo ==========

    //calcula a norma resíduo com os valores de A e x obtidos 
    #ifdef DEBUG
    printf("Calculando residuo...\n");
    #endif
    norma_residuo = calcResiduoSL(A, b, x, n, k, &tResiduo);
    #ifdef DEBUG
    printf("Norma do resíduo calculada em %.6es: %.6e\n", tResiduo, norma_residuo);
    #endif
    
    // ========== Impressão dos resultados ==========
    printf("%d\n", n);
    for (int i = 0; i < n; ++i)
        printf("%.16g ", x[i]);
    printf("\n");

    printf("%.8g\n", normaFinal);
    printf("%.16g\n", norma_residuo);
    tPrecond == 0.0 ? printf("Nao calculado\n") : printf("%.8g\n", tPrecond);
    printf("%.8g\n", tempoIter);
    printf("%.8g\n", tResiduo);
    // printf("Iterações: %d\n", iter);

    // ========== Libera memória ==========
    free(A);
    free(b);
    free(x);
    free(D);
    free(L);
    free(U);
    if (M) free(M); 

    return 0;
}
