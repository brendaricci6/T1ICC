#include "sislin.h"
#include "pcgc.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    //inicializa LIKIWD se definido
    LIKWID_MARKER_INIT;

    int n;          // dimensão do SL >10
    int k = 7;      // número de diagonais da matriz >1 e ímpar
    double omega;   // pré-condicionador
    int maxit;      // número máx. de iterações
    double epsilon; // erro aprox. absoluto máximo

    //variáveis para armazenar tempos de execução
    rtime_t tDLU = 0.0, tPrecond = 0.0, tempoIter = 0.0, tResiduo = 0.0;
    //variáveis para armazenar normas
    real_t normaFinal = 0.0, norma_residuo = 0.0;

    // ========== Leitura da entrada ============

    //lê n, k, omega, maxit, epsilon da entrada padrão (STDIN)
    int items_read = scanf("%d %lf %d %lf", &n, &omega, &maxit, &epsilon);

    //verifica se a leitura foi bem-sucedida
    if (items_read < 4) {
        printf("Erro: Não foi possível ler todos os 4 valores de entrada.\n");
        return 1;
    }

    //validação dos parâmetros
    if (n <= 10) {
        printf("Erro: dimensão deve ser > 10\n");
        return 1;
    }

    // =========== Geração do sistema ==========
    
    // Aloca matriz A inicializando com 0 (Layout V2: n*k)
    real_t *A = calloc(n * k, sizeof(real_t)); 
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

    // Aloca ASP e bsp para o sistema simétrico positivo
    real_t *ASP = calloc(n * k, sizeof(real_t));
    real_t *bsp = calloc(n, sizeof(real_t));

    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tGen);

    //calcula o tempo gasto 
    tGen = timestamp() - tGen;

    // ========== Decomposição DLU ===========

    real_t *D = malloc(n * sizeof(real_t));
    real_t *L = malloc(n * sizeof(real_t));
    real_t *U = malloc(n * sizeof(real_t));

    //função que gera o DLU 
    //calcula a decomposição DLU de A
    //armazena o tempo em tDLU
    geraDLU(ASP, n, k, D, L, U, &tDLU, epsilon);

    // ========== Geração do pré condicionador ===========
    
    //aloca vetor M que armazena o pré-condicionados
    real_t *M;
    //gera o pré-condicionador M usando D, L, U
    //o parâmetro omega e armazena o tempo em tPrecond
    if (omega != -1.0) {
        M = malloc(n * sizeof(real_t));
        if (!M) {
            printf("Erro de alocacao M\n");
            free(A); free(b); free(x); free(D); free(L); free(U); free(ASP); free(bsp);
            return 1;
        }
        geraPreCond(D, L, U, omega, n, k, M, &tPrecond, epsilon);
    } else {
        M = NULL;
    }

    // ========== Execução do método PCG ===========
    
    //iterações 
    int iter = 0;

    //mede o tempo de execução do GCG
    //executa o pcg 
    iter = gradienteConjugado(ASP, bsp, x, n, maxit, epsilon, M, &normaFinal, &tempoIter);
    
    // =========== Clculo do resíduo ==========

    //calcula a norma resíduo com os valores de A e x obtidos 
    norma_residuo = calcResiduoSL(A, b, x, n, k, &tResiduo);
    
    // ========== Impressão dos resultados ===========
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

    // ============Libera memória ==========
    free(A);
    free(b);
    free(x);
    free(ASP);
    free(bsp);
    free(D);
    free(L);
    free(U);
    if (M) free(M); 

    LIKWID_MARKER_CLOSE;

    return 0;
}