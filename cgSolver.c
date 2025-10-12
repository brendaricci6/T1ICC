#include "sislin.h"
#include "pcgc.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    int n;      //demenção do SL >10 
    int k;      //numero de diagonais da matriz >1 e impar
    double omega; //pre-condicionador
    int maxit;  //num máx de iterações
    double epsilon; //erro aprox absoluto máximo
    rtime_t t_pc, t_iter_total, t_iter_medio, t_residuo;
    real_t norma_erro, norma_residuo;

    //leitura dos valores 
    int items_read = scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon);

    //verifica se todos os valores foram lidos corretamente
    if (items_read < 5) {
        fprintf(stderr, "Erro: Não foi possível ler todos os 5 valores de entrada.\n");
        return 1; //retorna erro
    }

    //verifica dimensão > 10
    if( n <= 10 ){
        fprintf(stderr, "valor de dimensão inválido\n"); 
        return 1; 
    }

    //verifica diagonais > 1 e impar
    if (( k <= 1) || (k % 2 == 0)){
        fprintf(stderr, "valor de diagonais inválido\n"); 
        return 1; 
    }

    printf("\n--- Valores lidos ---\n");
    printf("Dimensão (n): %d\n", n);
    printf("Diagonais (k): %d\n", k);
    printf("Pré-condicionador (ω): %f\n", omega);
    printf("Max. Iterações (maxit): %d\n", maxit);
    printf("Tolerância (ε): %g\n", epsilon);
    printf("---------------------\n");
    
    //-------------------- Inicialização -----------------------
    srandom(20252); 

    //-------------------- Gera Sistema linear original -----------------------
    real_t *matrizA = NULL;
    real_t *vetorB = NULL;
    real_t *x = calloc(n, sizeof(real_t)); //vetor de respostas

    //printf("Criando um sistema aleatório %dx%d com uma matriz %d-diagonal.\n\n", n, n, k);

    //aloca matriz e vetor
    matrizA = malloc(n * n * sizeof(real_t));
    if (!matrizA) {
        fprintf(stderr, "Erro de alocação de memória para matrizA\n");
    return 1;
    }
    vetorB = malloc(n * sizeof(real_t));
    if (!vetorB) {
        fprintf(stderr, "Erro de alocação de memória para vetorB\n");
        free (matrizA);
    return 1;
    }

    // Gerando sistema Linear
    printf("Gerando sistema tridiagonal simétrico positivo...\n");
    rtime_t t_total = timestamp();
    criaKDiagonal(n, k, matrizA, vetorB);
    t_total = timestamp() - t_total;
    printf("Sistema gerado em %.6lfs.\n\n", t_total);

    //-------------------- tranf p sistema simétrico positivo -----------------------
    real_t *ASP = NULL;
    real_t *bsp = NULL;

    //aloca matriz e vetor
    ASP = malloc(n * n * sizeof(real_t));
    if (!ASP) {
        fprintf(stderr, "Erro de alocação de memória para ASP\n");
    return 1;
    }
    bsp = malloc(n * sizeof(real_t));
    if (!bsp) {
        fprintf(stderr, "Erro de alocação de memória para vetorB\n");
        free (ASP);
    return 1;
    }
    rtime_t algum_tempo;
    genSimetricaPositiva(matrizA, vetorB, n, k, ASP, bsp, &algum_tempo);

    free(matrizA);
    free(vetorB);

    //-------------------- prepara pré-cond -----------------------
    // Decomposicao DLU
    real_t *D;
    real_t *L;
    real_t *U;
    rtime_t tDLU;

    //rtime_t t_pc = timestamp(); 
    geraDLU(ASP, n, k, &D, &L, &U, &tDLU);
    printf("Decomposição DLU gerada em %.6lfs.\n", tDLU);


    // Gerando Pre-Condicionador
    real_t *matrizPreCond;
    rtime_t tPrecond;
    geraPreCond(D, L, U, omega, n, k, &matrizPreCond, &tPrecond);
    printf("Pré-condicionador gerado em %.6lfs.\n", tPrecond);

    //-------------------- sol com gradientes conjugados -----------------------
    printf("\nExecutando método de Gradientes Conjugados...\n");
    int iter = 0;
    real_t residuo = 0.0;

    rtime_t tCG = timestamp();
    // função a ser implementada
    iter = conjugateGradient(ASP, bsp, x, n, k, epsilon, maxit, &residuo, matrizPreCond);  
    tCG = timestamp() - tCG;

    printf("Gradientes Conjugados concluído em %.6lfs (%d iterações)\n", tCG, iter);
    printf("Resíduo final estimado: %.6e\n", residuo);

    // Residuo

    rtime_t tResiduo;
    norma_residuo = calcResiduoSL(ASP, bsp, x, n, k, &tResiduo);
    printf("Norma do resíduo calculada em %.6lfs: %.6e\n", tResiduo, norma_residuo);

    //-------------------- impressão dos resultados -----------------------

    t_pc = tPrecond;                     
    t_iter_total = tCG;                  
    t_iter_medio = (iter > 0) ? (tCG / iter) : 0.0; 
    t_residuo = tResiduo; 

    printf ("\n \n ################################## IMPRESSÃO DOS RESULTADOS ################################## \n \n");
        // Imprime a dimensão do sistema
    printf("%d\n", n);

    // Imprime o vetor solução 'x' com 16 dígitos de precisão [cite: 71]
    for (int i = 0; i < n; ++i) {
        printf("%.16g ", x[i]);
    }
    printf("\n");

    // Imprime a norma do erro, resíduo e os tempos com 8 dígitos de precisão [cite: 72]
    
    // norma: norma máxima do erro aproximado em x 
    printf("%.8g\n", norma_erro);
    
    // residuo: A norma euclidiana do resíduo 
    printf("%.8g\n", norma_residuo);
    
    // Tempo PC: tempo para calcular a matriz pré-condicionante 
    printf("%.8g\n", t_pc);
    
    // Tempo iter: Tempo médio para calcular uma iteração 
    printf("%.8g\n", t_iter_medio);
    
    // Tempo residuo: Tempo para calcular a norma euclidiana do resíduo 
    printf("%.8g\n", t_residuo);

    //-------------------- limpa as alocações -----------------------
    free(ASP);
    free(bsp);
    free(x);
    free(D);
    free(L);
    free(U);
    free(matrizPreCond);

    return 0;
}