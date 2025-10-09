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

    //leitura dos valores 
    int items_read = scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon);

    //verifica se todos os valores foram lidos corretamente
    if (items_read < 5) {
        fprintf(stderr, "Erro: Não foi possível ler todos os 5 valores de entrada.\n");
        return 1; //retorna erro
    }

    //verifica dimensão > 10
    if( n <= 10 ){
        printf("valor de dimensão inválido\n"); 
        return 1; 
    }

    //verifica diagonais > 1 e impar
    if (( k <= 1) || (k % 2 == 0)){
        printf("valor de diagonais inválido\n"); 
        return 1; 
    }

    printf("\n--- Valores lidos ---\n");
    printf("Dimensão (n): %d\n", n);
    printf("Diagonais (k): %d\n", k);
    printf("Pré-condicionador (ω): %f\n", omega);
    printf("Max. Iterações (maxit): %d\n", maxit);
    printf("Tolerância (ε): %g\n", epsilon);
    printf("---------------------\n");
    
    //chamar funções de sislin.c para: 
    //alorcar memória para a matriz A e o vetor b  

    real_t **matrizA = NULL;
    real_t *vetorB = NULL;
    real_t *x = calloc(n, sizeof(real_t));

    printf("Criando um sistema aleatório %dx%d com uma matriz %d-diagonal.\n\n", n, n, k);

    if (!alocaKDiagonal(n, &matrizA)) {
        return 1;
    }
    if (!alocaVetorB(n, &vetorB)) {
        liberaKDiagonal(n, matrizA);
        return 1;
    }

    // Gerando sistema Linear
    printf("Gerando sistema tridiagonal simétrico positivo...\n");
    rtime_t t_total = timestamp();
    criaKDiagonal(matrizA, vetorB, n, k);
    t_total = timestamp() - t_total;
    printf("Sistema gerado em %.6lfs.\n\n", t_total);


    // Decomposicao DLU
    real_t *D;
    real_t *L;
    real_t *U;
    rtime_t tDLU;

    geraDLU(matrizA, n, k, &D, &L, &U, &tDLU);
    printf("Decomposição DLU gerada em %.6lfs.\n", tDLU);


    // Gerando Pre-Condicionador
    real_t *matrizPreCond;
    rtime_t tPrecond;
    geraPreCond(D, L, U, omega, n, k, &matrizPreCond, &tPrecond);
    printf("Pré-condicionador gerado em %.6lfs.\n", tPrecond);

    // Gradientes Conjugados
    printf("\nExecutando método de Gradientes Conjugados...\n");
    int iter = 0;
    real_t residuo = 0.0;

    rtime_t tCG = timestamp();
    // função a ser implementada
    iter = conjugateGradient(matrizA, vetorB, x, n, k, epsilon, maxit, &residuo, matrizPreCond);  
    tCG = timestamp() - tCG;

    printf("Gradientes Conjugados concluído em %.6lfs (%d iterações)\n", tCG, iter);
    printf("Resíduo final estimado: %.6e\n", residuo);

    // Residuo

    rtime_t tResiduo;
    real_t rnorm = calcResiduoSL(matrizA, vetorB, x, n, k, &tResiduo);
    printf("Norma do resíduo calculada em %.6lfs: %.6e\n", tResiduo, rnorm);


    return 0; 
}