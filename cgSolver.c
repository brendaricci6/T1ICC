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
    //alorpar memória para a matriz A e o vetor b 
    //gerar o sistema linear inicial (usando criaKDiagonal) 
    //gerar o pré condionador (usando gera LU e leraPreCond)

    //implementar criaKDiagonal e chamr aqui 


    return 0; 
}