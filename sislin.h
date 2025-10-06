#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"

//funções criadas 
int alocaKDiagonal(int n, real_t ***matriz);
int alocaVetorB(int n, real_t **vetorA);
void liberaKDiagonal(int n, real_t **matriz); 
void liberaVetorB(real_t *vetorB);
void imprimeSistema(int n, real_t **A, real_t *B);
void criaKDiagonal(int n, int k, real_t ***A, real_t **B);

void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t **ASP, real_t **bsp, real_t *tempo);
void geraDLU (real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, real_t *tempo);
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, real_t *tempo);
real_t calcResiduoSL (real_t *A, real_t *b, real_t *X, int n, int k, real_t *tempo);


#endif // __SISLIN_H__

