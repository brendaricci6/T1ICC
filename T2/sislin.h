#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdlib.h>
#include "utils.h"

// Configuração do LIKWID (obrigatório pelo enunciado)
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

// Definições do Sistema
#define N_DIAG 7        // Total de diagonais
#define OFFSET_CENTER 3 // Onde fica a diagonal principal (índice 3)

typedef double real_t;
typedef double rtime_t;

// Funções do Sislin
void criaKDiagonal(int n, int k, real_t *A, real_t *b);
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, real_t *ASP, real_t *bsp, rtime_t *tempo);
void geraDLU(real_t *A, int n, int k, real_t *D, real_t *L, real_t *U, rtime_t *tempo, double eps);
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t *M, rtime_t *tempo, double eps);

// OP2: Cálculo do Resíduo
real_t calcResiduoSL(real_t *A, real_t *b, real_t *X, int n, int k, rtime_t *tempo);

// Debug
void imprimeSistema(int n, real_t *A, real_t *B);
void imprimeDiagonais(int n, real_t *A, real_t *B);

#endif