#!/bin/bash

# ==============================================================================
# Script de Benchmark para o Trabalho 2 - HPC
# Atualizado para compatibilidade com cgSolver.c
# ==============================================================================

# --- Configurações Iniciais ---
EXECUTABLE="./cgSolver"
RESULT_DIR="resultadosT1"
LIKWID_CMD="/home/soft/likwid/bin/likwid-perfctr"

# Parâmetros fixos do enunciado
K=7              # Número de diagonais (Sempre 7 conforme solicitado)
MAXIT=25         # Máximo de iterações
OMEGA=0.0        # Pré-condicionador (0.0 ativa cálculo normal, -1.0 desativa)
EPSILON=1.0e-9   # Tolerância

# Tamanhos de N solicitados
SIZES=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# Grupos do LIKWID a serem testados
# DICA: Verifique a disponibilidade com 'likwid-perfctr -a' na máquina alvo
GROUP_FLOPS="FLOPS_DP"
GROUP_L2="L2CACHE" 
GROUP_MEM="L3" # Atenção: Em algumas máquinas pode ser "MEM" ou "L3"

# --- Verificações ---

if ! command -v $LIKWID_CMD &> /dev/null; then
    echo "AVISO: LIKWID não encontrado no PATH."
    echo "O script tentará rodar, mas os passos com likwid falharão se não estiver instalado."
    # Não vamos sair com exit 1 para permitir testes de tempo sem likwid se necessário
fi

# Cria diretórios para armazenar os logs brutos
mkdir -p $RESULT_DIR/tempo
mkdir -p $RESULT_DIR/flops
mkdir -p $RESULT_DIR/l2
mkdir -p $RESULT_DIR/mem

# --- 1. Compilação ---
echo "--- 1. Compilando (Flags de Otimização + LIKWID) ---"

# Limpa executável anterior
rm -f $EXECUTABLE

# O script assume que existe um Makefile configurado corretamente
# para gerar o binário com nome 'cgSolver'
make

if [ $? -ne 0 ]; then
    echo "ERRO: Falha na compilação. Verifique seu Makefile."
    exit 1
fi
echo "Compilação SUCESSO."

# --- 2. Execução dos Testes ---

echo ""
echo "--- 2. Iniciando Benchmarks ---"
echo "Parâmetros Fixos: k=$K, omega=$OMEGA, maxit=$MAXIT, epsilon=$EPSILON"
echo "Logs serão salvos em ./$RESULT_DIR/"

for n in "${SIZES[@]}"; do
    echo "========================================"
    echo "Executando para N = $n"
    
    # CORREÇÃO PRINCIPAL AQUI:
    # O C espera: scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon);
    INPUT_STR="$n $K $OMEGA $MAXIT $EPSILON"

    # ---------------------------------------------------------
    # A. Tempo de Execução (Sem Likwid, apenas medição interna)
    # ---------------------------------------------------------
    echo "  [1/4] Medindo Tempo..."
    echo "$INPUT_STR" | $EXECUTABLE > "$RESULT_DIR/tempo/N_${n}.log"

    # ---------------------------------------------------------
    # B. Operações Aritméticas (FLOPS_DP)
    # ---------------------------------------------------------
    echo "  [2/4] Medindo FLOPS ($GROUP_FLOPS)..."
    # O comando likwid recebe o executável e passa o INPUT_STR via pipe
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_FLOPS -m $EXECUTABLE > "$RESULT_DIR/flops/N_${n}.log" 2>&1

    # ---------------------------------------------------------
    # C. Cache Miss L2
    # ---------------------------------------------------------
    echo "  [3/4] Medindo Cache L2 ($GROUP_L2)..."
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_L2 -m $EXECUTABLE > "$RESULT_DIR/l2/N_${n}.log" 2>&1

    # ---------------------------------------------------------
    # D. Banda de Memória
    # ---------------------------------------------------------
    echo "  [4/4] Medindo Memória ($GROUP_MEM)..."
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_MEM -m $EXECUTABLE > "$RESULT_DIR/mem/N_${n}.log" 2>&1

done

# --- 3. Consolidação Rápida ---

echo ""
echo "--- 3. Resumo da Execução ---"
echo "Todos os testes concluídos."
echo ""
echo "Resumo dos tempos de iteração (extraído dos logs):"
echo "N     | Tempo Iteração (s)"
echo "------|-------------------"
for n in "${SIZES[@]}"; do
    # O script assume que o tempo da iteração é o penúltimo valor numérico 
    # impresso pelo printf final do código C.
    # Ordem do printf no C: tPrecond -> tempoIter -> tResiduo
    
    # Verifica se o arquivo existe antes de tentar ler
    if [ -f "$RESULT_DIR/tempo/N_${n}.log" ]; then
        TIME_VAL=$(grep -v "Nao calculado" "$RESULT_DIR/tempo/N_${n}.log" | tail -n 2 | head -n 1)
        printf "%-5d | %s\n" "$n" "$TIME_VAL"
    else
        printf "%-5d | Erro/Não encontrado\n" "$n"
    fi
done

echo ""
echo "Script finalizado."
