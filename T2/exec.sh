#!/bin/bash

# ==============================================================================
# Script de Benchmark para o Trabalho 2 - HPC
# Autores: Brenda Pinheiro Ricci & João Pedro Vieira Santos
# ==============================================================================

# --- Configurações Iniciais ---
EXECUTABLE="./cgSolverOld"
RESULT_DIR="resultadosT1"
LIKWID_CMD="/home/soft/likwid/bin/likwid-perfctr"

# Parâmetros fixos do enunciado
K=7
MAXIT=25
OMEGA=0.0
EPSILON=1.0e-9

# Tamanhos de N solicitados no enunciado
SIZES=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# Grupos do LIKWID a serem testados
# Ajuste conforme a arquitetura (ex: MEM ou L3, L2CACHE ou CACHE)
# Verifique disponibilidade com 'likwid-perfctr -a'
GROUP_FLOPS="FLOPS_DP"
GROUP_L2="L2CACHE" 
GROUP_MEM="L3" # Ou "MEM" dependendo da arquitetura
# Nota: FLOPS_AVX pode ser adicionado se disponível na arquitetura

# --- Verificações ---

if ! command -v $LIKWID_CMD &> /dev/null; then
    echo "ERRO: LIKWID não encontrado. Execute em uma máquina com LIKWID instalado."
    exit 1
fi

# Cria diretórios para armazenar os logs brutos
mkdir -p $RESULT_DIR/tempo
mkdir -p $RESULT_DIR/flops
mkdir -p $RESULT_DIR/l2
mkdir -p $RESULT_DIR/mem

# --- 1. Compilação ---
echo "--- 1. Compilando (Flags de Otimização + LIKWID) ---"

# Limpa anterior
rm -f $EXECUTABLE

# Compilação conforme enunciado: -O3 -march=native -mavx -fopt-info-vec
# Adiciona -DLIKWID_PERFMON para ativar os marcadores no código C
make

if [ $? -ne 0 ]; then
    echo "ERRO: Falha na compilação."
    exit 1
fi
echo "Compilação SUCESSO."

# --- 2. Execução dos Testes ---

echo ""
echo "--- 2. Iniciando Benchmarks ---"
echo "Parâmetros: k=$K, maxit=$MAXIT, omega=$OMEGA"
echo "Logs serão salvos em ./$RESULT_DIR/"

for n in "${SIZES[@]}"; do
    echo "========================================"
    echo "Executando para N = $n"
    
    # Prepara a entrada para o programa (stdin)
    # Formato: N OMEGA MAXIT EPSILON
    INPUT_STR="$n $OMEGA $MAXIT $EPSILON"

    # ---------------------------------------------------------
    # A. Tempo de Execução (Medido internamente pelo utils.h)
    # Executamos sem overhead do likwid (apenas wrapper básico se necessário)
    # ---------------------------------------------------------
    echo "  [1/4] Medindo Tempo..."
    # O comando 'echo' alimenta o stdin do executável
    echo "$INPUT_STR" | $EXECUTABLE > "$RESULT_DIR/tempo/N_${n}.log"

    # ---------------------------------------------------------
    # B. Operações Aritméticas (FLOPS_DP)
    # ---------------------------------------------------------
    echo "  [2/4] Medindo FLOPS ($GROUP_FLOPS)..."
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

# --- 3. Consolidação (Exemplo Simples) ---

echo ""
echo "--- 3. Resumo da Execução ---"
echo "Todos os testes concluídos."
echo "Para gerar os gráficos, extraia os dados dos arquivos em $RESULT_DIR/."
echo ""
echo "Exemplo de extração rápida de TEMPO (op1 - Iteração):"
echo "N | Tempo Iteração (s)"
echo "--|-------------------"
for n in "${SIZES[@]}"; do
    # O script assume que o tempo da iteração é a PENÚLTIMA linha numérica impressa pelo seu código main
    # Formato do seu printf final:
    # ...
    # tPrecond
    # tempoIter  <-- Pegamos este
    # tResiduo
    
    TIME_VAL=$(grep -v "Nao calculado" "$RESULT_DIR/tempo/N_${n}.log" | tail -n 2 | head -n 1)
    echo "$n | $TIME_VAL"
done

echo ""
echo "Script finalizado."