#!/bin/bash

# --- Configurações Básicas ---
EXECUTABLE="./cgSolver"
RESULT_DIR="resultadosT1"
LIKWID_CMD="/home/soft/likwid/bin/likwid-perfctr"

# Parâmetros do Problema (Conforme Enunciado)
K=7              # Número de diagonais
MAXIT=25         # Limite de iterações do Gradiente Conjugado
OMEGA=0.0        # 0.0 = Pré-condicionador padrão / -1.0 = Sem PC
EPSILON=1.0e-9   # Tolerância para convergência

# Definição dos Tamanhos de Sistema (N)
SIZES=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# --- Configurações do LIKWID ---
# Grupos de performance a serem monitorados
GROUP_FLOPS="FLOPS_DP"
GROUP_L2="L2CACHE" 
GROUP_MEM="L3" 

# --- Verificação de Ambiente ---

# Testa se o LIKWID está acessível no caminho especificado
if ! command -v $LIKWID_CMD &> /dev/null; then
    echo "AVISO: Executável do LIKWID não encontrado."
    echo "O script continuará, mas as métricas de hardware falharão."
fi

# Preparação de Diretórios
# Cria a estrutura de pastas para organizar os logs por métrica
mkdir -p $RESULT_DIR/tempo
mkdir -p $RESULT_DIR/flops
mkdir -p $RESULT_DIR/l2
mkdir -p $RESULT_DIR/mem

# --- 1. Compilação do Projeto ---
echo "--- 1. Compilando (Otimizações + Marker API) ---"

# Garante uma compilação limpa
rm -f $EXECUTABLE

# Invoca o Makefile para gerar o binário
make

if [ $? -ne 0 ]; then
    echo "ERRO: Falha crítica na compilação via Makefile."
    exit 1
fi
echo "Compilação concluída com sucesso."

# --- 2. Execução dos Benchmarks ---

echo ""
echo "--- 2. Iniciando Bateria de Testes ---"
echo "Params: K=$K, OMEGA=$OMEGA, MAXIT=$MAXIT"
echo "Saída: Logs salvos em ./$RESULT_DIR/"

for n in "${SIZES[@]}"; do
    echo "========================================"
    echo "Processando N = $n"
    
    # Formatação da string de entrada para o scanf do C
    # Ordem: N K OMEGA MAXIT EPSILON
    INPUT_STR="$n $K $OMEGA $MAXIT $EPSILON"

    # A. Tempo de Execução
    # Execução pura sem overhead do LIKWID para medir 'Wall Clock'
    echo "  [1/4] Medindo Tempo..."
    echo "$INPUT_STR" | $EXECUTABLE > "$RESULT_DIR/tempo/N_${n}.log"

    # B. Operações Aritméticas (DP FLOPS)
    # Usa -C 0 (Core 0) e -m (Marker API)
    echo "  [2/4] Medindo FLOPS ($GROUP_FLOPS)..."
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_FLOPS -m $EXECUTABLE > "$RESULT_DIR/flops/N_${n}.log" 2>&1

    # C. Cache Miss (L2 ou L1)
    echo "  [3/4] Medindo Cache L2 ($GROUP_L2)..."
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_L2 -m $EXECUTABLE > "$RESULT_DIR/l2/N_${n}.log" 2>&1

    # D. Banda de Memória (L3/MEM)
    echo "  [4/4] Medindo Memória ($GROUP_MEM)..."
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_MEM -m $EXECUTABLE > "$RESULT_DIR/mem/N_${n}.log" 2>&1

done

# --- 3. Consolidação dos Resultados ---

echo ""
echo "--- 3. Resumo Preliminar (Tempo de Iteração) ---"
echo "N     | Tempo Iteração (s)"
echo "------|-------------------"

for n in "${SIZES[@]}"; do
    # Extração simples do tempo usando grep/tail
    # Assume que o tempo é impresso próximo ao final do log
    
    if [ -f "$RESULT_DIR/tempo/N_${n}.log" ]; then
        TIME_VAL=$(grep -v "Nao calculado" "$RESULT_DIR/tempo/N_${n}.log" | tail -n 2 | head -n 1)
        printf "%-5d | %s\n" "$n" "$TIME_VAL"
    else
        printf "%-5d | Erro/Log ausente\n" "$n"
    fi
done

echo ""
echo "Benchmark finalizado."