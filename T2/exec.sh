#!/bin/bash

# --- Configurações Básicas ---
EXECUTABLE="./cgSolverOld"
RESULT_DIR="resultadosT1"
LIKWID_CMD="/home/soft/likwid/bin/likwid-perfctr"

# Parâmetros do Problema
K=7              # Número de diagonais (Definido, mas não usado no input abaixo)
MAXIT=25         # Limite de iterações
OMEGA=0.0        # Fator de relaxamento
EPSILON=1.0e-9   # Tolerância

# Definição dos Tamanhos (N)
SIZES=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# --- Configurações do LIKWID ---
# Grupos de performance a serem monitorados
GROUP_FLOPS="FLOPS_DP"
GROUP_L2="L2CACHE" 
GROUP_MEM="L3"

# --- Verificação de Ambiente ---

# Testa presença do LIKWID
if ! command -v $LIKWID_CMD &> /dev/null; then
    echo "ERRO: Ferramenta LIKWID não encontrada no caminho especificado."
    exit 1
fi

# Preparação de Diretórios
mkdir -p $RESULT_DIR/{tempo,flops,l2,mem}

# --- 1. Compilação do Projeto ---
echo "--- 1. Compilando (Otimização + LIKWID Marker) ---"

rm -f $EXECUTABLE

# O Makefile deve tratar as flags -DLIKWID_PERFMON
make

if [ $? -ne 0 ]; then
    echo "ERRO: Falha na compilação via Makefile."
    exit 1
fi
echo "Compilação concluída com sucesso."

# --- 2. Execução dos Benchmarks ---

echo ""
echo "--- 2. Iniciando Benchmarks ---"
echo "Params: MAXIT=$MAXIT, OMEGA=$OMEGA, EPSILON=$EPSILON"
echo "Saída: ./$RESULT_DIR/"

for n in "${SIZES[@]}"; do
    echo "========================================"
    echo "Processando N = $n"
    
    # Formatação da entrada (Stdin)
    # Nota: A variável K foi definida mas não está sendo passada aqui.
    # Confirme se o 'cgSolverOld' precisa ou não do K.
    INPUT_STR="$n $OMEGA $MAXIT $EPSILON"

    # A. Tempo de Execução (Wall Clock interno)
    echo "  [1/4] Medindo Tempo..."
    echo "$INPUT_STR" | $EXECUTABLE > "$RESULT_DIR/tempo/N_${n}.log"

    # B. Operações Aritméticas (FLOPS_DP)
    echo "  [2/4] Medindo FLOPS ($GROUP_FLOPS)..."
    echo "$INPUT_STR" | $LIKWID_CMD -C 0 -g $GROUP_FLOPS -m $EXECUTABLE > "$RESULT_DIR/flops/N_${n}.log" 2>&1

    # C. Cache Miss (L2)
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
    # Extração baseada na estrutura do log (Penúltima linha)
    if [ -f "$RESULT_DIR/tempo/N_${n}.log" ]; then
        TIME_VAL=$(grep -v "Nao calculado" "$RESULT_DIR/tempo/N_${n}.log" | tail -n 2 | head -n 1)
        printf "%-5d | %s\n" "$n" "$TIME_VAL"
    else
        printf "%-5d | Erro/Log ausente\n" "$n"
    fi
done

echo ""
echo "Script finalizado."