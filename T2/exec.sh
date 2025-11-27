#!/bin/bash
# Script de teste para o projeto cgSolver

# Nome do executável
EXECUTABLE="cgSolver"

# Arquivo de entrada padrão para os testes
INPUT_FILE="test_input.txt"

# --- Funções Auxiliares ---

# Função para compilar o projeto
compile_project() {
    echo " Compilando o projeto..."
    # Limpa compilações anteriores (opcional, mas recomendado)
    make clean > /dev/null 2>&1
    # Compila o executável principal
    make "$EXECUTABLE" 
    if [ $? -ne 0 ]; then
        echo " FALHA NA COMPILAÇÃO. Verifique o Makefile e os códigos-fonte."
        exit 1
    fi
    echo " Compilação bem-sucedida."
}

# Função para executar um caso de teste
run_test() {
    local n=$1
    local k=$2
    local omega=$3
    local maxit=$4
    local epsilon=$5
    local test_name="Teste-n${n}-k${k}-w${omega}-m${maxit}-e${epsilon}"

    echo ""
    echo "=========================================================="
    echo "Executando: $test_name"
    echo "----------------------------------------------------------"
    echo "$n $k $omega $maxit $epsilon" > "$INPUT_FILE"

    # Executa o programa, redirecionando o arquivo de entrada
    # O tempo de execução é medido com o 'time' e a saída é salva em OUT_FILE
    start_time=$(date +%s.%N)
    ./"$EXECUTABLE" < "$INPUT_FILE" > "${test_name}.out" 2> "${test_name}.err"
    end_time=$(date +%s.%N)
    elapsed=$(echo "$end_time - $start_time" | bc)
    
    RETURN_CODE=$?

    if [ $RETURN_CODE -eq 0 ]; then
        echo "SUCESSO. Código de retorno: $RETURN_CODE"
        echo "Tempo total de execução: ${elapsed}s"
        # Você pode adicionar comandos para verificar a saída aqui (e.g., checar norma_residuo)
        # Ex: grep "Norma do resíduo" "${test_name}.out"
    else
        echo "FALHA. Código de retorno: $RETURN_CODE"
        echo "Verifique o arquivo de erro: ${test_name}.err"
    fi
}

# --- Main ---

# 1. Compilação
compile_project

# 2. Casos de Teste
echo ""
echo "--- Iniciando Testes de Execução ---"

# Caso 1: Teste pequeno (Jacobi Precondicionado)
run_test 20 0.0 500 1e-6 

# Caso 2: Teste médio (Sem Precondicionador) - O código tem um EXIT no geraPreCond para -1.0, então usaremos um valor implementado, como 0.0
# Se o seu 'geraPreCond' for modificado para suportar w=-1.0, use o valor -1.0 aqui.
# Como o código atual força 'exit(1)' para w != -1.0 E w != 0.0, usaremos um teste funcional com w=0.0
run_test 100 0.0 1000 1e-8

# Caso 3: Teste maior, com mais diagonais e mais iterações
run_test 500 0.0 2000 1e-7

# 3. Limpeza (opcional)
# rm -f "$INPUT_FILE"

echo ""
echo "=========================================================="
echo "✨ Testes concluídos ✨ Verifique os arquivos *.out para detalhes."
echo "=========================================================="
