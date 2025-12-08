import os
import re
import matplotlib.pyplot as plt

# ==============================================================================
# CONFIGURAÇÕES
# ==============================================================================
pastas_raiz = ['resultadosT1', 'resultadosT2']

# Mapeamento: Pasta -> (Texto para buscar no LIKWID, Rótulo do Eixo Y, Título)
config_metricas = {
    'flops': {
        'regex_key': r'DP MFLOP/s', 
        'ylabel': 'MFLOP/s',
        'title': 'Performance de Cálculo (FLOPS)'
    },
    'mem': {
        # Ajuste se sua versão do Likwid usar "Memory bandwidth [MBytes/s]"
        'regex_key': r'Memory bandwidth \[MBytes/s\]', 
        'ylabel': 'MBytes/s',
        'title': 'Banda de Memória'
    },
    'l2': {
        # Ajuste se sua versão usar "L2 miss ratio" ou "L2 miss rate"
        'regex_key': r'L2 miss ratio', 
        'ylabel': 'Miss Ratio',
        'title': 'Taxa de Miss na Cache L2'
    },
    'tempo': {
        'regex_key': r'Runtime \(RDTSC\) \[s\]',
        'ylabel': 'Segundos',
        'title': 'Tempo de Execução'
    }
}

# ==============================================================================
# FUNÇÃO DE PARSER (AJUSTADA PARA SEU LOG)
# ==============================================================================
def extrair_valor_likwid(conteudo, chave_regex):
    """
    Busca uma linha no formato de tabela do LIKWID:
    |      DP MFLOP/s      |    4363.7043 |
    """
    # Regex explica:
    # \|             -> Procura uma barra vertical literal
    # \s+            -> Espaços em branco
    # chave_regex    -> O nome da métrica (ex: DP MFLOP/s)
    # \s+\|\s+       -> Espaços, barra vertical, espaços
    # ([\d\.eE\+\-]+)-> O VALOR (captura números, pontos e notação científica ex: 1.2e-05)
    padrao = r"\|\s+" + chave_regex + r"\s+\|\s+([\d\.eE\+\-]+)"
    
    match = re.search(padrao, conteudo)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return 0.0
    return None

# ==============================================================================
# PROCESSAMENTO
# ==============================================================================
dados_finais = {} # Estrutura: { 'flops': {'T1': [], 'T2': []}, ... }

for metrica, config in config_metricas.items():
    dados_finais[metrica] = {}
    
    for raiz in pastas_raiz:
        caminho_pasta = os.path.join(raiz, metrica)
        dados_finais[metrica][raiz] = []
        
        if not os.path.exists(caminho_pasta):
            print(f"Aviso: Pasta {caminho_pasta} não encontrada. Pulando.")
            continue

        arquivos = [f for f in os.listdir(caminho_pasta) if f.startswith("N_") and f.endswith(".log")]
        
        for arq in arquivos:
            # Extrai N do nome do arquivo (N_64.log -> 64)
            try:
                n_val = int(re.search(r"N_(\d+)", arq).group(1))
            except:
                continue

            caminho_arq = os.path.join(caminho_pasta, arq)
            with open(caminho_arq, 'r', encoding='utf-8', errors='ignore') as f:
                conteudo = f.read()
            
            # Extrai o valor usando a Regex do LIKWID
            valor = extrair_valor_likwid(conteudo, config['regex_key'])
            
            # FALLBACK: Se não achar tabela LIKWID (às vezes 'tempo' não tem likwid)
            # Tenta pegar o número "solto" na 3ª linha de escalares se for Tempo
            if valor is None and metrica == 'tempo':
                 # Tenta achar um float pequeno isolado no início do arquivo
                 # (Baseado no seu log: 9.0956688e-05)
                 match_loose = re.search(r"[\n\r](\d+\.\d+e?-?\d*)[\n\r]", conteudo)
                 if match_loose:
                     valor = float(match_loose.group(1))

            if valor is not None:
                dados_finais[metrica][raiz].append((n_val, valor))
        
        # Ordena por N
        dados_finais[metrica][raiz].sort(key=lambda x: x[0])

# ==============================================================================
# PLOTAGEM
# ==============================================================================
fig, axs = plt.subplots(2, 2, figsize=(15, 10))
axs = axs.flatten()

ordem_plot = ['flops', 'tempo', 'l2', 'mem']

for i, metrica in enumerate(ordem_plot):
    ax = axs[i]
    config = config_metricas.get(metrica, {})
    
    dados_metrica = dados_finais.get(metrica, {})
    
    tem_dados = False
    for raiz, valores in dados_metrica.items():
        if valores:
            tem_dados = True
            X = [v[0] for v in valores]
            Y = [v[1] for v in valores]
            # Limpa nome para legenda (resultadosT1 -> T1)
            label = raiz.replace('resultados', '')
            ax.plot(X, Y, marker='o', linestyle='-', label=label)
    
    if tem_dados:
        ax.set_title(config.get('title', metrica))
        ax.set_xlabel('Tamanho N (log)')
        ax.set_ylabel(config.get('ylabel', 'Valor'))
        ax.set_xscale('log') # Escala Log no X é essencial para N=64 a 20000
        ax.grid(True, which="both", ls="-", alpha=0.4)
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'Sem dados encontrados', ha='center')
        ax.set_title(metrica)

plt.tight_layout()
plt.savefig('graficos_desempenho.png')
print("Gráfico salvo: graficos_desempenho.png")
plt.show()