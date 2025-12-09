import os
import re
import matplotlib.pyplot as plt

# ==============================================================================
# CONFIGURAÇÕES
# ==============================================================================
pastas_raiz = ['resultadosT1', 'resultadosT2']

# Configuração visual dos gráficos
config_graficos = {
    'tempo_op1': {
        'ylabel': 'Tempo Médio (ms)',
        'title': 'Tempo de Execução OP1 (PCG)',
        'log_y': False, # ALTERADO: Eixo Y linear
        'log_x': True   # Eixo X logarítmico (padrão)
    },
    'tempo_op2': {
        'ylabel': 'Tempo (ms)',
        'title': 'Tempo de Execução OP2 (Resíduo)',
        'log_y': False, # ALTERADO: Eixo Y linear
        'log_x': True
    },
    'mem': {
        'regexs': [r'Memory bandwidth \[MBytes/s\]', r'L3 bandwidth \[MBytes/s\]'],
        'ylabel': 'MBytes/s',
        'title': 'Banda de Memória',
        'log_y': False,
        'log_x': True
    },
    'l2': {
        'regexs': [r'L2 miss ratio', r'L2 miss rate', r'data cache miss ratio'], 
        'ylabel': 'Miss Ratio',
        'title': 'Cache Miss L2',
        'log_y': False,
        'log_x': True
    },
    'flops': {
        # Agora buscamos duas métricas diferentes neste grupo
        'metricas_extras': {
            'Total': [r'DP MFLOP/s', r'MFLOP/s'],
            'AVX':   [r'AVX DP MFLOP/s', r'AVX MFLOP/s']
        },
        'ylabel': 'MFLOP/s',
        'title': 'Performance Aritmética (FLOPS)',
        'log_y': False,
        'log_x': True
    }
}

# ==============================================================================
# FUNÇÕES DE EXTRAÇÃO
# ==============================================================================
def extrair_likwid(conteudo, lista_regex):
    """Busca valor na tabela LIKWID."""
    for chave in lista_regex:
        padrao = r"\|\s+" + re.escape(chave).replace(r'\\', '') + r"\s+\|\s+([\d\.eE\+\-]+)"
        match = re.search(padrao, conteudo)
        if match:
            return float(match.group(1))
    return None

def extrair_tempos_c(conteudo):
    """Extrai os tempos do printf do C (últimas linhas)."""
    linhas = [l.strip() for l in conteudo.strip().split('\n') if l.strip()]
    if len(linhas) < 2: return None, None
    try:
        # Assume que os últimos prints são: ... -> tempoIter -> tResiduo
        t_op2_s = float(linhas[-1]) # Último
        t_op1_s = float(linhas[-2]) # Penúltimo
        return (t_op1_s * 1000.0), (t_op2_s * 1000.0) # Converte para ms
    except ValueError:
        return None, None

# ==============================================================================
# LEITURA E PROCESSAMENTO
# ==============================================================================
# Estrutura de dados para plotagem
# flops agora tem subcategorias: 'T1_Total', 'T1_AVX', etc.
dados = {
    'tempo_op1': {'T1': [], 'T2': []},
    'tempo_op2': {'T1': [], 'T2': []},
    'mem':       {'T1': [], 'T2': []},
    'l2':        {'T1': [], 'T2': []},
    'flops':     {'T1_Total': [], 'T1_AVX': [], 'T2_Total': [], 'T2_AVX': []}
}

for raiz in pastas_raiz:
    label_base = raiz.replace('resultados', '') # T1 ou T2
    
    # 1. TEMPOS
    path_tempo = os.path.join(raiz, 'tempo')
    if os.path.exists(path_tempo):
        arquivos = sorted([f for f in os.listdir(path_tempo) if f.startswith("N_")])
        for arq in arquivos:
            try:
                n = int(re.search(r'N_(\d+)', arq).group(1))
            except: continue
            
            with open(os.path.join(path_tempo, arq), 'r') as f: content = f.read()
            t1, t2 = extrair_tempos_c(content)
            if t1 is not None:
                dados['tempo_op1'][label_base].append((n, t1))
                dados['tempo_op2'][label_base].append((n, t2))

    # 2. LIKWID (Mem, L2)
    for key in ['mem', 'l2']:
        path = os.path.join(raiz, key)
        if not os.path.exists(path): continue
        for arq in sorted([f for f in os.listdir(path) if f.startswith("N_")]):
            try: n = int(re.search(r'N_(\d+)', arq).group(1))
            except: continue
            with open(os.path.join(path, arq), 'r') as f: content = f.read()
            
            val = extrair_likwid(content, config_graficos[key]['regexs'])
            if val is not None: dados[key][label_base].append((n, val))

    # 3. LIKWID (FLOPS com AVX)
    path_flops = os.path.join(raiz, 'flops')
    if os.path.exists(path_flops):
        for arq in sorted([f for f in os.listdir(path_flops) if f.startswith("N_")]):
            try: n = int(re.search(r'N_(\d+)', arq).group(1))
            except: continue
            with open(os.path.join(path_flops, arq), 'r') as f: content = f.read()
            
            # Extrai Total
            val_total = extrair_likwid(content, config_graficos['flops']['metricas_extras']['Total'])
            if val_total is not None:
                dados['flops'][f"{label_base}_Total"].append((n, val_total))
            
            # Extrai AVX
            val_avx = extrair_likwid(content, config_graficos['flops']['metricas_extras']['AVX'])
            if val_avx is not None:
                dados['flops'][f"{label_base}_AVX"].append((n, val_avx))

    # Ordenação
    for k in dados:
        if isinstance(dados[k], list): dados[k].sort(key=lambda x: x[0])
        else: 
            for subk in dados[k]: dados[k][subk].sort(key=lambda x: x[0])

# ==============================================================================
# PLOTAGEM
# ==============================================================================
fig = plt.figure(figsize=(14, 16)) # Aumentei altura para caber bem
gs = fig.add_gridspec(3, 2)

posicoes = {
    'tempo_op1': gs[0, 0], 'tempo_op2': gs[0, 1],
    'mem': gs[1, 0],       'l2': gs[1, 1],
    'flops': gs[2, :]
}

for key, pos in posicoes.items():
    ax = fig.add_subplot(pos)
    conf = config_graficos[key]
    
    # Lógica especial para FLOPS (tem 4 linhas)
    if key == 'flops':
        colors = {'T1': 'tab:blue', 'T2': 'tab:orange'}
        for subk, lista in dados[key].items():
            if not lista: continue
            versao = subk.split('_')[0] # T1 ou T2
            tipo = subk.split('_')[1]   # Total ou AVX
            
            estilo = '-' if tipo == 'Total' else '--'
            marcador = 'o' if tipo == 'Total' else 'x'
            
            X, Y = zip(*lista)
            ax.plot(X, Y, marker=marcador, linestyle=estilo, 
                    color=colors[versao], label=f"{versao} {tipo}")
    else:
        # Lógica padrão (T1 e T2)
        for versao in ['T1', 'T2']:
            lista = dados[key][versao]
            if lista:
                X, Y = zip(*lista)
                ax.plot(X, Y, marker='o', label=versao)

    ax.set_title(conf['title'], fontsize=12, fontweight='bold')
    ax.set_ylabel(conf['ylabel'])
    ax.set_xlabel('Tamanho N')
    ax.grid(True, which="both", ls="--", alpha=0.4)
    
    # Escalas
    if conf.get('log_x', True): ax.set_xscale('log')
    if conf.get('log_y', False): ax.set_yscale('log')
    
    ax.legend()

plt.tight_layout()
plt.savefig('graficos_finais_v3.png')
print("Gráficos gerados em 'graficos_finais_v3.png'")
plt.show()