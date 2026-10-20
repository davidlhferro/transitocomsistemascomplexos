"""
otimizacao_antes_depois.py
==========================

Percolação dinâmica de tráfego na malha viária de Fortaleza e teste de
intervenção "antes vs. depois".

Inspirado em: D. Li, B. Fu, Y. Wang, G. Lu, Y. Berezin, H. E. Stanley,
S. Havlin, "Percolation transition in dynamical traffic network with
evolving critical bottlenecks", PNAS 112, 669 (2015).

Ideia física em uma frase
-------------------------
Cada via e tem uma velocidade relativa r_e = v_e / v_e^livre em [0, 1].
Fixado um "nível de exigência" q, dizemos que a via está ABERTA se r_e >= q.
Aumentando q de 0 a 1 fechamos as vias da mais lenta para a mais rápida
(percolação de ligações). O componente gigante S1(q) colapsa em q_c, onde o
segundo maior cluster S2(q) tem seu pico. q_c alto = cidade que continua
conectada mesmo exigindo velocidades altas = cidade mais resiliente.

Pipeline
--------
 1. Malha: OSM (corredor Benfica-Centro-Aldeota) via osmnx, ou malha sintética
 2. Velocidades de rush: função BPR + carga ~ betweenness + ruído correlacionado
 3. Varredura de percolação: S1(q), S2(q) e q_c (algoritmo de Newman-Ziff)
 4. Red bonds (pontes do cluster crítico) vs. edge betweenness
 5. Intervenção: "destravar" k <= 1% das arestas (r_e -> 1) com 5 estratégias
    e medir o deslocamento Delta q_c (dentro e fora da amostra)

Uso
---
    python codigo/otimizacao_antes_depois.py              # tenta OSM, senão sintético
    python codigo/otimizacao_antes_depois.py --sintetico  # força a malha sintética
    python codigo/otimizacao_antes_depois.py --dias 60 --orcamento 0.005
"""

import argparse
import os
import time

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr


# =============================================================================
# 0. PARÂMETROS (todos com significado físico)
# =============================================================================

# Recorte Benfica-Centro-Aldeota: (oeste, sul, leste, norte) em graus
BBOX = (-38.548, -3.752, -38.490, -3.718)

# Velocidade de fluxo livre por classe de via (km/h), seguindo o CTB (Art. 61):
# trânsito rápido 80, arterial 60, coletora 40, local 30.
VEL_LIVRE = {"motorway": 80, "trunk": 60, "primary": 60, "secondary": 50,
             "tertiary": 40, "residential": 30, "unclassified": 30,
             "living_street": 20}

# Capacidade relativa (via local = 1). Avenidas têm mais faixas.
CAPACIDADE = {"motorway": 3.0, "trunk": 2.5, "primary": 2.0, "secondary": 1.6,
              "tertiary": 1.3}                    # demais classes -> 1.0
ARTERIAIS = {"motorway", "trunk", "primary", "secondary"}

# Função BPR (Bureau of Public Roads): t = t0 [1 + alpha (V/C)^beta]
# => r = v/v_livre = t0/t = 1 / (1 + alpha x^beta), com x = V/C.
ALPHA_BPR = 0.15
BETA_BPR = 4.0

DEMANDA = 2.0         # intensidade global do rush (multiplica a carga)
CARGA_FUNDO = 0.6     # carga local (garagens, entregas) que toda via recebe
SIGMA = 0.45          # amplitude do ruído log-normal na carga
KAPPA_ART = 0.5       # arteriais flutuam (1+kappa) vezes mais: semáforos,
                      # paradas de ônibus, conversões, incidentes
PASSOS_DIFUSAO = 12   # comprimento de correlação ~ sqrt(passos) quarteirões
FRAC_CORRELACIONADA = 0.7   # fração da variância do ruído que é espacial
FRAC_CRONICA = 0.5    # fração da variância que é crônica (igual todo dia)

Q_GRADE = np.linspace(0.0, 1.0, 201)   # grade de q para as curvas S1, S2


# =============================================================================
# 1. MALHA VIÁRIA
# =============================================================================

def malha_osm(arquivo_cache, arquivo_xml=None):
    """Baixa (ou lê do cache) a malha 'drive' do recorte e devolve um grafo
    simples não-direcionado com atributos t0 (s), classe, x, y (metros).
    Se arquivo_xml for dado (recorte feito com codigo/recortar_osm.sh),
    lê dele em vez de consultar o Overpass."""
    import osmnx as ox

    if os.path.exists(arquivo_cache):
        G = ox.load_graphml(arquivo_cache)
    elif arquivo_xml is not None:
        G = ox.graph_from_xml(arquivo_xml, bidirectional=False)
        G = ox.truncate.truncate_graph_bbox(G, BBOX)   # corta as pontas
        ox.save_graphml(G, arquivo_cache)
    else:
        # Servidor oficial do Overpass e dois espelhos públicos (se um cair)
        for url in ["https://overpass-api.de/api",
                    "https://overpass.kumi.systems/api",
                    "https://overpass.private.coffee/api"]:
            ox.settings.overpass_url = url
            try:
                G = ox.graph_from_bbox(BBOX, network_type="drive")
                break
            except Exception as erro:
                print(f"[aviso] {url} falhou ({type(erro).__name__})")
        else:
            raise ConnectionError("nenhum servidor Overpass respondeu")
        ox.save_graphml(G, arquivo_cache)

    G = ox.truncate.largest_component(G, strongly=False)
    G = ox.project_graph(G)            # coordenadas em metros (UTM)
    # Avenidas de pista dupla geram 2 nós por cruzamento: fundimos nós a < 12 m
    G = ox.consolidate_intersections(G, tolerance=12, rebuild_graph=True,
                                     dead_ends=False)
    G = ox.add_edge_speeds(G, hwy_speeds=VEL_LIVRE, fallback=30)
    G = ox.add_edge_travel_times(G)
    return multidigrafo_para_grafo(G)


def multidigrafo_para_grafo(G):
    """MultiDiGraph do osmnx -> nx.Graph simples.
    Projeção não-direcionada: uma ligação {u,v} existe se há via em algum
    sentido (conectividade fraca). Entre vias paralelas fica a mais rápida."""
    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        if u == v:
            continue                               # laços não conectam nada
        classe = d.get("highway", "residential")
        if isinstance(classe, list):
            classe = classe[0]
        classe = classe.replace("_link", "")       # alças herdam a classe
        t0 = float(d["travel_time"])
        nome = d.get("name", "")
        if isinstance(nome, list):
            nome = nome[0]
        if H.has_edge(u, v) and H[u][v]["t0"] <= t0:
            continue
        H.add_edge(u, v, t0=t0, classe=classe, nome=str(nome))
    for n, d in G.nodes(data=True):
        if n in H:
            H.nodes[n]["x"] = float(d["x"])
            H.nodes[n]["y"] = float(d["y"])
    H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
    return nx.convert_node_labels_to_integers(H)


def malha_sintetica(L=45, semente=7):
    """Reticulado quadrado diluído imitando um bairro planejado:
    avenidas (primary) a cada 9 quadras, coletoras (tertiary) entre elas,
    8% das ruas locais removidas (quadras fechadas, becos)."""
    rng = np.random.default_rng(semente)
    H = nx.Graph()
    for (i, j) in nx.grid_2d_graph(L, L).nodes:
        H.add_node((i, j), x=100.0 * j + rng.normal(0, 8),
                   y=100.0 * i + rng.normal(0, 8))
    for (a, b) in nx.grid_2d_graph(L, L).edges:
        linha = a[0] if a[0] == b[0] else None     # rua horizontal: linha fixa
        coluna = a[1] if a[1] == b[1] else None    # rua vertical: coluna fixa
        k = linha if linha is not None else coluna
        if k % 9 == 4:
            classe = "primary"
        elif k % 9 == 0:
            classe = "tertiary"
        else:
            classe = "residential"
            if rng.random() < 0.08:
                continue
        comprimento = 100.0 * (1 + 0.2 * rng.random())
        H.add_edge(a, b, classe=classe,
                   t0=comprimento / (VEL_LIVRE[classe] / 3.6))
    H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
    return nx.convert_node_labels_to_integers(H)


# =============================================================================
# 2. MODELO ESTOCÁSTICO DE VELOCIDADES NO RUSH
# =============================================================================

def carga_estrutural(H, arestas, amostras=500, semente=0):
    """Fluxo médio em cada via ~ edge betweenness (todos vão de todo lugar a
    todo lugar pelo caminho mais rápido em fluxo livre). Normalizado: <b> = 1.
    Usa 'amostras' origens sorteadas para não custar O(N*E)."""
    k = min(amostras, H.number_of_nodes())
    bc = nx.edge_betweenness_centrality(H, k=k, weight="t0", seed=semente)
    b = np.array([bc[(u, v)] if (u, v) in bc else bc[(v, u)]
                  for u, v in arestas])
    return b / b.mean()


def ruido_correlacionado(P, arestas, rng):
    """Campo gaussiano espacialmente correlacionado (média 0, variância 1).
    Ruído branco nos nós é 'difundido' por passeio aleatório preguiçoso:
    phi <- (phi + P phi)/2. Após m passos a correlação alcança ~sqrt(m) quadras
    (núcleo de calor discreto). Somamos uma parte branca nas arestas."""
    phi = rng.standard_normal(P.shape[0])
    for _ in range(PASSOS_DIFUSAO):
        phi = 0.5 * (phi + P @ phi)
    phi = (phi - phi.mean()) / phi.std()
    eta_espacial = 0.5 * (phi[arestas[:, 0]] + phi[arestas[:, 1]])
    eta_espacial = (eta_espacial - eta_espacial.mean()) / eta_espacial.std()
    eta_local = rng.standard_normal(len(arestas))
    c = FRAC_CORRELACIONADA
    return np.sqrt(c) * eta_espacial + np.sqrt(1 - c) * eta_local


def velocidade_relativa(carga, capacidade, arterial, eta):
    """r_e = 1 / (1 + alpha x_e^beta)   (BPR), com carga x_e log-normal:
    x_e = D (x0 + b_e/c_e) exp(sigma_e eta_e - sigma_e^2/2).
    O termo -sigma^2/2 garante <exp(...)> = 1: o ruído não muda a carga média."""
    sigma = SIGMA * (1 + KAPPA_ART * arterial)
    x = DEMANDA * (CARGA_FUNDO + carga / capacidade) \
        * np.exp(sigma * eta - sigma**2 / 2)
    return 1.0 / (1.0 + ALPHA_BPR * x**BETA_BPR)


# =============================================================================
# 3. PERCOLAÇÃO: q_c E LIGAÇÃO CRÍTICA (Newman-Ziff + union-find)
# =============================================================================

def limiar_critico(N, u, v, r):
    """Abre as ligações da MAIS RÁPIDA para a MAIS LENTA (q descendo de 1 a 0)
    juntando clusters com union-find.

    Fato-chave (handout T1, Exemplo 3): S2 só diminui quando os dois MAIORES
    clusters se fundem, e nesse instante S2 = tamanho do menor deles. Logo
        max_q S2(q) = max sobre fusões de min(s_a, s_b),
    e a ligação que realiza esse máximo é a LIGAÇÃO CRÍTICA. Seu r é q_c.

    De brinde integra a ROBUSTEZ  R = integral_0^1 S1(q) dq / N  (Schneider,
    Moreira, Andrade, Havlin & Herrmann, PNAS 2011). Entre dois valores
    consecutivos de r nada muda, então S1(q) é uma escada e a integral é soma.

    Retorna (q_c, índice da ligação crítica, pico de S2 em nós, R)."""
    pai = list(range(N))
    tam = [1] * N

    def raiz(i):
        while pai[i] != i:
            pai[i] = pai[pai[i]]          # compressão de caminho (halving)
            i = pai[i]
        return i

    melhor, critica = 0, -1
    S1, robustez, q_anterior = 1, 0.0, 1.0
    for e in np.argsort(-r).tolist():
        robustez += S1 * (q_anterior - r[e])   # degrau da escada S1(q)
        q_anterior = r[e]
        a, b = raiz(u[e]), raiz(v[e])
        if a == b:
            continue                       # ligação interna: fecha um ciclo
        if tam[a] < tam[b]:
            a, b = b, a                    # b é o cluster menor
        if tam[b] > melhor:
            melhor, critica = tam[b], e
        pai[b] = a                         # união por tamanho
        tam[a] += tam[b]
        S1 = max(S1, tam[a])
    robustez += S1 * q_anterior            # de q = 0 até o menor r
    return r[critica], critica, melhor, robustez / N


def curvas_S1_S2(N, arestas, r):
    """S1(q) e S2(q) (frações de N) na grade Q_GRADE, via componentes conexas
    (scipy, em C). Usado só para os gráficos."""
    S1, S2 = np.zeros(len(Q_GRADE)), np.zeros(len(Q_GRADE))
    for i, q in enumerate(Q_GRADE):
        aberta = r >= q
        A = coo_matrix((np.ones(aberta.sum()),
                        (arestas[aberta, 0], arestas[aberta, 1])), shape=(N, N))
        _, rotulo = connected_components(A, directed=False)
        tamanhos = np.sort(np.bincount(rotulo))[::-1]
        S1[i] = tamanhos[0] / N
        S2[i] = tamanhos[1] / N if len(tamanhos) > 1 else 0.0
    return S1, S2


# =============================================================================
# 4. RED BONDS: pontes do cluster crítico, ordenadas pela massa que desconectam
# =============================================================================

def red_bonds(arestas, r, q_c):
    """No ponto crítico (ligações com r >= q_c, incluindo a crítica), o cluster
    gigante é quase uma árvore de 'blobs'. Uma red bond é uma PONTE desse
    cluster: cortá-la separa o cluster em dois pedaços. Nem toda ponte é red
    bond no sentido de Coniglio: pontes que só seguram um galho pendurado
    (dangling end) soltam um pedaço pequeno. Por isso ordenamos pela massa.

    Teoria (handout T1, Problema 8): no ponto crítico há só ~L^(3/4) red bonds.
    Massa de corte = tamanho do menor pedaço. Calculamos todas de uma vez:
    contraímos cada blob (componente 2-aresta-conexo) num nó e ganhamos uma
    árvore cujas arestas são as pontes; a massa de uma ponte é o tamanho da
    subárvore abaixo dela (ou o complemento). Custo O(N + E).

    Retorna índices das red bonds, da maior para a menor massa de corte."""
    G = nx.Graph()
    for e in np.where(r >= q_c)[0]:
        G.add_edge(arestas[e, 0], arestas[e, 1], id=e)
    C = G.subgraph(max(nx.connected_components(G), key=len))
    pontes = list(nx.bridges(C))

    sem_pontes = nx.Graph(C)
    sem_pontes.remove_edges_from(pontes)
    blob = {}
    for i, comp in enumerate(nx.connected_components(sem_pontes)):
        for n in comp:
            blob[n] = i
    massa_sub = np.bincount(list(blob.values())).astype(float)

    T = nx.Graph()                                # árvore de blobs
    T.add_nodes_from(range(len(massa_sub)))
    for a, b in pontes:
        T.add_edge(blob[a], blob[b], id=C[a][b]["id"])
    pai_T = nx.dfs_predecessors(T, 0)
    pos_ordem = list(nx.dfs_postorder_nodes(T, 0))  # filhos antes dos pais
    for n in pos_ordem:
        if n != 0:
            massa_sub[pai_T[n]] += massa_sub[n]
    total = massa_sub[0]

    lista = []
    for n in pos_ordem:
        if n != 0:
            massa = min(massa_sub[n], total - massa_sub[n])
            lista.append((massa, T[n][pai_T[n]]["id"]))
    lista.sort(reverse=True)
    return np.array([e for _, e in lista], dtype=int)


# =============================================================================
# 5. INTERVENÇÃO: destravar k ligações (r_e -> 1) e medir Delta q_c
# =============================================================================

def com_intervencao(N, u, v, r, escolhidas):
    """(q_c, R) depois de destravar as vias escolhidas."""
    r_novo = r.copy()
    r_novo[escolhidas] = 1.0      # via destravada: fluida mesmo no rush
    q_c, _, _, robustez = limiar_critico(N, u, v, r_novo)
    return q_c, robustez


def guloso_ligacoes_criticas(N, u, v, r, k_max):
    """Estratégia cirúrgica (iterativa): acha a ligação crítica, destrava,
    recalcula q_c, acha a NOVA ligação crítica, ... k_max vezes.
    Devolve q_c e R após cada passo e a lista de ligações escolhidas."""
    r_novo = r.copy()
    qcs, Rs, escolhidas = [], [], []
    for passo in range(k_max + 1):
        q_c, e, _, robustez = limiar_critico(N, u, v, r_novo)
        qcs.append(q_c)
        Rs.append(robustez)
        if passo < k_max:
            r_novo[e] = 1.0
            escolhidas.append(e)
    return np.array(qcs), np.array(Rs), np.array(escolhidas, dtype=int)


def jaccard(X, Y):
    """Sobreposição entre dois conjuntos de vias: |X inter Y| / |X uniao Y|."""
    X, Y = set(np.asarray(X).tolist()), set(np.asarray(Y).tolist())
    return len(X & Y) / max(1, len(X | Y))


# =============================================================================
# 6. GRÁFICOS (estilo Phys. Rev. E)
# =============================================================================

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

# Paleta categórica validada (ordem fixa) + cinza neutro para o "antes"
COR = {"antes": "#52514e", "aleatoria": "#e87ba4", "lentas": "#eda100",
       "betweenness": "#2a78d6", "redbonds": "#1baf7a", "guloso": "#eb6834"}
MARCADOR = {"aleatoria": "v", "lentas": "s", "betweenness": "D",
            "redbonds": "^", "guloso": "o"}
NOME = {"aleatoria": "Random", "lentas": "Slowest (brute force)",
        "betweenness": "Top betweenness", "redbonds": "Red bonds (one shot)",
        "guloso": "Critical bonds (iterative)"}


def fig_transicao(S1, S2, qc_medio, titulo, arquivo):
    fig, ax = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
    ax[0].plot(Q_GRADE, S1.mean(0), color=COR["antes"], lw=1.5)
    ax[0].fill_between(Q_GRADE, S1.mean(0) - S1.std(0), S1.mean(0) + S1.std(0),
                       color=COR["antes"], alpha=0.2, lw=0)
    ax[0].set_ylabel(r"$S_1(q)/N$")
    ax[1].plot(Q_GRADE, S2.mean(0), color=COR["guloso"], lw=1.5)
    ax[1].set_ylabel(r"$\langle S_2(q)\rangle/N$")
    ax[1].set_xlabel(r"$q$ (required relative speed)")
    for a in ax:
        a.axvline(qc_medio, color="k", ls="--", lw=0.8)
    ax[0].text(qc_medio + 0.02, 0.85, rf"$\langle q_c\rangle = {qc_medio:.3f}$")
    ax[0].set_title(titulo, fontsize=8)
    fig.savefig(arquivo)
    plt.close(fig)


def fig_antes_depois(curvas, qcs, frac, arquivo):
    fig, ax = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
    for nome, (S1, S2) in curvas.items():
        rot = "Before" if nome == "antes" else NOME[nome]
        lw = 2.0 if nome in ("antes", "guloso") else 1.2
        ax[0].plot(Q_GRADE, S1.mean(0), color=COR[nome], lw=lw,
                   label=rf"{rot} ($q_c$={np.mean(qcs[nome]):.3f})")
        ax[1].plot(Q_GRADE, S2.mean(0), color=COR[nome], lw=lw)
    ax[0].set_ylabel(r"$S_1(q)/N$")
    ax[1].set_ylabel(r"$\langle S_2(q)\rangle/N$")
    ax[1].set_xlabel(r"$q$")
    ax[0].legend(fontsize=6, loc="lower left", handlelength=1.5)
    ax[0].set_title(rf"Before vs. after unlocking {100 * frac:.1f}% of roads",
                    fontsize=8)
    fig.savefig(arquivo)
    plt.close(fig)


def fig_orcamento(ks, E, dqc, dqc_plano, dR, dR_plano, arquivo):
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True, sharey="row")
    paineis = [(0, 0, dqc), (0, 1, dqc_plano), (1, 0, dR), (1, 1, dR_plano)]
    for lin, col, dados in paineis:
        for nome in ["aleatoria", "lentas", "betweenness", "redbonds", "guloso"]:
            m = dados[nome].mean(0)
            s = dados[nome].std(0) / np.sqrt(len(dados[nome]))   # erro padrão
            ax[lin, col].errorbar(100 * ks / E, m, yerr=s, color=COR[nome],
                                  marker=MARCADOR[nome], ms=3.5, lw=1.2,
                                  capsize=2, label=NOME[nome])
        ax[lin, col].axhline(0, color="k", lw=0.5)
    ax[0, 0].set_title("Real-time control (same day)", fontsize=8)
    ax[0, 1].set_title(r"Fixed plan (train days $\to$ test days)", fontsize=8)
    ax[0, 0].set_ylabel(r"$\Delta q_c$")
    ax[1, 0].set_ylabel(r"$\Delta R$,  with $R=\int_0^1 S_1(q)\,dq\,/N$")
    for col in (0, 1):
        ax[1, col].set_xlabel(r"Budget $k/E$ (%)")
    alcas, rotulos = ax[0, 0].get_legend_handles_labels()
    fig.legend(alcas, rotulos, loc="lower center", ncol=5, fontsize=7,
               bbox_to_anchor=(0.5, -0.04))
    fig.savefig(arquivo)
    plt.close(fig)


def fig_mapa(H, arestas, r, criticas, top_bc, q_c, arquivo):
    xy = np.array([[H.nodes[n]["x"], H.nodes[n]["y"]] for n in range(len(H))])
    seg = xy[arestas]                               # (E, 2, 2)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    # Cinza: claro = fluido, escuro = congestionado (r baixo)
    cinza = LinearSegmentedColormap.from_list("cinza", ["#111111", "#d4d4d4"])
    lc = LineCollection(seg, array=r, cmap=cinza, norm=plt.Normalize(0, 1),
                        linewidths=0.6)
    ax.add_collection(lc)
    ax.add_collection(LineCollection(seg[top_bc], colors=COR["betweenness"],
                                     linewidths=2.2, label="Top edge betweenness"))
    ax.add_collection(LineCollection(seg[criticas], colors=COR["guloso"],
                                     linewidths=2.8, label="Critical bonds"))
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(lc, ax=ax, shrink=0.7, label=r"$r_e = v_e/v_e^{\rm free}$")
    ax.legend(loc="upper right", fontsize=7, frameon=True)
    ax.set_title(rf"One rush-hour realization, $q_c$ = {q_c:.3f}", fontsize=8)
    fig.savefig(arquivo)
    plt.close(fig)


# =============================================================================
# 7. PROGRAMA PRINCIPAL
# =============================================================================

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--osm-xml", default=None,
                   help="recorte .osm local (ver codigo/recortar_osm.sh)")
    p.add_argument("--sintetico", action="store_true",
                   help="usar malha sintética (sem internet)")
    p.add_argument("--dias", type=int, default=40,
                   help="número de realizações do rush (dias)")
    p.add_argument("--orcamento", type=float, default=0.01,
                   help="fração máxima de arestas destravadas")
    p.add_argument("--semente", type=int, default=2015)
    args = p.parse_args()

    raiz_proj = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pasta_fig = os.path.join(raiz_proj, "figuras")
    pasta_dados = os.path.join(raiz_proj, "dados")
    os.makedirs(pasta_fig, exist_ok=True)
    os.makedirs(pasta_dados, exist_ok=True)
    rng = np.random.default_rng(args.semente)
    t_inicio = time.time()

    # ---- 1. Malha ----------------------------------------------------------
    H, rotulo = None, "sintetico"
    if not args.sintetico:
        try:
            H = malha_osm(os.path.join(pasta_dados, "benfica_centro_aldeota.graphml"),
                          args.osm_xml)
            rotulo = "fortaleza"
        except Exception as erro:
            print(f"[aviso] OSM indisponível ({type(erro).__name__}); "
                  "usando malha sintética.")
    if H is None:
        H = malha_sintetica()
    N, E = H.number_of_nodes(), H.number_of_edges()
    arestas = np.array(list(H.edges()), dtype=int)
    u, v = arestas[:, 0].tolist(), arestas[:, 1].tolist()
    print(f"Malha '{rotulo}': N = {N} cruzamentos, E = {E} vias, "
          f"<k> = {2 * E / N:.2f}")

    # ---- 2. Ingredientes fixos do modelo de velocidades --------------------
    classe = [H[a][b]["classe"] for a, b in arestas]
    capacidade = np.array([CAPACIDADE.get(c, 1.0) for c in classe])
    arterial = np.array([c in ARTERIAIS for c in classe], dtype=float)
    print("Calculando edge betweenness (carga estrutural)...")
    b = carga_estrutural(H, [tuple(e) for e in arestas], semente=args.semente)
    A = nx.to_scipy_sparse_array(H, nodelist=range(N), weight=None, format="csr")
    P = A.multiply(1.0 / A.sum(axis=1)[:, None]).tocsr()   # P = D^-1 A

    k_max = max(1, int(args.orcamento * E))
    ks = np.unique(np.linspace(0, k_max, 11).astype(int))
    top_bc = np.argsort(-b)[:k_max]            # estático: só depende do mapa

    # ---- 3. Gera os "dias" de rush ------------------------------------------
    # Desordem CONGELADA (quenched): igual todo dia (buraco, escola, terminal).
    # Desordem RECOZIDA (annealed): sorteada de novo a cada dia (chuva, batida).
    R = args.dias
    eta_cronico = ruido_correlacionado(P, arestas, rng)
    psi = FRAC_CRONICA
    r_dias = []
    for _ in range(R):
        eta = np.sqrt(psi) * eta_cronico \
            + np.sqrt(1 - psi) * ruido_correlacionado(P, arestas, rng)
        r_dias.append(velocidade_relativa(b, capacidade, arterial, eta))
    dias_treino = range(R // 2)          # 1a metade: "histórico" da CTAFOR
    dias_teste = range(R // 2, R)        # 2a metade: dias futuros

    estrategias = ["aleatoria", "lentas", "betweenness", "redbonds", "guloso"]
    dqc = {s: np.zeros((R, len(ks))) for s in estrategias}   # Delta q_c
    dR = {s: np.zeros((R, len(ks))) for s in estrategias}    # Delta R
    qc0, R0 = np.zeros(R), np.zeros(R)
    curvas = {s: (np.zeros((R, len(Q_GRADE))), np.zeros((R, len(Q_GRADE))))
              for s in ["antes"] + estrategias}
    qcs_kmax = {s: np.zeros(R) for s in ["antes"] + estrategias}
    gul_por_dia, rb_por_dia = [], []

    # ---- 4a. Controle em TEMPO REAL: escolhe e testa no mesmo dia ----------
    print(f"Simulando {R} dias, orçamento até k = {k_max} vias "
          f"({100 * k_max / E:.2f}% de E)...")
    for d in range(R):
        r = r_dias[d]
        q_c, _, pico, robustez = limiar_critico(N, u, v, r)
        qc0[d], R0[d] = q_c, robustez

        # listas ordenadas de cada estratégia (as k primeiras são destravadas)
        qcs_gul, Rs_gul, gul = guloso_ligacoes_criticas(N, u, v, r, k_max)
        rb = red_bonds(arestas, r, q_c)[:k_max]
        gul_por_dia.append(gul)
        rb_por_dia.append(rb)
        listas = {"aleatoria": rng.choice(E, k_max, replace=False),
                  "lentas": np.argsort(r)[:k_max],
                  "betweenness": top_bc,
                  "redbonds": rb,
                  "guloso": gul}

        for j, k in enumerate(ks):
            for s in estrategias:
                if s == "guloso":                          # já calculado
                    q_novo, R_novo = qcs_gul[k], Rs_gul[k]
                else:
                    q_novo, R_novo = com_intervencao(N, u, v, r, listas[s][:k])
                dqc[s][d, j] = q_novo - q_c
                dR[s][d, j] = R_novo - robustez

        # curvas S1, S2 antes vs depois, com o orçamento cheio
        curvas["antes"][0][d], curvas["antes"][1][d] = curvas_S1_S2(N, arestas, r)
        qcs_kmax["antes"][d] = q_c
        for s in estrategias:
            r_novo = r.copy()
            r_novo[listas[s]] = 1.0
            curvas[s][0][d], curvas[s][1][d] = curvas_S1_S2(N, arestas, r_novo)
            qcs_kmax[s][d] = q_c + dqc[s][d, -1]
        print(f"  dia {d + 1:2d}/{R}: q_c = {q_c:.3f}, R = {robustez:.3f}, "
              f"pico S2 = {pico:4d} nós | críticas: "
              f"Dq_c = {dqc['guloso'][d, -1]:+.3f}, DR = {dR['guloso'][d, -1]:+.4f}")

    # ---- 4b. PLANO FIXO: escolhe com dias de treino, testa em dias novos ---
    # Pontuação de Borda: via escolhida em 1o lugar ganha 1, em último ~0;
    # somamos sobre os dias de treino. Vias que são gargalo TODO dia sobem.
    def borda(listas_por_dia):
        pontos = np.zeros(E)
        for d in dias_treino:
            lista = listas_por_dia[d]
            pontos[lista] += 1.0 - np.arange(len(lista)) / k_max
        return np.argsort(-pontos, kind="stable")[:k_max]

    r_treino = np.mean([r_dias[d] for d in dias_treino], axis=0)
    plano = {"aleatoria": rng.choice(E, k_max, replace=False),
             "lentas": np.argsort(r_treino)[:k_max],   # mais lentas em média
             "betweenness": top_bc,
             "redbonds": borda(rb_por_dia),
             "guloso": borda(gul_por_dia)}
    dqc_plano = {s: np.zeros((len(dias_teste), len(ks))) for s in estrategias}
    dR_plano = {s: np.zeros((len(dias_teste), len(ks))) for s in estrategias}
    for i, d in enumerate(dias_teste):
        for j, k in enumerate(ks):
            for s in estrategias:
                q_novo, R_novo = com_intervencao(N, u, v, r_dias[d], plano[s][:k])
                dqc_plano[s][i, j] = q_novo - qc0[d]
                dR_plano[s][i, j] = R_novo - R0[d]

    # ---- 4c. Relatório -----------------------------------------------------
    rho, _ = spearmanr(b, np.mean(r_dias, axis=0))
    persist = [jaccard(gul_por_dia[d], gul_por_dia[d + 1]) for d in range(R - 1)]
    print("\n================ RESULTADOS ================")
    print(f"<q_c> antes = {qc0.mean():.4f} +- {qc0.std():.4f} (desvio entre dias)")
    print(f"<R>   antes = {R0.mean():.4f} +- {R0.std():.4f}")
    print(f"Spearman(betweenness, <r>)          = {rho:+.3f}")
    print(f"Jaccard(red bonds, top BC)          = "
          f"{np.mean([jaccard(x, top_bc) for x in rb_por_dia]):.3f}")
    print(f"Jaccard(ligações críticas, top BC)  = "
          f"{np.mean([jaccard(x, top_bc) for x in gul_por_dia]):.3f}")
    print(f"Jaccard(críticas dia d, dia d+1)    = {np.mean(persist):.3f}")
    for titulo, X, Y in [("Delta q_c", dqc, dqc_plano), ("Delta R", dR, dR_plano)]:
        print(f"\n{titulo} com k = {k_max} vias ({100 * k_max / E:.2f}% de E):")
        print(f"{'estratégia':30s} {'tempo real':>17s} {'plano fixo':>17s}")
        for s in estrategias:
            x, y = X[s][:, -1], Y[s][:, -1]
            print(f"{NOME[s]:30s} {x.mean():+.4f} +- {x.std() / np.sqrt(len(x)):.4f}"
                  f"  {y.mean():+.4f} +- {y.std() / np.sqrt(len(y)):.4f}")
    print("(+- = erro padrão da média sobre os dias)")

    # Onde ficam as vias do plano fixo? (nomes do OpenStreetMap)
    for s in ["guloso", "betweenness"]:
        print(f"\nPlano fixo, {NOME[s]} (10 primeiras):")
        for e in plano[s][:10]:
            a, c = arestas[e]
            nome = H[a][c].get("nome", "") or "(sem nome)"
            print(f"   {nome:40s} {H[a][c]['classe']:12s} <r> = {r_treino[e]:.2f}")

    # ---- 5. Figuras --------------------------------------------------------
    tag = rotulo
    nome_malha = "Fortaleza OSM" if rotulo == "fortaleza" else "synthetic grid"
    fig_transicao(curvas["antes"][0], curvas["antes"][1], qc0.mean(),
                  f"Rush-hour percolation ({nome_malha}, N={N}, {R} days)",
                  os.path.join(pasta_fig, f"{tag}_fig1_transicao.png"))
    fig_antes_depois(curvas, qcs_kmax, k_max / E,
                     os.path.join(pasta_fig, f"{tag}_fig2_antes_depois.png"))
    fig_orcamento(ks, E, dqc, dqc_plano, dR, dR_plano,
                  os.path.join(pasta_fig, f"{tag}_fig3_orcamento.png"))
    fig_mapa(H, arestas, r_dias[0], gul_por_dia[0], top_bc, qc0[0],
             os.path.join(pasta_fig, f"{tag}_fig4_mapa.png"))
    print(f"\nFiguras salvas em {pasta_fig}/  ({time.time() - t_inicio:.0f} s)")


if __name__ == "__main__":
    main()
