"""
teste_algoritmos.py
===================
Confere os algoritmos rápidos de otimizacao_antes_depois.py contra
FORÇA BRUTA (recalcular os componentes do zero para cada valor de q).

    python codigo/teste_algoritmos.py

Se tudo estiver certo, imprime três linhas "OK".
"""

import numpy as np
import networkx as nx
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from otimizacao_antes_depois import limiar_critico, red_bonds, malha_sintetica


def tamanhos_dos_clusters(N, arestas, abertas):
    """Tamanhos dos clusters (ordem decrescente) só com as arestas abertas."""
    A = coo_matrix((np.ones(abertas.sum()),
                    (arestas[abertas, 0], arestas[abertas, 1])), shape=(N, N))
    _, rotulo = connected_components(A, directed=False)
    return np.sort(np.bincount(rotulo))[::-1]


def teste_pico_S2_e_robustez(rng):
    """Pico de S2 = max min(s_a, s_b) (Handout 1, D5) e R = integral de S1."""
    for tentativa in range(30):
        if tentativa % 2:
            H = nx.convert_node_labels_to_integers(nx.grid_2d_graph(12, 12))
        else:
            H = nx.gnm_random_graph(150, 220, seed=tentativa)
        N = H.number_of_nodes()
        arestas = np.array(list(H.edges()))
        u, v = arestas[:, 0].tolist(), arestas[:, 1].tolist()
        r = rng.random(len(arestas))
        q_c, e_crit, pico, R = limiar_critico(N, u, v, r)

        # força bruta: S1 e S2 em cada degrau q = r_e
        pico_bruto, R_bruto, q_ant, S1_ant = 0, 0.0, 1.0, 1
        for q in np.sort(r)[::-1]:
            R_bruto += S1_ant * (q_ant - q)
            q_ant = q
            t = tamanhos_dos_clusters(N, arestas, r >= q)
            S1_ant = t[0]
            pico_bruto = max(pico_bruto, t[1] if len(t) > 1 else 0)
        R_bruto = (R_bruto + S1_ant * q_ant) / N

        # logo ANTES da ligação crítica entrar, S2 está no pico
        t = tamanhos_dos_clusters(N, arestas, r > r[e_crit])
        assert pico == pico_bruto == t[1]
        assert abs(R - R_bruto) < 1e-9
    print("OK: pico de S2, ligação crítica e R conferem com força bruta (30 grafos)")


def teste_R_monotono(rng):
    """Destravar uma via (r -> 1) nunca diminui R (Handout 1, P4d)."""
    H = malha_sintetica(L=20)
    N = H.number_of_nodes()
    arestas = np.array(list(H.edges()))
    u, v = arestas[:, 0].tolist(), arestas[:, 1].tolist()
    r = rng.random(len(arestas))
    R0 = limiar_critico(N, u, v, r)[3]
    for _ in range(200):
        r2 = r.copy()
        r2[rng.integers(len(arestas))] = 1.0
        assert limiar_critico(N, u, v, r2)[3] >= R0 - 1e-12
    print("OK: R nunca diminui ao destravar uma via (200 tentativas)")


def teste_red_bonds(rng):
    """Cada red bond é ponte do cluster crítico; massas em ordem decrescente."""
    H = malha_sintetica(L=20)
    N = H.number_of_nodes()
    arestas = np.array(list(H.edges()))
    u, v = arestas[:, 0].tolist(), arestas[:, 1].tolist()
    r = rng.random(len(arestas))
    q_c = limiar_critico(N, u, v, r)[0]
    rb = red_bonds(arestas, r, q_c)

    G = nx.Graph()
    G.add_edges_from(arestas[r >= q_c].tolist())
    C = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    massas = []
    for e in rb:
        a, b = arestas[e]
        C.remove_edge(a, b)
        assert not nx.has_path(C, a, b)          # é mesmo uma ponte
        massas.append(len(nx.node_connected_component(C, a)))
        massas[-1] = min(massas[-1], C.number_of_nodes() - massas[-1])
        C.add_edge(a, b)
    assert massas == sorted(massas, reverse=True)
    assert len(rb) == len(list(nx.bridges(C)))   # achou todas as pontes
    print(f"OK: {len(rb)} red bonds, todas pontes, maiores massas = {massas[:5]}")


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    teste_pico_S2_e_robustez(rng)
    teste_R_monotono(rng)
    teste_red_bonds(rng)
