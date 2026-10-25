# Plano de Melhorias

Roteiro para levar o protótipo ao nível de um artigo. As fases estão em ordem de **custo-benefício**: primeiro o que é barato e fortalece a defesa, depois o que pode mudar as conclusões.

Cada fase diz **o que fazer**, **onde mexer no código**, **como saber que terminou** e **qual handout estudar antes**.

---

## Visão geral

| Fase | Melhoria | Por quê | Esforço | Handout |
|---|---|---|---|---|
| 0 | Testes de robustez rápidos | Blindar a conclusão principal contra críticas óbvias | 1–2 dias | P1, P3 |
| 1 | Análise de tamanho finito em Fortaleza | É o que um físico estatístico mais valoriza | 3–5 dias | **T4** |
| 2 | Grafo dirigido (mão única) | **72% das ruas são de mão única**; a análise atual superestima a resiliência | 4–6 dias | **T6** |
| 3 | Reroteamento (equilíbrio de Wardrop) | Pode inverter o resultado do "plano fixo" | 5–8 dias | **T5** |
| 4 | Dados reais de velocidade | Transforma "mecanismo" em "Fortaleza" | depende do acesso | P1 |
| 5 | Quão bom é o guloso? | Responde "e o ótimo?" | 2–3 dias | T3 |
| 6 | Artigo (Frente 2) | Entregável final | 1–2 semanas | todos |

**Descoberta que mudou as prioridades:** com o sentido das ruas respeitado (componentes fortemente conexos), o limiar em Fortaleza cai de q_c ≈ 0,48 para **q_c ≈ 0,25** (handout T6, Seção 5). Isso coloca a Fase 2 acima do que estava no sumário anterior.

---

## Fase 0: testes de robustez rápidos

**Objetivo:** mostrar que "ligações críticas vencem em tempo real" não é um artefato das escolhas do modelo.

| Teste | Mudança no código | Resultado esperado | O que seria preocupante |
|---|---|---|---|
| Intervenção parcial: `r → r + 0,2` em vez de `r → 1` | novo argumento `--delta` em `com_intervencao` e `guloso_ligacoes_criticas` | ganhos menores, mesma ordem das estratégias | ligações críticas perdendo para o betweenness em tempo real |
| Orçamento menor: 0,2% | já existe: `--orcamento 0.002` | vantagem relativa das críticas ainda maior | — |
| Varredura de ψ (fração crônica): 0,1 até 0,9 | laço externo sobre `FRAC_CRONICA` | ΔR_plano / ΔR_tempo_real cresce de ≈0 até 1 | não crescer |
| Embaralhar velocidades | `rng.permutation(r)` antes de `limiar_critico` | q_c sobe de 0,475 para ≈0,50 | não mudar (correlações irrelevantes) |
| Teste estatístico pareado | trocar o erro padrão por diferenças dia a dia | significância maior | — |
| 160 dias em vez de 40 | `--dias 160` | barras de erro pela metade | — |

**Pronto quando:** existir uma tabela "teste → resultado" no P3, e uma figura de ΔR_plano/ΔR_tempo_real em função de ψ.

---

## Fase 1: análise de tamanho finito em Fortaleza

**Objetivo:** medir ν e d_f na malha real e dizer se ela está na classe da percolação 2D.

1. **Recortes concêntricos.** Criar `codigo/tamanho_finito.py`, que gera recortes quadrados de lado 1, 2, 4 e 8 km centrados na Praça do Ferreira, a partir do extrato do Ceará (`recortar_osm.sh` com o BBOX como parâmetro).
2. **Medidas por recorte**, com r uniforme (percolação pura) e depois com o modelo de trânsito:
   - largura σ[q_c] → 1/ν;
   - ⟨S₂^max⟩ → d_f;
   - número de pontes com massa de corte > 10% → 1/ν (estimativa independente);
   - ⟨q_c(L)⟩ → extrapolação para L → ∞.
3. **Ajustes log-log com bootstrap** (T4, Seção 4) e **colapso de dados** de S₁/N.
4. **Validar antes na rede quadrada** (já feito no T4: 1/ν = 0,728 e d_f = 1,889).

**Pronto quando:** existir uma figura com 4 painéis (largura, S₂^max, red bonds e colapso) e uma tabela de expoentes com barras de erro.

**Risco:** só cabem 3 a 4 tamanhos, e os bairros mudam de morfologia com L (T4, Problema 7). Diga isso explicitamente. Comparar "Fortaleza" com "grade de mesmo N" é um controle útil.

---

## Fase 2: grafo dirigido (mão única)

**Objetivo:** refazer toda a análise respeitando os sentidos.

1. **Malha dirigida:** nova função `malha_osm_dirigida`, igual a `malha_osm` mas guardando um `nx.DiGraph` (um arco por sentido). O código do T6 já faz isso; é só portar.
2. **Velocidade por arco:** o modelo de velocidades passa a sortear r por arco. A betweenness já é dirigida se o grafo for.
3. **Percolação de SCCs:** como o union-find não serve (T6, Problema 6), usar uma grade de 400 valores de q com `connected_components(connection='strong')`. Custa décimos de segundo por dia.
4. **Ligação crítica dirigida:** localizar a queda de S₂ na grade e testar os arcos do intervalo um a um (T6, solução 6c).
5. **Rodar as 5 estratégias de novo** e comparar com os resultados não-dirigidos.

**Pronto quando:** existir uma tabela "não-dirigido vs. dirigido" para q_c, R e ΔR de cada estratégia, e um mapa mostrando quantas das ligações críticas dirigidas são de mão única.

**Pergunta científica nova:** quais ruas de mão única, se voltassem a ser de mão dupla, mais aumentariam q_c? É uma "intervenção" diferente de destravar, e muito barata para a prefeitura.

---

## Fase 3: reroteamento (equilíbrio de Wardrop)

**Objetivo:** substituir a carga "tudo ou nada" por um equilíbrio de usuário.

1. **MSA** (T5, Seção 5) em nova função `carga_equilibrio`:
   - começa com a carga atual (betweenness);
   - a cada iteração calcula t_e = t⁰_e / r_e(x_e), refaz os caminhos mínimos (com amostragem de origens, como hoje) e mistura com peso 1/(n+1);
   - para quando a carga muda menos de 1%.
2. **Equilíbrio por dia:** o ruído do dia entra na capacidade, e o equilíbrio é recalculado.
3. **Intervenção com reroteamento:** depois de destravar as vias, recalcular o equilíbrio. Isso testa o paradoxo de Braess na prática (T5, Seção 4).

**Pronto quando:** existir uma comparação "sem vs. com reroteamento" da concentração de carga (por exemplo, o coeficiente de Gini da carga) e da tabela de estratégias.

**Risco:** o custo computacional sobe de 10 a 30 vezes. Comece com menos dias e menos origens.

**O que pode mudar:** a vitória do betweenness no plano fixo vem do corredor Padre Valdevino estar sempre travado. Com reroteamento, a carga se espalha, e esse resultado pode desaparecer.

---

## Fase 4: dados reais de velocidade

**Objetivo:** trocar r_e simulado por r_e medido.

Fontes a investigar (confirmar disponibilidade e licença):
- **GPS dos ônibus** (Etufor): velocidades em vias com linhas de ônibus. Perguntar ao Prof. José Soares, que como Cientista-Chefe de Dados do Ceará deve conhecer o acesso.
- **Portal de dados abertos da Prefeitura de Fortaleza** e **AMC/CTAFOR** (laços indutivos nos semáforos).
- **Waze for Cities** (programa de parceria com prefeituras).

Passos: (1) casar as medições com as vias do OSM (*map matching*); (2) calcular r_e = v_e / percentil 95 por intervalo de 15 minutos, como em Li *et al.*; (3) rodar a percolação ao longo do dia e ver q_c(t), com a transição no rush.

**Pronto quando:** existir uma figura de q_c(t) ao longo de um dia real.

---

## Fase 5: quão bom é o guloso?

1. **Redes pequenas** (grade 8×8): busca exaustiva dos melhores k = 2 e 3 arcos e comparação com o guloso.
2. **Fortaleza:** recozimento simulado a partir da solução gulosa, com ΔR como objetivo.
3. **Cota superior:** o ganho de R de destravar *todas* as vias com r < q_c dá um teto.

**Pronto quando:** existir a frase "o guloso alcança X% do melhor encontrado".

---

## Fase 6: artigo

Estrutura em revtex (PRE / PNAS):
1. Introdução: Li *et al.*, resiliência urbana e controle.
2. Dados e modelo: a malha OSM, o modelo de velocidades (P1) e a versão dirigida (Fase 2).
3. Métodos: varredura exata, o teorema do pico de S₂ e as red bonds (P2).
4. Resultados: transição, tamanho finito (Fase 1), estratégias, tempo real vs. plano fixo.
5. Discussão: gargalos que evoluem, implicações para a CTAFOR, limitações.
6. Material suplementar: testes de robustez (Fase 0) e validação dos algoritmos.

**Figura principal sugerida:** mapa de Fortaleza com as ligações críticas mais o painel ΔR × orçamento.

---

## Ordem recomendada

```
Fase 0 ──► Fase 1 ──► Fase 2 ──► Fase 3 ──► Fase 6
   (paralelo: Fase 4 = pedir dados já; Fase 5 = quando sobrar tempo)
```

Peça os dados (Fase 4) **logo na primeira reunião**. Eles demoram, e o resto do plano não depende deles.
