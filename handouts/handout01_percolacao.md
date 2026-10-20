# Handout 1 — Percolação em Redes: limiar crítico, expoentes e o segundo maior cluster

> **Projeto:** Transições de fase de percolação dinâmica e gargalos críticos na malha viária de Fortaleza
> **Pré-requisitos:** probabilidade básica (independência, valor esperado), séries geométricas, derivadas, expansão de Taylor.
> **Como estudar:** leia a Parte A inteira; faça as Deduções Guiadas (Parte C) **no papel, sem olhar o gabarito**; depois a Lista de Problemas (Parte D). O gabarito completo está na Parte E.
> **Ligação com o código:** cada seção indica a função de `codigo/otimizacao_antes_depois.py` onde a ideia aparece.

---

## Sumário

- A. Intuição física
- B. Formalismo
  - B.1 Percolação de ligações
  - B.2 Parâmetro de ordem, tamanho médio dos clusters e $S_2$
  - B.3 Expoentes críticos e leis de escala
  - B.4 Tamanho finito: por que $S_2$ tem um pico
  - B.5 Red bonds
  - B.6 Do tráfego à percolação: o protocolo de Li et al. (2015)
- C. Deduções guiadas (D1–D5)
- D. Lista de problemas (P1–P5)
- E. Gabarito detalhado
- F. Referências

---

## A. Intuição física

Pegue um tabuleiro quadriculado enorme. Cada aresta do tabuleiro (cada "rua" entre dois cruzamentos) está **aberta** com probabilidade $p$ e **fechada** com probabilidade $1-p$, independentemente das outras. Pergunta: existe um caminho de ruas abertas que atravessa o tabuleiro inteiro?

- Com $p$ pequeno, só existem ilhas pequenas de ruas abertas. Ninguém atravessa a cidade.
- Com $p$ perto de 1, quase tudo está aberto: existe um **cluster gigante** que ocupa uma fração finita da cidade.
- Entre os dois existe um valor **exato** $p_c$ em que o comportamento muda **abruptamente** (no limite de tabuleiro infinito). Para a rede quadrada, $p_c = 1/2$ exatamente.

Isso é uma **transição de fase**, igualzinha à da água virando gelo ou de um ímã perdendo a magnetização no ponto de Curie, só que **geométrica**: não há energia nem temperatura, só conectividade. O papel da magnetização é feito pela fração de nós no cluster gigante.

Três imagens para guardar:

1. **Abaixo de $p_c$:** muitas ilhas, a maior delas é pequena (cresce só como $\ln N$).
2. **Exatamente em $p_c$:** os clusters são **fractais**, cheios de buracos em todas as escalas, e o maior deles é um "esqueleto" frágil, cheio de **pontes únicas** (as *red bonds*): cortar uma só delas parte o cluster em dois.
3. **Acima de $p_c$:** um continente (o componente gigante) e ilhotas. O continente engole as ilhas conforme $p$ cresce.

**A ponte com o trânsito.** Numa hora de rush, cada via $e$ tem uma velocidade relativa $r_e = v_e / v_e^{\text{livre}} \in [0,1]$. Escolha um nível de exigência $q$ e declare a via **funcional** se $r_e \ge q$. Subir $q$ é exatamente o mesmo que **diminuir** $p$: fechamos vias da mais lenta para a mais rápida. Existe então um $q_c$ em que a cidade "se parte" em pedaços desconectados de vias fluidas. Quanto **maior** $q_c$, mais exigente pode ser o motorista e a cidade ainda continua conectada: **$q_c$ mede a resiliência da cidade**.

E o ponto que faz o projeto funcionar: em $q_c$ o cluster crítico é pendurado por poucas red bonds. **Consertar poucas vias, as certas, muda $q_c$ muito mais do que consertar muitas vias, as erradas.**

---

## B. Formalismo

### B.1 Percolação de ligações

Seja $G=(V,E)$ um grafo com $N=|V|$ nós. Na **percolação de ligações** cada aresta é mantida (aberta) independentemente com probabilidade $p$. Os **clusters** são os componentes conexos do subgrafo de arestas abertas.

- $n_s(p)$: número de clusters de tamanho $s$ (nós) **por nó** da rede.
- $s\,n_s$: probabilidade de que um nó escolhido ao acaso pertença a um cluster (finito) de tamanho $s$.
- $S_1$, $S_2$: tamanhos (em nós) do maior e do segundo maior cluster.

(Na **percolação de sítios** os nós é que são abertos com probabilidade $p$. Os limiares são diferentes, mas os expoentes críticos são os mesmos: mesma **classe de universalidade**.)

### B.2 Parâmetro de ordem, tamanho médio dos clusters e $S_2$

**Parâmetro de ordem.** A probabilidade de um nó pertencer ao cluster infinito:

$$P_\infty(p) = \lim_{N\to\infty} \frac{S_1}{N}, \qquad P_\infty = 0 \ \text{ para } p<p_c, \qquad P_\infty>0 \ \text{ para } p>p_c .$$

Como todo nó está em algum cluster: $\sum_s s\,n_s + P_\infty = 1$ (soma só sobre clusters finitos).

**Tamanho médio do cluster** (a "susceptibilidade"). Escolha um nó ao acaso que **não** esteja no gigante. O tamanho médio do cluster ao qual ele pertence é

$$\chi(p) = \frac{\sum_s s^2 n_s}{\sum_s s\, n_s}.$$

Por que $s^2$? Porque escolher um **nó** ao acaso favorece clusters grandes: a chance de cair num cluster de tamanho $s$ é $\propto s\,n_s$, e aí o tamanho é $s$. Daí $\sum s \cdot (s n_s)$. Em $p_c$, $\chi \to \infty$: é o análogo exato da susceptibilidade magnética que diverge no ponto de Curie.

**O segundo maior cluster $S_2$.** Numa simulação finita, $\chi$ é ruidoso e exige excluir o gigante "à mão". O truque usado por Li et al. (2015), que vem de Margolina, Herrmann & Stauffer (1982), é olhar $S_2$:

- abaixo de $p_c$ todos os clusters são pequenos, então $S_2$ é pequeno;
- acima de $p_c$ o gigante engole os rivais, então $S_2$ é pequeno de novo;
- **em $p_c$** os dois maiores clusters são ambos fractais de tamanho $\sim L^{d_f}$ e estão prestes a se fundir: **$S_2$ tem um pico**.

Definimos, numa amostra, $p_c$ (ou $q_c$) como **a posição do pico de $S_2$**. É isso que a função `limiar_critico` calcula.

### B.3 Expoentes críticos e leis de escala

Perto de $p_c$, tudo vira lei de potência. Com $\Delta = p - p_c$:

| Grandeza | Comportamento | 2D | Campo médio ($d\ge 6$, Erdős–Rényi) |
|---|---|---|---|
| Parâmetro de ordem | $P_\infty \sim \Delta^{\beta}$, $\Delta>0$ | $\beta = 5/36$ | $1$ |
| Susceptibilidade | $\chi \sim \lvert\Delta\rvert^{-\gamma}$ | $\gamma = 43/18$ | $1$ |
| Comprimento de correlação | $\xi \sim \lvert\Delta\rvert^{-\nu}$ | $\nu = 4/3$ | $1/2$ |
| Distribuição de tamanhos em $p_c$ | $n_s \sim s^{-\tau}$ | $\tau = 187/91$ | $5/2$ |
| Tamanho de corte | $s_\xi \sim \lvert\Delta\rvert^{-1/\sigma}$ | $\sigma = 36/91$ | $1/2$ |
| Dimensão fractal do cluster | $S_1(p_c) \sim L^{d_f}$ | $d_f = 91/48$ | $4$ |
| Dimensão das red bonds | $L_{\text{red}} \sim L^{1/\nu}$ | $1/\nu = 3/4$ | $2$ |

Os expoentes não são independentes. As **relações de escala** (todas verificáveis com a tabela):

$$\beta = \frac{\tau-2}{\sigma}, \qquad \gamma = \frac{3-\tau}{\sigma}, \qquad d\,\nu = \gamma + 2\beta \ \ (\text{hiperescala}), \qquad d_f = d - \frac{\beta}{\nu}.$$

**Universalidade.** Esses números não dependem dos detalhes (rede quadrada, triangular, ligações, sítios, malha viária com quarteirões tortos): dependem só da dimensão $d$ e de as correlações serem de curto alcance. Uma malha viária plana com desordem de curto alcance deve estar na classe da **percolação 2D**. Correlações de **longo** alcance (decaindo como $r^{-a}$ com $a < 2/\nu = 3/2$) mudam a classe (Weinrib, 1984): no nosso modelo o ruído é correlacionado só até alguns quarteirões, então a classe 2D vale.

### B.4 Tamanho finito: por que $S_2$ tem um pico (e onde)

Numa cidade de tamanho linear $L$ (com $N \sim L^2$ cruzamentos), a transição não é abrupta: é suavizada numa janela onde $\xi \gtrsim L$. A **hipótese de escala de tamanho finito** diz que

$$\frac{S_1}{N} = L^{-\beta/\nu}\, \mathcal{F}\!\left[(p-p_c)\,L^{1/\nu}\right], \qquad S_2^{\max} \sim L^{d_f}, \qquad p_c(L) - p_c \sim L^{-1/\nu}.$$

Consequências práticas para nós:

- $S_2^{\max}/N \sim L^{d_f - 2} = L^{-5/48}$ decai **muito** devagar. Para $L=45$: $45^{-5/48} \approx 0{,}67$. Por isso o pico de $S_2$ é bem visível (dezenas de % de $N$ em cidades reais).
- O $q_c$ medido numa amostra flutua de dia para dia com largura $\sim L^{-1/\nu} = L^{-3/4}$. Por isso o script faz média sobre muitos dias.
- **Cuidado:** média dos picos $\ne$ pico da média. $\langle q_c \rangle$ (média dos picos de cada dia) e a posição do pico de $\langle S_2(q) \rangle$ não precisam coincidir; na Fig. 1 do script elas diferem visivelmente. Saiba explicar por quê.

### B.5 Red bonds

No cluster crítico, uma **red bond** (Coniglio, 1982) é uma ligação **singularmente conectada**: todo caminho entre dois pontos distantes do cluster (por exemplo, entre os dois lados da cidade) passa por ela. Toda red bond é uma **ponte** (*bridge*) do cluster: removê-la desconecta o cluster.

**Cuidado: nem toda ponte é red bond.** O cluster crítico está cheio de galhos pendurados (*dangling ends*), ligados ao resto por uma única ligação. Essa ligação é uma ponte, mas cortá-la só solta um galho pequeno; ela não está no **backbone** (o caminho de transporte entre pontos distantes). Para separar as duas coisas usamos a **massa de corte** de cada ponte: o tamanho do pedaço **menor** que se solta. Red bonds de verdade separam dois pedaços macroscópicos (massa de corte grande); pontes de galhos têm massa de corte pequena. **Teste numérico** (rede quadrada, $r$ uniforme, média de 20 amostras no pico de $S_2$):

| $L$ | massa do cluster crítico | total de pontes | pontes com massa de corte $> 10\%$ do cluster |
|---|---|---|---|
| 20 | 255 | 135 | 20,8 |
| 40 | 1017 | 492 | 31,9 |
| razão | $4{,}0 \approx 2^{2{,}0}$ | $3{,}65 \approx 2^{1{,}87}$ | $1{,}53 \approx 2^{0{,}62}$ |

O total de pontes cresce aproximadamente como a massa do cluster (expoentes 1,87 e 2,0; o esperado assintoticamente é $d_f = 91/48 \approx 1{,}9$): é dominado pelos galhos. As pontes "macroscópicas" crescem como $L^{0{,}62}$, perto do $L^{3/4}$ de Coniglio (com $L$ tão pequeno, correções de escala são esperadas). É um ótimo exercício de análise de tamanho finito para repetir com $L = 20, 40, 80, 160$.

O cluster crítico tem a estrutura "**links–nós–blobs**": blobs (pedaços com laços redundantes) ligados em série por red bonds. O resultado notável é que o número de red bonds escala como $L_{\text{red}} \sim L^{1/\nu} = L^{3/4}$ em 2D, enquanto a massa do cluster escala como $L^{91/48}$. A **fração** de vias que são red bonds vai a zero:

$$\frac{L_{\text{red}}}{E} \sim \frac{L^{3/4}}{L^{2}} \to 0.$$

Esta é a justificativa física da hipótese do projeto: **poucas vias (bem menos que 1%) controlam a conectividade global no ponto crítico.** Você vai demonstrar o $L^{1/\nu}$ no Problema 5.

No código: `red_bonds` acha **todas** as pontes do cluster crítico e as ordena pela massa de corte, de modo que as red bonds verdadeiras ficam no topo da lista e as pontes de galhos ficam no fim.

### B.6 Do tráfego à percolação: o protocolo de Li et al. (2015)

1. Para cada via $e$ e cada instante, $r_e = v_e / v_e^{\max}$ ($v_e^{\max}$ = velocidade de fluxo livre, p.ex. percentil 95 do dia).
2. Para cada $q \in [0,1]$: mantenha as vias com $r_e \ge q$.
3. Meça $S_1(q)$ e $S_2(q)$. O pico de $S_2$ dá $q_c$.
4. As vias que, ao serem "removidas" logo acima de $q_c$, partem o gigante são os **gargalos** (*bottlenecks*). Melhorar a velocidade só delas sobe $q_c$.

**Dicionário:** fração de vias abertas $p(q) = \Pr(r_e \ge q) = 1 - F(q)$, onde $F$ é a função de distribuição acumulada de $r$. Se os $r_e$ fossem **independentes**, a malha seria percolação comum com $p = 1 - F(q)$ e

$$\boxed{\,q_c = F^{-1}(1 - p_c)\,}$$

(Dedução D4). Na rede quadrada, $q_c$ seria **a mediana** de $r$. Na cidade real os $r_e$ são correlacionados (congestionamento se espalha), e o desvio de $q_c$ em relação a $F^{-1}(1-p_c)$ mede o efeito das correlações.

**Algoritmo de Newman–Ziff.** Em vez de recalcular componentes para cada $q$ (custo $O(E)$ por valor de $q$), ordenamos as vias da mais rápida para a mais lenta e as **adicionamos uma a uma** com *union-find* (`limiar_critico`). Cada adição custa praticamente $O(1)$, e a varredura **inteira** em $q$ custa $O(E \log E)$ (o $\log$ vem da ordenação).

---

## C. Deduções guiadas

> Preencha as lacunas $\boxed{\ ?\ }$ no papel. Só depois confira na Parte E.

### D1. A cadeia 1D: tudo exato

Uma cadeia infinita de nós, ligações abertas com probabilidade $p$.

**(a)** Um cluster de exatamente $s$ nós é uma sequência de $s-1$ ligações abertas cercada por duas ligações fechadas. Por nó da rede:
$$n_s = \boxed{\ ?\ } .$$

**(b)** Verifique que $\sum_{s\ge1} s\,n_s = 1$ para $p<1$. Use $\sum_{s\ge 1} s x^{s-1} = 1/(1-x)^2$. O que isso diz sobre $P_\infty$ para $p<1$? Conclua que $p_c = \boxed{\ ?\ }$.

**(c)** Mostre que
$$\chi = \frac{\sum s^2 n_s}{\sum s n_s} = \frac{1+p}{1-p}.$$
Dica: $\sum_{s\ge1} s^2 x^{s-1} = \dfrac{1+x}{(1-x)^3}$. Leia o expoente: $\chi \sim (p_c - p)^{-\gamma}$ com $\gamma = \boxed{\ ?\ }$.

**(d)** A probabilidade de dois nós a distância $r$ estarem no mesmo cluster é $g(r) = p^r = e^{-r/\xi}$. Ache $\xi(p)$ e expanda perto de $p=1$ para obter $\nu = \boxed{\ ?\ }$.

### D2. A árvore de Bethe (rede sem laços, cada nó com $z$ vizinhos)

**(a) Limiar.** Chegando num nó por uma ligação, há $z-1$ ligações novas para continuar, cada uma aberta com probabilidade $p$. O número médio de "filhos" é $\boxed{\ ?\ }$. Um processo de ramificação sobrevive se o número médio de filhos for $>1$. Logo
$$p_c = \frac{1}{z-1}.$$

**(b) Parâmetro de ordem para $z=3$.** Seja $t$ a probabilidade de que um **ramo** (a ligação que sai de um nó numa direção específica, mais tudo o que está além dela) **não** leve ao infinito. O ramo é finito se a ligação está fechada **ou** se está aberta e os $z-1 = 2$ sub-ramos do nó seguinte são ambos finitos:
$$t = \boxed{\ ?\ } .$$
Mostre que $t$ satisfaz $p t^2 - t + (1-p) = 0$, com raízes $t=1$ e $t = (1-p)/p$. Para $p>1/2$ a raiz física é a menor. Um nó está no cluster infinito se **algum** de seus 3 ramos for infinito:
$$P_\infty = 1 - t^3 = 1 - \left(\frac{1-p}{p}\right)^3 .$$
Expanda com $p = \tfrac12 + \varepsilon$ e mostre que $P_\infty \approx 12\,\varepsilon$, logo $\beta = \boxed{\ ?\ }$.

**(c) Susceptibilidade.** Para $p<p_c$, seja $B$ o número médio de nós alcançados por um ramo. Justifique $B = p\,[\,1 + (z-1) B\,]$ e conclua
$$\chi = 1 + zB = \frac{1+p}{1-(z-1)p},$$
que diverge em $p_c$ com $\gamma = \boxed{\ ?\ }$.

### D3. Hiperescala e dimensão fractal em 2D

Usando $\beta = 5/36$, $\gamma = 43/18$, $\nu = 4/3$:

**(a)** Verifique $d\nu = \gamma + 2\beta$ para $d=2$.
**(b)** Calcule $d_f = d - \beta/\nu$ e confirme $91/48$.
**(c)** Justifique $d_f = d - \beta/\nu$ com o argumento: dentro de uma caixa de lado $\xi$, a massa do gigante é $P_\infty\, \xi^d$, e essa massa é a de um fractal de tamanho $\xi$, ou seja, $\xi^{d_f}$.

### D4. Do tráfego para a percolação: $q_c = F^{-1}(1-p_c)$

Suponha os $r_e$ independentes e identicamente distribuídos com distribuição acumulada $F(q) = \Pr(r_e \le q)$, contínua.

**(a)** A probabilidade de uma via estar aberta no nível $q$ é $p(q) = \boxed{\ ?\ }$.
**(b)** A cidade está no limiar quando $p(q_c) = p_c$. Isole $q_c$.
**(c)** Para a rede quadrada ($p_c = 1/2$), mostre que $q_c$ é a **mediana** de $r$.
**(d)** Se $r = h(x)$ com $h$ **estritamente decrescente** (mais carga $x$, menos velocidade), mostre que $\text{mediana}(r) = h(\text{mediana}(x))$. Por que isso **não** vale para a média?

### D5. Por que o pico de $S_2$ marca uma única ligação (o teorema do código)

Adicione as ligações uma a uma (da maior para a menor $r$), juntando clusters. Quando a ligação $e$ junta clusters de tamanhos $s_a \ge s_b$, chame $m_e = s_b = \min(s_a, s_b)$.

**(a)** Mostre que, no instante dessa fusão, $S_2 \ge m_e$. (Dica: logo antes da fusão existem **dois** clusters distintos com tamanho $\ge s_b$.)
**(b)** Mostre que $S_2$ só pode **diminuir** quando os dois maiores clusters se fundem. Para isso, examine os três casos possíveis de fusão: (i) dois clusters que não são nem o 1º nem o 2º maior; (ii) o 2º maior com um menor; (iii) o 1º com o 2º.
**(c)** Conclua que
$$\max_{q} S_2(q) = \max_{e}\, m_e,$$
e que a ligação $e^\ast$ que realiza o máximo é aquela cuja adição destrói o pico. Seu $r_{e^\ast}$ é o $q_c$.

(É exatamente esta a linha `if tam[b] > melhor: melhor, critica = tam[b], e` do código.)

---

## D. Lista de problemas

### P1. Dualidade e o $p_c = 1/2$ da rede quadrada (heurístico)

A rede quadrada $\mathcal{L}$ tem uma rede dual $\mathcal{L}^\ast$ (nós nos centros das faces). Cada ligação de $\mathcal{L}$ cruza exatamente uma ligação de $\mathcal{L}^\ast$. Declare uma ligação dual **aberta** se a ligação original que ela cruza está **fechada**.

**(a)** Se as ligações de $\mathcal{L}$ estão abertas com probabilidade $p$, com que probabilidade estão abertas as ligações de $\mathcal{L}^\ast$?
**(b)** Num quadrado $L\times L$, argumente que **ou** existe caminho aberto da esquerda para a direita em $\mathcal{L}$, **ou** existe caminho aberto de cima para baixo em $\mathcal{L}^\ast$, mas nunca os dois (desenhe!).
**(c)** Supondo que haja um único limiar, conclua que $p_c = 1/2$.
**(d)** A malha sintética do script retira 8% das ruas locais. Para ela, o limiar de ligações é maior ou menor que $1/2$? Explique fisicamente.

### P2. O $q_c$ de uma cidade com carga log-normal (cálculo explícito)

Modelo simplificado de uma única classe de via (sem correlações), igual ao do script:
$$r = \frac{1}{1+\alpha x^{4}}, \qquad x = \bar{x}\, e^{\sigma \eta - \sigma^2/2}, \qquad \eta \sim \mathcal{N}(0,1).$$

**(a)** Mostre que $\langle e^{\sigma\eta - \sigma^2/2} \rangle = 1$. Por que o termo $-\sigma^2/2$ é fisicamente necessário?
**(b)** Para uma rede com limiar $p_c$, mostre que
$$q_c = \left[\,1 + \alpha\, \bar{x}^{4}\, e^{4\sigma\,\Phi^{-1}(p_c) - 2\sigma^2}\right]^{-1},$$
onde $\Phi$ é a distribuição acumulada da normal padrão.
**(c)** Calcule $q_c$ para $\alpha = 0{,}15$, $\bar{x} = 1{,}5$, $\sigma = 0{,}5$, com $p_c = 0{,}5$ e com $p_c = 0{,}55$. ($\Phi^{-1}(0{,}55) \approx 0{,}1257$.)
**(d)** Se a demanda sobe 10% ($\bar{x} \to 1{,}65$), quanto cai $q_c$ (para $p_c = 0{,}5$)? Interprete: por que um aumento pequeno de demanda causa uma queda grande da resiliência?
**(e)** Aumentar $\sigma$ (mais flutuação), mantendo a carga média, aumenta ou diminui $q_c$ para $p_c = 1/2$? Isso é intuitivo?

### P3. Laços atrapalham? Bethe vs. rede quadrada

**(a)** Para a rede quadrada ($z = 4$), a aproximação de árvore de Bethe dá $p_c = 1/(z-1)$. Calcule e compare com o exato $1/2$.
**(b)** Explique por que, na rede real, uma ligação que fecha um **laço** não ajuda a aumentar o cluster, e por que isso faz o $p_c$ verdadeiro ser **maior** que a estimativa de árvore.
**(c)** Malhas viárias têm $\langle k \rangle \approx 2{,}5$–$3{,}5$ e muitos laços curtos (quarteirões). Fortaleza tem regiões de grade regular (Aldeota) e regiões irregulares com muitos becos sem saída. Qual tipo de região você espera ter $q_c$ maior, mantendo a mesma distribuição de velocidades? Justifique com as ideias de (b).

### P4. O paradoxo: melhorar uma via pode **diminuir** $q_c$

Considere três caminhos (cadeias) de nós:
- $A$: $a_1 - a_2 - a_3 - a_4 - a_5$ com $r(a_5a_4) = 0{,}99$, $r(a_4a_3) = 0{,}98$, $r(a_3a_2) = 0{,}97$, $r(a_2a_1) = 0{,}96$;
- $B$: $b_1 - \dots - b_5$ com $r(b_1b_2) = 0{,}99$, $r(b_2b_3)=0{,}98$, $r(b_3b_4)=0{,}97$, $r(b_4b_5)=0{,}96$;
- $C$: $c_1 - c_2 - c_3 - c_4$ com $r = 0{,}99,\ 0{,}98,\ 0{,}97$;

e duas pontes: $e = a_5 b_1$ com $r_e = 0{,}80$ e $f = b_5 c_1$ com $r_f = 0{,}50$.

**(a)** Usando D5, calcule $\max S_2$ e $q_c$.
**(b)** Agora "destrave" a via $e$: $r_e \to 1$. Recalcule $\max S_2$ e $q_c$.
**(c)** Mostre que, apesar disso, $S_1(q)$ **não diminuiu para nenhum** $q$. Em quais $q$ ele aumentou?
**(d)** Prove em geral: se aumentamos o $r$ de qualquer conjunto de vias, $S_1(q)$ não diminui para nenhum $q$. Conclua que a **robustez**
$$R = \frac{1}{N}\int_0^1 S_1(q)\, dq$$
(Schneider, Moreira, Andrade, Havlin & Herrmann, 2011) é **monótona**, enquanto $q_c$ não é. Por isso o script reporta as duas.

### P5. Por que há tão poucas red bonds: o argumento de Coniglio via fórmula de Russo

Seja $\Pi_L(p)$ a probabilidade de existir um caminho aberto atravessando um quadrado $L\times L$ ("evento $\mathcal{A}$"). Uma ligação $e$ é **pivotal** para $\mathcal{A}$ numa configuração se mudar o estado de $e$ (aberta ↔ fechada) muda se $\mathcal{A}$ ocorre ou não.

**Fórmula de Russo** (admita): para eventos crescentes,
$$\frac{d\Pi_L}{dp} = \langle \#\{\text{ligações pivotais}\} \rangle .$$

**(a)** Mostre que "$e$ é pivotal" **não depende** do estado de $e$ (depende só das outras ligações). Conclua que $\langle \#\{\text{pivotais abertas}\}\rangle = p\,\langle \#\{\text{pivotais}\}\rangle$.
**(b)** Argumente que, quando $\mathcal{A}$ ocorre, as pivotais **abertas** são exatamente as red bonds do cluster que atravessa.
**(c)** Por escala de tamanho finito, $\Pi_L(p) = \Phi\!\left[(p - p_c) L^{1/\nu}\right]$ para alguma função $\Phi$. Calcule $d\Pi_L/dp$ em $p_c$ e conclua que
$$\langle L_{\text{red}} \rangle \sim L^{1/\nu}.$$
**(d)** Para a malha sintética do script ($L = 45$, $E \approx 3700$), estime a ordem de grandeza do número de red bonds e a fração de vias que elas representam. Compare com o orçamento de 1% usado no script.

---

## E. Gabarito detalhado

### D1 (cadeia 1D)

**(a)** $n_s = (1-p)^2\, p^{\,s-1}$. São $s-1$ ligações internas abertas ($p^{s-1}$) e duas ligações das pontas fechadas ($(1-p)^2$). "Por nó" porque, ancorando o cluster pela sua ponta esquerda, cada nó da rede é ponta esquerda de no máximo um cluster.

**(b)** $\sum_s s\,n_s = (1-p)^2 \sum_s s\,p^{s-1} = (1-p)^2 \cdot \dfrac{1}{(1-p)^2} = 1$. Todo nó está num cluster finito, então $P_\infty = 0$ para qualquer $p<1$. Logo $p_c = 1$: em 1D basta **uma** ligação fechada para cortar a cadeia.

**(c)** $\sum_s s^2 n_s = (1-p)^2 \dfrac{1+p}{(1-p)^3} = \dfrac{1+p}{1-p}$; o denominador é 1. Então $\chi = \dfrac{1+p}{1-p} \approx \dfrac{2}{p_c - p}$, ou seja, $\gamma = 1$.

**(d)** $p^r = e^{r \ln p}$, logo $\xi = -1/\ln p$. Com $p = 1 - \delta$: $\ln p \approx -\delta$, então $\xi \approx 1/\delta = (p_c - p)^{-1}$ e $\nu = 1$.

### D2 (Bethe)

**(a)** $(z-1)p$. Sobrevive se $(z-1)p > 1$, ou seja, $p_c = 1/(z-1)$.

**(b)** $t = (1-p) + p\,t^2$: ligação fechada (prob. $1-p$) **ou** ligação aberta (prob. $p$) e os dois sub-ramos independentes finitos (prob. $t^2$). Rearranjando, $p t^2 - t + (1-p) = 0$. Fórmula de Bhaskara:
$$t = \frac{1 \pm \sqrt{1 - 4p(1-p)}}{2p} = \frac{1 \pm \lvert 2p - 1\rvert}{2p} \in \left\lbrace 1,\ \frac{1-p}{p} \right\rbrace .$$
Para $p > 1/2$, $(1-p)/p < 1$ é a raiz física (a probabilidade de extinção é a **menor** raiz, resultado clássico de processos de ramificação). Então $P_\infty = 1 - t^3$. Com $p = \tfrac12 + \varepsilon$:
$$\frac{1-p}{p} = \frac{\tfrac12 - \varepsilon}{\tfrac12 + \varepsilon} \approx 1 - 4\varepsilon, \qquad t^3 \approx 1 - 12\varepsilon, \qquad P_\infty \approx 12\,\varepsilon .$$
Linear em $\varepsilon$: $\beta = 1$ (campo médio).

**(c)** O ramo contribui se a ligação estiver aberta (prob. $p$): ganha 1 nó (o vizinho) mais os $z-1$ sub-ramos do vizinho, cada um com média $B$. Logo $B = p\,[1 + (z-1)B]$, então $B = \dfrac{p}{1-(z-1)p}$ e
$$\chi = 1 + zB = \frac{1 - (z-1)p + zp}{1-(z-1)p} = \frac{1+p}{1-(z-1)p}.$$
Diverge como $(p_c - p)^{-1}$: $\gamma = 1$.

### D3 (hiperescala)

**(a)** $\gamma + 2\beta = \tfrac{43}{18} + \tfrac{10}{36} = \tfrac{43}{18} + \tfrac{5}{18} = \tfrac{48}{18} = \tfrac{8}{3}$ e $d\nu = 2 \cdot \tfrac43 = \tfrac83$. ✓

**(b)** $\beta/\nu = \tfrac{5}{36}\cdot\tfrac{3}{4} = \tfrac{5}{48}$; $d_f = 2 - \tfrac{5}{48} = \tfrac{91}{48} \approx 1{,}896$. ✓

**(c)** $P_\infty\, \xi^d \sim \xi^{d_f}$ e $P_\infty \sim \Delta^\beta \sim \xi^{-\beta/\nu}$ (pois $\Delta \sim \xi^{-1/\nu}$). Então $\xi^{d - \beta/\nu} = \xi^{d_f}$.

### D4 ($q_c = F^{-1}(1-p_c)$)

**(a)** $p(q) = \Pr(r_e \ge q) = 1 - F(q)$.
**(b)** $1 - F(q_c) = p_c \Rightarrow q_c = F^{-1}(1 - p_c)$.
**(c)** $p_c = 1/2 \Rightarrow F(q_c) = 1/2$: é a definição de mediana.
**(d)** Se $h$ é decrescente, $r \le h(m)$ se e só se $x \ge m$. Com $m$ = mediana de $x$: $\Pr(r \le h(m)) = \Pr(x \ge m) = 1/2$, logo $h(m)$ é a mediana de $r$. Para a média isso falha porque $h$ não é linear: $\langle h(x) \rangle \ne h(\langle x \rangle)$ (desigualdade de Jensen). Mediana (e qualquer quantil) comuta com funções monótonas; a média não.

### D5 (o teorema do pico de $S_2$)

**(a)** Logo antes da fusão os clusters $a$ e $b$ são distintos e ambos têm tamanho $\ge s_b$. Então o segundo maior cluster tem pelo menos $s_b$ nós: $S_2 \ge s_b = m_e$.

**(b)** Sejam $C_1 \ge C_2 \ge C_3 \ge \dots$ os tamanhos antes da fusão.
- (i) Fundir dois clusters fora do top-2: os dois maiores continuam lá e um cluster novo, maior, aparece. $S_2$ não diminui.
- (ii) Fundir $C_2$ com um cluster menor: $C_2$ cresce. Se ultrapassar $C_1$, o novo $S_2$ é o antigo $C_1 \ge C_2$. $S_2$ não diminui.
- (iii) Fundir $C_1$ com $C_2$: o novo $S_2$ é $C_3 \le C_2$. **Só aqui $S_2$ pode diminuir.**

(Fusão de $C_1$ com um cluster que não é $C_2$: o top-2 só ganha massa, $S_2$ não muda.)

**(c)** Pelo item (a), $\max_q S_2 \ge \max_e m_e$. Ao final do processo (grafo conexo) $S_2 = 0$, então o pico tem de ser destruído por alguma fusão. Pelo item (b), essa fusão é do tipo (iii), e nela $m_e = C_2 = S_2^{\max}$. Logo $\max_q S_2 = \max_e m_e$, e a ligação que realiza esse máximo é exatamente a que destrói o pico. Como adicionamos ligações em ordem decrescente de $r$, essa fusão acontece em $q = r_{e^\ast}$, que é o $q_c$. $\blacksquare$

**Comentário para a lousa:** o algoritmo é $O(E\,\alpha(N))$ (praticamente linear) e calcula **exatamente** o pico de $S_2$, sem grade em $q$. O teste `codigo/teste_algoritmos.py` confere isso contra força bruta.

### P1 (dualidade)

**(a)** $1 - p$.
**(b)** Um caminho dual de cima a baixo "corta" todas as ligações originais que ele cruza (estão fechadas), separando esquerda de direita; e vice-versa. Topologicamente, num retângulo, exatamente um dos dois cruzamentos ocorre (é o mesmo argumento do jogo de tabuleiro Hex: alguém sempre ganha, e nunca os dois).
**(c)** Para o quadrado, os dois eventos têm a mesma probabilidade por simetria quando $p = 1 - p$, ou seja, em $p = 1/2$. Se o limiar fosse $p_c < 1/2$, para $p_c < p < 1/2$ teríamos travessia original quase certa e, como a dual tem probabilidade $1 - p > 1/2 > p_c$, travessia dual também quase certa: contradição com (b). Simetricamente para $p_c > 1/2$. Logo $p_c = 1/2$ (Sykes & Essam, 1964; prova rigorosa: Kesten, 1980).
**(d)** Maior que $1/2$. Ruas já removidas não podem conectar nada; cada rua restante precisa estar aberta com probabilidade maior para compensar. (Na aproximação mais simples, a probabilidade efetiva de uma ligação da rede quadrada estar aberta é $0{,}92\,p$ para as ruas locais, então seria preciso $0{,}92\,p \gtrsim 1/2$.)

### P2 (cidade log-normal)

**(a)** $\langle e^{\sigma \eta} \rangle = e^{\sigma^2/2}$ (função geradora de momentos da normal: complete o quadrado na integral gaussiana). Então $\langle e^{\sigma\eta - \sigma^2/2} \rangle = 1$. Sem o termo, aumentar o ruído aumentaria a **carga média**, misturando dois efeitos diferentes (mais carros vs. mais flutuação). Com ele, $\sigma$ controla só a flutuação.

**(b)** $r$ é decrescente em $\eta$. $r \le q \iff \eta \ge \eta^\ast(q)$. Então $F(q) = 1 - \Phi(\eta^\ast)$. A condição $F(q_c) = 1 - p_c$ dá $\Phi(\eta^\ast) = p_c$, ou seja, $\eta^\ast = \Phi^{-1}(p_c)$. Substituindo $x = \bar{x}\, e^{\sigma\eta^\ast - \sigma^2/2}$ em $r$ e usando $x^4 = \bar{x}^4 e^{4\sigma\eta^\ast - 2\sigma^2}$, sai a fórmula.

**(c)** Para $p_c = 0{,}5$: $x^\ast = 1{,}5\, e^{-0{,}125} = 1{,}3237$; $\alpha x^{\ast 4} = 0{,}4606$; $q_c = 1/1{,}4606 = 0{,}685$.
Para $p_c = 0{,}55$: $x^\ast = 1{,}5\, e^{0{,}5 \cdot 0{,}1257 - 0{,}125} = 1{,}4096$; $\alpha x^{\ast 4} = 0{,}5922$; $q_c = 0{,}628$.
Uma rede pior conectada (maior $p_c$) tem resiliência menor, como esperado.

**(d)** $\bar{x} = 1{,}65$: $x^\ast = 1{,}4561$, $\alpha x^{\ast 4} = 0{,}674$, $q_c = 0{,}597$. Queda de $0{,}088$, ou seja, $\approx 13\%$ de resiliência perdida para 10% de demanda a mais. O motivo é o expoente $\beta_{\text{BPR}} = 4$: $d\ln(\alpha x^4)/d\ln \bar{x} = 4$, então a **velocidade** é muito mais sensível que a **carga** perto da saturação. Em termos do $q_c$: $d q_c / d\ln\bar{x} = -4\, q_c (1 - q_c) \approx -0{,}86$.

**(e)** Para $p_c = 1/2$, $\eta^\ast = 0$ e $x^{\ast 4} \propto e^{-2\sigma^2}$: aumentar $\sigma$ **diminui** $x^\ast$ e **aumenta** $q_c$. Contraintuitivo? Não: com carga média fixa, um log-normal mais largo tem **mediana menor** (a média é puxada pela cauda de poucas vias muito carregadas). Como $q_c$ depende da mediana, a "via típica" fica mais leve. Moral: $q_c$ é controlado pelos **quantis** da distribuição de velocidades, não pela média. Isso só vale sem correlações; com correlações espaciais o resultado pode mudar (teste no script!).

### P3 (laços)

**(a)** $1/3 \approx 0{,}333 < 0{,}5$.
**(b)** Numa árvore, toda ligação aberta conecta um nó **novo**. Numa rede com laços, uma ligação pode ligar dois nós que já estavam no mesmo cluster: ela é "desperdiçada" para a conectividade (redundante). Precisamos então de mais ligações abertas para atingir a mesma conectividade: $p_c$ sobe. No algoritmo de Newman–Ziff, são as ligações em que `a == b` (o `continue`).
**(c)** A grade. Cuidado para não confundir com (b): lá comparamos redes com o **mesmo** $z$ (árvore vs. rede com laços). Aqui o que muda é o grau. Uma região com muitos becos sem saída tem $\langle k \rangle$ baixo e estrutura de árvore: cada beco está pendurado por uma única ligação (uma ponte), e qualquer via lenta no caminho desconecta tudo o que está além. A grade tem $\langle k \rangle \approx 4$ e caminhos alternativos em toda parte, então seu $p_c$ é menor e, pela relação $q_c = F^{-1}(1-p_c)$ de D4, seu $q_c$ é **maior**. (Ótima pergunta para levar ao Prof. José Soares: *como $q_c$ varia entre bairros com morfologias diferentes?*)

### P4 (paradoxo)

**(a)** As ligações entram em ordem decrescente de $r$. As cadeias crescem um nó por vez ($m_e = 1$). Em $q = 0{,}80$, $e$ funde $A$ (5) com $B$ (5): $m_e = 5$. Em $q = 0{,}50$, $f$ funde $A\cup B$ (10) com $C$ (4): $m_f = 4$. Então $\max S_2 = 5$ e $q_c = 0{,}80$.

**(b)** Com $r_e = 1$, $e$ entra primeiro, juntando $a_5$ e $b_1$ (dois nós isolados, $m = 1$). Depois $A\cup B$ cresce um nó por vez pelas duas pontas ($m = 1$). Em $q = 0{,}50$, $f$ funde 10 com 4: $m = 4$. Agora $\max S_2 = 4$ e $q_c = 0{,}50$. **$q_c$ caiu de $0{,}80$ para $0{,}50$!**

**(c)** Antes: para $q \in (0{,}80, 0{,}96]$, $S_1 \le 5$. Depois: nesse intervalo, $S_1$ chega a $10$. Para $q \le 0{,}80$ as duas situações são idênticas. $S_1$ só aumentou (para $q$ entre $0{,}80$ e $0{,}99$). A transição ficou mais **suave** (menos abrupta), e o "pico" passou a ser o evento seguinte. A rede **melhorou**, mas o estimador $q_c$ = pico de $S_2$ "pulou" para outro evento.

**(d)** Aumentar $r_e$ só pode **acrescentar** arestas ao conjunto $\lbrace e : r_e \ge q \rbrace$, para todo $q$. Acrescentar arestas a um grafo nunca diminui o maior componente: $S_1(q)$ é não-decrescente, ponto a ponto. Integrando, $R$ é não-decrescente. Já $q_c$ é a posição de um máximo (argmax), que pode pular de forma descontínua (item b). Por isso o script reporta $\Delta q_c$ (o observável de Li et al.) **e** $\Delta R$ (monótono e robusto). Nas simulações, $\Delta q_c$ de um único dia às vezes sai negativo; $\Delta R$ nunca (o teste `codigo/teste_algoritmos.py` verifica isso).

### P5 (red bonds e Russo)

**(a)** Pivotalidade de $e$ é definida comparando as duas configurações "$e$ aberta" e "$e$ fechada" com todo o resto igual: é uma função só das outras ligações. Como o estado de $e$ é independente das outras, $\Pr(e \text{ pivotal e aberta}) = p \Pr(e \text{ pivotal})$. Somando sobre $e$: $\langle \#\text{pivotais abertas} \rangle = p\, \langle \#\text{pivotais} \rangle$.

**(b)** Se $\mathcal{A}$ ocorre e $e$ é pivotal e aberta, fechar $e$ destrói **toda** travessia: $e$ está em todo caminho que atravessa. É a definição de red bond (ligação singularmente conectada) do cluster que atravessa. Reciprocamente, uma red bond aberta é pivotal.

**(c)** $\dfrac{d\Pi_L}{dp} = \Phi'\!\left[(p-p_c)L^{1/\nu}\right] L^{1/\nu}$; em $p = p_c$, $\dfrac{d\Pi_L}{dp} = \Phi'(0)\,L^{1/\nu}$. Pelos itens (a), (b) e Russo:
$$\langle L_{\text{red}} \rangle \approx p_c\, \frac{d\Pi_L}{dp}\Big|_{p_c} \sim L^{1/\nu} = L^{3/4} \quad (2D).$$
(Coniglio, 1982.) Em campo médio $\nu = 1/2$ dá $L^2$, mas aí "$L$" deve ser lido como a distância química.

**(d)** $L^{3/4} = 45^{0{,}75} \approx 17$ red bonds (a menos de uma constante de ordem 1). Fração: $17/3700 \approx 0{,}5\%$. O orçamento de 1% (37 vias) é da mesma ordem que o número de red bonds. Isso é consistente com a Fig. 3 do script: o ganho do método iterativo é muito rápido nos primeiros $\sim 0{,}2$–$0{,}4\%$ do orçamento e depois desacelera.

---

## F. Referências

1. D. Stauffer & A. Aharony, *Introduction to Percolation Theory*, 2ª ed., Taylor & Francis (1994). **Caps. 1–4 são a base deste handout.**
2. D. Li, B. Fu, Y. Wang, G. Lu, Y. Berezin, H. E. Stanley, S. Havlin, *Percolation transition in dynamical traffic network with evolving critical bottlenecks*, PNAS **112**, 669 (2015).
3. C. M. Schneider, A. A. Moreira, J. S. Andrade Jr., S. Havlin, H. J. Herrmann, *Mitigation of malicious attacks on networks*, PNAS **108**, 3838 (2011). (Medida de robustez $R$.)
4. M. E. J. Newman & R. M. Ziff, *Efficient Monte Carlo algorithm and high-precision results for percolation*, Phys. Rev. Lett. **85**, 4104 (2000).
5. A. Coniglio, *Cluster structure near the percolation threshold*, J. Phys. A **15**, 3829 (1982). (Red bonds e $L^{1/\nu}$.)
6. A. Margolina, H. J. Herrmann, D. Stauffer, *Size of largest and second largest cluster in random percolation*, Phys. Lett. A **93**, 73 (1982).
7. M. F. Sykes & J. W. Essam, *Exact critical percolation probabilities for site and bond problems in two dimensions*, J. Math. Phys. **5**, 1117 (1964).
8. A. Weinrib, *Long-range correlated percolation*, Phys. Rev. B **29**, 387 (1984).
9. R. Cohen, K. Erez, D. ben-Avraham, S. Havlin, *Resilience of the Internet to random breakdowns*, Phys. Rev. Lett. **85**, 4626 (2000). (Ponte para o Handout 2: critério de Molloy–Reed.)
