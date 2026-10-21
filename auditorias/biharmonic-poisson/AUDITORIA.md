# Auditoria — "On a biharmonic scalar equation with singular nonlinearity"

> **Arquivo temporário.** Esta auditoria fica no repositório só enquanto o
> artigo é revisado. Apagar a pasta `auditorias/biharmonic-poisson/` depois.

- **Manuscrito:** De Araujo, Ercole, Faria, Salirrosas — versão `DeAraujo_Faria_Manuel (new) (Copy).pdf`, 22 páginas (trechos em vermelho = revisão).
- **Problema:** $\Delta^2 u = \lambda u^{-\alpha} + u^p$ em $\Omega$, $u>0$, $u=\Delta u=0$ em $\partial\Omega$, $0<\alpha<1$.
- **Escopo:** todas as contas e cadeias de desigualdades das pp. 1–21, coerência entre resumo, enunciados e provas, referências.
- **Data:** 21/10/2026.

## Resumo

Não encontrei erro matemático que invalide algum resultado. Todas as contas
conferidas fecham (lista na seção 4). Os pontos abaixo são, em ordem de
importância:

1. **Uma melhora real:** a hipótese $N\le 3$ do Teorema 1.1 pode ser retirada (seção 1).
2. **Incoerência de enunciado:** o Corolário 1 diz convergência *fraca*; o resumo e a prova dão convergência *forte*. E a subsequência é desnecessária (seção 2.1).
3. Simplificações e ajustes de redação/notação (seções 2 e 3).

---

## 1. Sugestão principal: a priori estimativa em qualquer dimensão

O Teorema 1.1 (pp. 2, 10–11) usa a identidade de energia (3.8), Young e
$H^2\hookrightarrow L^\infty$, o que obriga $N\le 3$. Um argumento direto pela
representação de Green dá a cota $L^\infty$ **para todo $N$**, com constante
explícita.

Seja $u$ solução clássica positiva, $M=\|u\|_\infty$, $F=\lambda u^{-\alpha}+u^p$ e
$\psi(x)=\int_\Omega G(x,y)\phi_0(y)\,dy$. Como o próprio artigo mostra (p. 10)
que $z=-\Delta u\in H^1_0$ resolve $-\Delta z=F$ fracamente com $F\in H^{-1}$,
vale $u=\iint G(x,y)G(y,\xi)F(\xi)\,d\xi\,dy$.

- **Cota inferior.** $u^{-\alpha}\ge M^{-\alpha}$, logo $u\ge \lambda M^{-\alpha}\psi$. Por (2.6), $\phi_0\ge\varphi_1/\lambda_1$, então $\psi\ge\varphi_1/\lambda_1^2$ e
  $$u\ \ge\ \frac{\lambda}{\lambda_1^2M^{\alpha}}\,\varphi_1 .$$
- **Cota superior.** Daí $u^{-\alpha}\le \lambda^{-\alpha}\lambda_1^{2\alpha}M^{\alpha^2}\varphi_1^{-\alpha}$. Pelo Lema 2.1,
  $$M\ \le\ \underbrace{\lambda^{1-\alpha}\lambda_1^{2\alpha}C(\alpha)\|\phi_0\|_\infty}_{A}\,M^{\alpha^2}\;+\;\underbrace{\|\psi\|_\infty}_{B}\,M^{p}.$$
- Como $\alpha^2<1$ e $p<1$: $\;M\le\max\{(2A)^{1/(1-\alpha^2)},\,(2B)^{1/(1-p)}\}$.

Isso dá a cota uniforme e a barreira $u\ge\delta\varphi_1$ (Observação 2) sem
restrição de dimensão. A conta via energia (3.8) continua necessária para o
Corolário 1, então ela fica; só o Teorema 1.1 e o resumo mudam
("when $N\le3$" sai).

## 2. Pontos a corrigir

### 2.1 Corolário 1: fraca × forte, e subsequência (pp. 1, 3, 20–21)
- O enunciado (p. 3) diz "converges **weakly** in $H^2\cap H^1_0$"; o resumo (p. 1) diz "converge **strongly**" e a prova (p. 21) termina com convergência forte. Trocar o enunciado para "strongly".
- O limite é sempre o mesmo, $v=\varphi_1/\|\varphi_1\|_2$. Pelo argumento de subsequência de subsequência, **a família inteira** $v_p$ converge quando $p\to1^-$. Dá para tirar "up to a subsequence" do enunciado e "along a subsequence" do resumo.
- O enunciado ainda diz "positive **weak** solution"; a prova mostra que $v$ é suave e clássica. Pode dizer isso.

### 2.2 Lema 4.1: o limite existe (pp. 15–16)
O mesmo argumento aplicado a **qualquer** sequência $p_n\to1^-$ com
$K_{p_n}^{1-p_n}\to L'$ mostra $L'\ge\lambda_1^{2\alpha}$, $L'\ge\|\phi_0\|_\infty^2$
e, se $L'>\lambda_1^{2\alpha}$, $L'=\|\phi_0\|_\infty^2$. Logo todo ponto de
acumulação vale $\max\{\lambda_1^{2\alpha},\|\phi_0\|_\infty^2\}$ e
$$\lim_{p\to1^-}K_p^{1-p}=\max\{\lambda_1^{2\alpha},\|\phi_0\|_\infty^2\}.$$
(4.4) é verdadeira, só mais fraca que o necessário. No caso $\lambda_1\le1$ o
texto escreve "(4.4) gives $\lim K_p^{1-p}=\|\phi_0\|^2_\infty$", que exige
justamente essa passagem implícita; com a forma acima fica direto.

### 2.3 Termo redundante em $\min\{K_1,K_2\}$ (pp. 8, 11, 13)
Por (2.7), $\lambda_1\|\phi_0\|_\infty\ge1$, então $K_2=K_1/(\lambda_1\|\phi_0\|_\infty)\le K_1$. Assim:
- (3.5): $\min\{K_1,K_2\}=K_2$;
- (3.11) e $\lambda^*$ (p. 13): o termo $K_1$ no mínimo é redundante.

### 2.4 Referência cruzada errada (p. 11)
Observação 2: "the comparison argument following Lemma 2.2" → é a **Observação 1** (p. 5).

### 2.5 Colisões de notação
- p. 10: o $\epsilon$ da desigualdade de Young é o mesmo símbolo do parâmetro de regularização. Usar $\eta$.
- p. 14: a constante $C$ do valor mínimo $f(t_0)$ colide com $C(\alpha)$ e com o $C$ genérico. Usar $C_{\alpha,p}$.
- p. 5: $\delta_\lambda$ na prova vs. $\delta$ no enunciado do Lema 2.2.
- Regularizações diferentes: $(u^\alpha+\epsilon)^{-1}$ na Seção 3 e $(u+\epsilon)^{-\alpha}$ na Seção 4. O texto da p. 19 explica, mas unificar evitaria o comentário.

## 3. Sugestões de organização e redação (opcionais)

1. **Identidade de energia como lema.** A justificativa de (3.8) (p. 10) está dentro da subseção $0<p<1$, mas é usada na não existência ($p>1$, p. 13) e no Corolário 1 (p. 20). Vale virar um lema na Seção 2 ("toda solução clássica positiva satisfaz $z\in W^{1,\infty}\cap H^1_0$ e (3.8)"), válido para todo $p>0$.
2. **Lema 2.2.** O lema em si (para (2.8)) não é usado; o que se usa é o argumento, aplicado a (1.1) (Observações 1 e 2). Pode-se enunciar direto para (1.1).
3. **Parágrafo "Uniqueness" (p. 11).** Está dentro da prova do Teorema 1.1, mas não prova nada. Melhor como observação após o teorema, ou retirar.
4. **Teorema 1.2.** Comentar a lacuna entre $\lambda^*$ (existência) e $\lambda_{ne}$ (não existência), ex.: definir $\Lambda=\sup\{\lambda:\text{existe solução}\}$ e notar $\lambda^*\le\Lambda\le\lambda_{ne}$ como questão aberta.
5. **Proposição 1.** A restrição $\|\nabla u\|_\infty\le M_\Omega K_p$ em $\mathcal S$ não é necessária para Schauder (basta fechado, convexo e imagem relativamente compacta). É inofensiva; pode ficar ou sair.
6. p. 2: "torsion function" aparece na introdução mas $\phi_0$ não recebe esse nome na Seção 2. Acrescentar "(the torsion function)" em (2.5).
7. p. 4: "G is nonnegative" — é positiva; nada muda.
8. "Corollary 1" destoa da numeração "Theorem 1.1–1.3"; usar "Corollary 1.4".

## 4. O que foi conferido e fecha

| Trecho | Verificação |
|---|---|
| (2.6), (2.7) | $\varphi_1=\lambda_1G[\varphi_1]\le\lambda_1\phi_0$; $1=\|\varphi_1\|_\infty\le\lambda_1\|\phi_0\|_\infty$. Desigualdade estrita (p. 16) via $-\Delta(\lambda_1\phi_0-\varphi_1)=\lambda_1(1-\varphi_1)\ge0$, $\not\equiv0$ ✓ |
| (2.9) | $w=-\Delta u/\lambda_1$ dá o sistema com os fatores $1/\lambda_1$ corretos ✓ |
| Lema 2.2, Obs. 1 | $\Delta^2\omega\ge\lambda/(K^\alpha+1)-\delta\lambda_1^2\ge0$; dois princípios do máximo; constante $\lambda/(\lambda_1^2K^\alpha)$ para (1.1) ✓ |
| Lema 3.1 | Barreira $b=A(d-d^s)$, $s=2-\alpha$: $-\Delta b=A[s(s-1)d^{-\alpha}-(1-sd^{s-1})\Delta d]$ ✓; escolha de $\rho$ e $A$ independentes de $n$ ✓; estimativa de gradiente interior com $r=d/2$ dá $C(1+d^{1-\alpha}+d)$ ✓; Schauder global, Arzelà–Ascoli, bootstrap ✓ |
| (3.2)–(3.5) | (3.2) ⇔ $K_1^p\|\phi_0\|\le K_1/(2\|\phi_0\|)$; (3.4) ⇔ $\lambda^{1-\alpha}C(\alpha)/\tau^\alpha\le K_1/(2\|\phi_0\|)$; intervalo de $\tau$ não vazio pois $t^{1/\alpha}/(t^\alpha+1)\to\infty$ ✓ |
| Lema 3.2 | Quatro cotas da invariância de $\mathcal A$ ✓ |
| Existência $0<p<1$ | Compacidade $W^{2,q}\hookrightarrow C^{1,\gamma}$, continuidade, Schauder; $\lambda/(u_\epsilon^\alpha+\epsilon)\le\lambda^{1-\alpha}\tau^{-\alpha}\varphi_1^{-\alpha}$ ✓ |
| (3.8) | $z\in W^{1,\infty}\cap H^1_0$, Hardy $\int d^{-\alpha}|\psi|\le C\|\nabla\psi\|_2$, $F\in H^{-1}$ ✓ |
| Young (p. 10) | Constantes com expoentes $\tfrac{2}{1\mp\alpha}$, $\tfrac{2}{1\pm p}$ conferidas termo a termo ✓ |
| A priori ($N\le3$) | $\|u\|_2\le\lambda_1^{-1}\|\Delta u\|_2$, $H^2\hookrightarrow L^\infty$ ✓; Obs. 2 ✓ |
| (3.9)–(3.11), $\epsilon_0$ | $\lambda_1^2(K_1^\alpha+\epsilon)<1$ para $\epsilon<\epsilon_0$ ✓; Lema 3.3 ✓ |
| Não existência | (3.16)–(3.18) ✓; $f'(t)=t^{-\alpha-2}((p-1)t^{p+\alpha}-(\alpha+1)\lambda)$ ✓; $t_0$, $f(t_0)$ e $\lambda_{ne}=(\lambda_1^2/C)^{(p+\alpha)/(p-1)}$ conferidos ✓ |
| Lema 4.1 | $h$ crescente após $t_*=\|\phi_0\|^{(p+1)/(1-p)}$ ($h'(t_*)=1-p>0$) ✓; (4.7) ✓; análise de $L$ ✓ (ver 2.2) |
| (4.9), Obs. 3 | $k_p^{p-1}=\lambda_1^2$; $\lambda_1k_p=\lambda_1^{-(1+p)/(1-p)}\le\|\phi_0\|^{(1+p)/(1-p)}<K_p$ ✓ |
| Lemas 4.2, 4.3, Prop. 1 | $\varphi_1^p\ge\varphi_1$; cotas $k_p\varphi_1\le T_\epsilon u\le K_p\phi_0$; $\|w\|_\infty\le K_p$; $k_p\varphi_1\in\mathcal S$ ✓ |
| Teorema 1.3 | $\lambda_1^{-2/(1-p)}\varphi_1$ é subsolução ($\varphi_1\le1$); passagem ao limite via Lema 3.1 ✓ |
| Corolário 1 | Hölder com $|\Omega|^{(1+\alpha)/2}$, $|\Omega|^{(1-p)/2}$ ✓; liminf/limsup de $\|u_p\|_2^{-(1-p)}$ ✓; Rayleigh, Euler–Lagrange, simplicidade, convergência forte ✓ |

## 5. Tipografia e referências

- Endereços (p. 22): vírgula dupla em "Viçosa,," e "Minas Gerais,,".
- "Hernandez" → "Hernández" (texto e ref. [11]).
- Ref. [14]: nome completo do autor, as outras usam iniciais. Padronizar.
- Ref. [5] (Chen–Feng, *Anal. Math. Phys.* 14 (2024), "1"): conferir número do artigo/páginas.
- Todas as 15 referências são citadas no texto ✓.
