# Auditoria — *Small positive solutions for singular prescribed mean curvature equations with convection*

Autores: A. L. A. de Araujo e K. F. L. Madalena · versão auditada: PDF de 17 páginas (outubro/2026).

## Veredito

A estratégia (sub/supersolução na escala $t_\lambda=\lambda^{1/(1+\alpha)}$, truncamento, contração do gradiente e remoção da contração por estimativa $C^{1,\beta}$ uniforme) é sólida e o Teorema 1.1 parece **verdadeiro**. Mas a prova, como está escrita, **não fecha**:

| # | Onde | Problema | Gravidade |
|---|------|----------|-----------|
| G1 | Lema 2.2, (11), (13)–(14), p. 5–6, p. 15 | A desigualdade (11) é **falsa**; a cadeia de estimativas da p. 5 tem um erro algébrico | Grave (corrigível) |
| G2 | Lema 5.4 | Escolha **circular** das constantes: $C_0 \to H \to C \to C_0$ | Grave (corrigível) |
| G3 | Lema 5.4, p. 14 | A barreira inferior $w\ge c_1 d$ **uniforme em $\lambda$** não é provada | Grave (corrigível) |
| G4 | Prop. 4.2–4.3 | A Prop. 4.3 exige $R\ge\|\nabla\overline u_\lambda\|_\infty$, mas o texto diz "todo $R>0$" | Grave (corrigível) |

Nenhum dos quatro exige mudar o enunciado do teorema; as correções estão abaixo.

---

## Problemas graves

### G1. O Lema 2.2 (11) é falso

**Afirmação:** $\tfrac12 \le (s\,a(s^2))'/a(s^2) \le 1$ para $s>0$; em (14), o quociente seria $>5/8$ em $[1,\sqrt2)$.

**Contraexemplo.** Para $v=s^2\in[1,2)$, com $c=7b$,

$$
\frac{(s\,a(s^2))'}{a(s^2)}=\frac{a(v)+2va'(v)}{a(v)}=\frac{5v^2-12v+11}{v^2-4v+11}=:R(v).
$$

$R'(v)=0 \iff v^2-11v+11=0$, logo o mínimo ocorre em $v_*=(11-\sqrt{77})/2\approx1{,}1125$, onde $R(v_*)\approx 0{,}4929<\tfrac12$. Em $v=6/5$: $R=3{,}8/7{,}64\approx0{,}497$. Portanto (11) e (14) são falsas.

**Origem do erro (p. 5).** A linha
$$b(s^2-2)(s^2-2+4s^2)+c \;\ge\; b(s^2-2)\,4b\,s^2+c$$
tem um $b$ a mais: o correto é $4b\,s^2(s^2-2)$. O resto da cadeia ($\ge c-16b^2$) herda o erro. A cota usada para o numerador, $(7-\sqrt2)/(8\sqrt2)\approx0{,}494$, é **maior** que o mínimo real $c-\tfrac{16}{5}b=\tfrac{19}{40\sqrt2}\approx0{,}336$ (atingido em $v=6/5$).

**Impacto.**
- Positividade $a+2va'>0$ (eq. (10)): continua verdadeira — o mínimo é $19/(40\sqrt2)>0$.
- (25) ($\tfrac32\le t^2\phi/\Phi\le2$): numericamente **verdadeira** (ínfimo $\approx1{,}650$), mas a prova cita (13)/(11) e, portanto, está sem justificativa.
- p. 15: "o autovalor radial está em $[\tfrac12 a, a]$" — falso pelo mesmo motivo.
- Constante de elipticidade $\nu$ em (41): continua uniforme em $\lambda$, mas deve ser declarada corretamente.

**Correção sugerida.** Em $[1,2)$: numerador $=b(5v^2-12v+11)\ge \tfrac{19}{5}b$ (mínimo em $v=6/5$) e denominador $\le b+c=8b$. Logo
$$\frac{19}{40}\le\frac{(s\,a(s^2))'}{a(s^2)}\le1,\qquad s>0,$$
e o autovalor radial satisfaz $a+2va'\ge \nu:=\tfrac{19}{40\sqrt2}$. Trocar $\tfrac12$ por $\tfrac{19}{40}$ em (11), (19); trocar $\tfrac58$ por $\tfrac{19}{40}$ em (14); em (25) usar $1+\tfrac{19}{40}=\tfrac{59}{40}$ (suficiente para $\Delta_2$), ou provar diretamente o valor $3/2$. Na p. 15, usar $[\nu, a]$.

Script de verificação: [`verifica_lema22.py`](verifica_lema22.py).

### G2. Circularidade na escolha de $C_0$ (Lema 5.4)

No Lema 5.4, $R_\lambda=C_0t_\lambda$ e
$$t^{-1}f_{R_\lambda}(x,tw,t\nabla w)\le a\,t^{r_1-1}c_2^{r_1}+b\,C_0^{r_2}t^{r_2-1}=:H.$$
A constante $C$ dos Lemas 5.1–5.2 depende de $H$, que depende de $C_0$. A frase final "Increasing $C_0$ once and for all, if necessary" altera $H$ e, portanto, $C$ — o argumento é circular.

**Correção (ordem não circular):**
1. Fixe $H:=1$.
2. Com $c_1,c_2$ uniformes (ver G3), $\nu,L$ (G1) e $H=1$, os Lemas 5.1–5.2 dão $C_*$, $\beta$ independentes de $\lambda$ **e de $C_0$**.
3. Defina $C_0:=\max\{C_*,G\}$ ($G=2R$ de (30); ver G4).
4. Escolha $\lambda_0\le\Lambda$ tal que $a\,t_\lambda^{r_1-1}c_2^{r_1}+bC_0^{r_2}t_\lambda^{r_2-1}\le1$ e $C_0t_\lambda<1$ para $\lambda\le\lambda_0$ (possível pois $r_1,r_2>1$).
5. Então $\|w\|_{C^{1,\beta}}\le C_*$, logo $\|u\|_{C^1}\le C_*t_\lambda\le C_0t_\lambda=R_\lambda$.

### G3. Barreira inferior uniforme em $\lambda$ não provada

A p. 14 diz que "as cotas de ordem da Prop. 4.3, divididas por $t_\lambda$, dão barreiras independentes de $\lambda$". A superior sim: $\overline u_\lambda/t_\lambda=\beta\le B$. A inferior **não**: o Lema 2.5 dá $u_{0,\lambda}\ge C_1 d$ com $C_1=C_1(\lambda)$, e $c_1$ precisa ser independente de $\lambda$ (é usado em (43) e em $C_F$ no Lema 5.2).

**Correção.** $v_\lambda:=u_{0,\lambda}/t_\lambda$ resolve $-\operatorname{div}\mathcal A_\lambda(\nabla v)=v^{-\alpha}$. Seja $\psi=\varphi_1$ (autofunção principal do Laplaciano, $\psi\ge k\,d$). Como $|D\mathcal A_\lambda|\le L$,
$$-\operatorname{div}\mathcal A_\lambda(\nabla(\varepsilon\psi))\le \varepsilon\,L\,N\|D^2\psi\|_\infty\le(\varepsilon\|\psi\|_\infty)^{-\alpha}\le(\varepsilon\psi)^{-\alpha}$$
para $\varepsilon\le\varepsilon_0$ independente de $\lambda$. Testando com $(\varepsilon\psi-v_\lambda)^+$ (como na Prop. 4.3; $s\mapsto s^{-\alpha}$ decrescente) obtém-se $v_\lambda\ge\varepsilon_0 k\,d=:c_1d$.

### G4. A Prop. 4.3 não vale para todo $R>0$

No conjunto $\{u>\overline u_\lambda\}$ a equação (39) tem $f_R(x,\overline u_\lambda,\nabla\overline u_\lambda)=f(x,\overline u_\lambda,\pi_R\nabla\overline u_\lambda)$, mas a desigualdade de supersolução usa $f(x,\overline u_\lambda,\nabla\overline u_\lambda)$. Como não há monotonia de $f$ em $\xi$, a subtração só funciona se $\pi_R(\nabla\overline u_\lambda)=\nabla\overline u_\lambda$, isto é, $R\ge\|\nabla\overline u_\lambda\|_\infty$ ($\le Gt_\lambda$).

**Correção.** Na Prop. 4.3, supor $R\ge Gt_\lambda$ (por exemplo $R\ge1$, já que $Gt_\lambda\le1$). No Lema 5.4, exigir $C_0\ge G$ (já incluído em G2, passo 3). O lado da subsolução não precisa disso (usa só $f\ge0$).

---

## Problemas moderados

- **M1. Lema 2.7 (princípio de comparação) está corrompido.** "$0<q<1$" — $q$ não é definido; "$t\mapsto t^{-\alpha}/t^{1/2}$" não corresponde a nenhuma hipótese visível. Sugestão: remover o lema e usar o argumento direto com $(u_{0,\lambda}-\overline u_\lambda)^+$, que o próprio texto já menciona e que é válido (suporte compacto, pois $\overline u_\lambda\ge t_\lambda>0$ em $\overline\Omega$). Há também uma linha solta "." após o lema.
- **M2. Pseudomonotonia (Prop. 4.2).** Para $\int(T_\lambda u_n)^{-\alpha}(u_n-u)\to0$ a limitação do Lema 4.1 não basta; falta o argumento de compacidade: dividir em $\{d<\delta\}$ (Hardy: $\le C\delta^{1-\alpha}\|\nabla(u_n-u)\|_2$) e $\{d\ge\delta\}$ ($\le\delta^{-\alpha}\|u_n-u\|_1\to0$). Para o termo de convecção não é "convergência dominada": basta $|f_R|\le$ const. e $u_n\to u$ em $L^2$.
- **M3. Regularidade no Lema 2.5.** "Standard regularity for singular uniformly elliptic equations" sem referência; citar [19, Thm B.1] (como na Seção 5). Também vale conferir **as hipóteses exatas** de [19, Thm B.1] (forma do operador, exigência $0\le u\le Cd$, dependência das constantes), pois todo o Lema 5.2 depende disso. Não tive acesso ao texto de [19] nesta auditoria.
- **M4. Lema 5.1, p. 13.** "if necessary, increasing $K$ once at the initial step" é desnecessário (as três quartas partes já bastam) e, se feito, invalidaria (47)–(48), que dependem de $K$. Remover.
- **M5. Dimensão.** O texto assume $N>2$, mas o argumento só usa $\delta_*=N/\sqrt2-1>0$, isto é, $N\ge2$. Ou enunciar para $N\ge2$, ou justificar a restrição.
- **M6. Conflitos de notação** (alguns tornam fórmulas ambíguas):
  - $a,b$: constantes em (2) **e** função $a(\cdot)$ e $b=1/(8\sqrt2)$ em (4);
  - $\alpha$: expoente **e** subsolução na Def. 2.1 (aparece $\lambda/\alpha^\alpha$);
  - $\beta$: expoente de Hölder **e** supersolução (Def. 2.1, (30));
  - $\lambda$: parâmetro **e** variável da norma de Luxemburg / (22);
  - $R$: raio em (30), raio de truncamento $\pi_R$, e $R_\lambda$;
  - $M,B$: Seção 3 ($M=R^2+1$, $B=M-\rho^2$) **e** Lema 5.1 ($\|w\|_\infty\le M$, $q=As-Bs^{2-\alpha}$);
  - $\rho$, $C_1$, $\delta$ também reutilizados (Seção 3 vs. 5; Lema 2.5 vs. (45); p. 14).
- **M7.** Após (32): "If $a=0$ or $b=0$ the entry is omitted", mas (2) exige $a,b>0$. Ou permitir $a,b\ge0$ em (2), ou remover a frase.
- **M8.** O Lema 5.4 fala em "*the* solution of (40)"; unicidade não foi provada — usar "a solution supplied by Proposition 4.2".

## Pequenos problemas de texto

- p. 3: "Lemma **??**" (referência quebrada).
- p. 13: "HÃ¶lder"; p. 17: "ItajubÃ¡" (codificação UTF-8).
- p. 2: "...gradient-dependent nonlinearity. Is the main issue addressed here." (frase fragmentada).
- p. 4: "can to be obtained"; p. 6: "Moroever", "called as Luxemburg norm"; p. 8: "follwing"; p. 10: "from the a priori".
- p. 7: "(25) follows from (13)" — deveria citar (11)/(19) (e ver G1).
- Resquícios de revisão: "We now **correct** the supersolution construction" (p. 8), "without any circular dependence" (p. 12), "coupled choice (49)" (p. 12).
- Notação inconsistente: $\|\nabla\beta\|^2$ (Def. 2.1) vs. $|\nabla u|^2$; "$div$" em itálico em (3), (28), (39) vs. romano em outras.
- Lema 2.2: o caso $s\in[2,\infty)$ é redundante com $[\sqrt2,2)$ (em ambos $a=c$).
- A subseção de Orlicz–Sobolev (pp. 6–7) só é usada para $W_0^{1,\Phi}=H_0^1$; pode ser bem encurtada.

## Referências

- **Não citadas no texto (11):** [8], [16], [20], [21], [25], [26], [28], [30], [32], [33], [34]. Em particular [8] (Brenier, transporte ótimo) não tem relação com o artigo. Citar ou remover.
- [11] Figueiredo–Rădulescu: volume "**1923** (2021)" está errado (Ann. Mat. Pura Appl. 2021 é vol. 200); conferir volume e páginas.
- [9] "Moutain" → "Mountain".
- [23] "Lopéz" → "López".
- [28] "...quasilinear elliptic problem **an** oscillatory potential" → "**with** an oscillatory potential".
- [17] Kufner–John–Fučík, *Function Spaces*: conferir editora/ano (edição original Noordhoff/Academia, 1977).
- [12], [14], [17], [29]: formatação diferente das demais (títulos sem itálico).

## O que foi verificado e está correto

- Escala: $\lambda(t_\lambda w)^{-\alpha}=t_\lambda w^{-\alpha}$.
- $a\in C^1([0,\infty))$: continuidade e derivada coincidem em $v=1$ ($a=1/\sqrt2$, $a'=-1/(4\sqrt2)$) e em $v=2$; $c\le a\le1$.
- $a+2va'>0$ (com a correção de G1); $D\mathcal A_\lambda$ tem autovalores em $[19/(40\sqrt2),1]$, uniformes em $\lambda$; logo (41) vale com $\nu=19/(40\sqrt2)$, $L=1$.
- (21)–(23); (25) numericamente (ínfimo $\approx1{,}650$); condição $\Delta_2$.
- Cálculo (35)–(36), (37)–(38) e a escolha de $\Lambda$ em (32): $\overline u_\lambda=t_\lambda\beta$ é de fato supersolução e $|\nabla\overline u_\lambda|\le1$.
- Lema 4.1 (Hardy) e o lado da subsolução na Prop. 4.3.
- Lema 5.1: (45), (50)–(53) e a escolha (46)–(49) de $K,\rho,A,B$; a comparação final.
- Reescalonamento (59) e a cota de $t^{-1}f_{R_\lambda}$.
- Seção 6: com $C_0t_\lambda<1$, $a(|\nabla u|^2)=(1+|\nabla u|^2)^{-1/2}$ e (3) coincide com (1).

---

## Correções implementadas em `mean-curvature-singular_AMS.tex`

- **G1:** Lema 2.2 reescrito com $\tfrac{19}{40}$ e a nova cota (12) $\nu_0=\tfrac{19}{40\sqrt2}$ para o autovalor radial; cadeia da p. 5 corrigida; (14)/(16) corrigida; Observação 2.3 (o $\tfrac12$ não vale); Observação 2.4 (autovalores de $D\mathcal A$, monotonia forte (19)); (29) agora com prova e cota $\tfrac{59}{40}$.
- **G2:** Lema 5.5 reescrito em 5 passos, com ordem não circular ($H=1$ → $C_\sharp$ → $C_0=\max\{C_\sharp,G\}$ → $\lambda_0$).
- **G3:** novo Lema 5.4 (barreira inferior uniforme $u_{0,\lambda}\ge c_1t_\lambda d$, via $\varepsilon_0\varphi_1$).
- **G4:** Prop. 4.3 agora supõe $R\ge Gt_\lambda$; (38) inclui $\|\nabla\overline u_\lambda\|_\infty\le Gt_\lambda$.
- **M1–M8:** princípio de comparação (Lema 2.9) reescrito com prova curta; argumento de compacidade do termo singular na Prop. 4.2; citações de regularidade no Lema 2.7; "aumentar $K$" removido do Lema 5.1; $N\ge2$; notação renomeada ($k_1,k_2$; $\overline\psi,\underline\psi$; $\vartheta,R_\Omega,r_\Omega,M_0,B_0$; $\mu$ na norma de Luxemburg; $E$ no Lema 5.1; $m_1,m_2$ no Lema 2.7); frase sobre $a=0$ ou $b=0$ removida.
- **Texto:** "Lemma ??", acentos (inputenc utf8), erros de digitação e resquícios de revisão.
- **Referências:** removidas as não citadas e sem relação com o artigo — Brenier, Kesavan, Stampacchia–Lewy, Omari–Zanolin, Tan–Fang. As não citadas mas relacionadas passaram a ser citadas: del Pino e Gui–Lin (regularidade singular), Montenegro (princípio do máximo forte), Ripoll, Tersenov e Lorca–Montenegro 2012 (estimativas de altura/gradiente). Corrigidas [FR] (vol. 200, 2213–2233), Mountain, López, Kufner–John–Fučík.
