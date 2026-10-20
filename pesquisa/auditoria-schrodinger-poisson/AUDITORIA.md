# Auditoria — *Schrödinger–Poisson systems with singular terms and general nonlinearities of critical exponential growth*

> **Arquivo temporário.** Esta pasta guarda a auditoria só enquanto o manuscrito
> está sendo revisado; pode ser apagada depois.

Autores: A. L. A. de Araujo, H. R. de Assis, L. F. O. Faria — manuscrito de 28 páginas
(versão "Copy (Copy)").
Data da auditoria: 2026-10-20.

## Resumo

- **Nenhum erro matemático encontrado.** Conferi linha a linha todas as contas
  das Seções 2–5 (lista abaixo). Os Teoremas 1.3 e 1.5 seguem dos argumentos
  apresentados.
- **A Seção 6 foi reproduzida exatamente** com um código independente
  (`reproduz_secao6.py`): 39 iterações, mesmos resíduos e mesmos máximos.
- Há **um ponto de escopo que vale a pena resolver** (a hipótese `N ≥ 3` não
  é usada em lugar nenhum; com `N = 2` o exemplo numérico passaria a ser um caso
  do Teorema 1.5), **uma melhoria fácil** no resultado de unicidade, e alguns
  ajustes menores de redação e referências.

| # | Prioridade | Assunto | Onde |
|---|---|---|---|
| 1 | Média | `N ≥ 3` não é usado; incluir `N = 2` cobriria o exemplo da Seção 6 | §1, §2, §5, §6 |
| 2 | Média | A unicidade vale em toda a classe da Def. 1.2 e só precisa de `r ≥ 1` | Teo. 1.5, §5.5 |
| 3 | Baixa | Corolário 2.3 e Teorema 2.8 nunca são usados | §2 |
| 4 | Baixa | Conferir o enunciado exato de [12, Theorem B.1] | §5, §5.1, §5.4 |
| 5 | Baixa | "Crescimento crítico" é tratado só no regime de norma pequena — vale dizer isso | Resumo, §1, §4 |
| 6 | Baixa | Letra `η` usada para o sinal do acoplamento e para funções-teste | Lema 2.2, Prop. 3.3 |
| 7 | Baixa | Formatação de referências | Referências |
| 8 | Baixa | Disponibilizar o código da Seção 6 | Data availability |

---

## 1. A hipótese `N ≥ 3` não é usada (prioridade média)

`N ≥ 3` aparece no resumo, na §1 (linha após (P)), no início da §2, no início
da §5 e na §6, mas **nenhum passo das demonstrações depende dela**. Todas as
ferramentas usadas valem para `N ≥ 2`:

- Trudinger–Moser (`α_2 = 4π`), Hardy em `W_0^{1,2}`, Lieberman [15],
  princípio do máximo forte / Hopf [18], Picone;
- a desigualdade de monotonicidade `c_N‖u−v‖^N ≤ ∫(|∇u|^{N−2}∇u − |∇v|^{N−2}∇v)·∇(u−v)`
  vale para expoente `≥ 2`, logo para `N = 2` (com igualdade);
- em (40), o expoente `N − 2 + bγ` continua positivo quando `N = 2` (fica `bγ > 0`);
- `r > N − 1` vira `r > 1`, e (18) continua fazendo sentido.

**Consequência.** Para `N = 2` o `N`-Laplaciano é o Laplaciano. O exemplo (59)
tem `N = 2`, `r = 2`, `α = 1`, `γ = 1/2`, e a condição (35) dá

```
(rN' + γ − 1)(1 − γ)/(N − 1) = (4 + 1/2 − 1)·(1/2) = 1.75 > 1.
```

Ou seja: se `N = 2` for permitido, (59) é exatamente um caso de (P2) coberto
pelo Teorema 1.5 (para `λ < λ*`), e a frase "this experiment does not
constitute a numerical verification of the theorems" deixaria de ser necessária
(resta só saber se `λ = 10⁻²` está abaixo de `λ*`, que não é explícito).

**Sugestão.** Trocar `N ≥ 3` por `N ≥ 2` em todo o texto, ou explicar por que
`N ≥ 3` é exigido (por exemplo, motivação física). Se mantiverem `N ≥ 3`, um
referee provavelmente vai perguntar a mesma coisa.

## 2. Unicidade: o resultado pode ser enunciado de forma mais forte (prioridade média)

A demonstração em §5.5 (identidades (54)–(58)) **não usa** que as soluções
estão em `C¹(Ω̄)`, nem a condição (35). Ela usa só:

- a formulação fraca da Definição 1.2 (que já garante que `u − v` e
  `φ_u − φ_v` são funções-teste admissíveis e que todas as integrais
  `∫φ_u u^r`, `∫φ_u v^r`, `∫φ_v u^r`, `∫φ_u u^{r−1}v`, … são finitas);
- a desigualdade de Young `u^{r−1}v ≤ ((r−1)/r)u^r + (1/r)v^r`, que vale para `r ≥ 1`;
- `α = 0`.

**Sugestão.** Enunciar: *para `α = 0`, `r ≥ 1` e todo `λ > 0`, (P2) tem no
máximo uma solução fraca (Definição 1.2) com `u > 0`*. É um enunciado mais
forte e a prova já existe. O texto atual ("in this class", "Since r > 1") pode
ser mantido como está se preferirem não mexer.

## 3. Resultados enunciados e não usados (prioridade baixa)

- **Corolário 2.3** não é citado em nenhum lugar (a §4 usa o Lema 2.2 direto).
- **Teorema 2.8** (comparação tipo Díaz–Saa) não é citado em nenhum lugar; a
  unicidade da §5.5 é feita por outro argumento. Se for removido, a referência
  [11] fica sem citação e também pode sair.
- A última afirmação do Lema 2.2 (`ω_n → ω` em `L^{p'}`) também não é usada.

Conferi as contas do Teorema 2.8 (gradientes de `η₁`, `η₂`, as cotas
`0 ≤ η₁ ≤ u₁`, `η₂ ≤ ‖u₁/u₂‖_∞^{N−1}u₁` e a identidade de Picone) e estão
corretas — a questão é só de economia do texto.

## 4. Citação de [12, Theorem B.1] (prioridade baixa)

O texto parafraseia o resultado como: se `−Δ_N z = h`, `0 ≤ h ≤ A d^{−γ}` e
`0 ≤ z ≤ B d`, então `z ∈ C^{1,β}(Ω̄)` com estimativa dependendo só de
`N, Ω, γ, A, B`. Recomendo conferir no artigo original se:

- as hipóteses são exatamente essas (em particular se `0 ≤ z ≤ Bd` é hipótese
  do teorema ou é usada só na demonstração);
- a dependência dos constantes é a indicada, pois a §5.4 precisa dela para ter
  `‖u_ε‖_{C^{1,β}} ≤ C_λ` uniforme em `ε`.

Os usos no artigo (Lema 5.1 para `V`, §5.4 para `u_ε`) satisfazem as
hipóteses parafraseadas: `0 ≤ φ₀^{−γ} ≤ Cd^{−γ}`, `0 ≤ V ≤ C₁d`, e
`λ^{N−1} ≤ −Δ_N u_ε ≤ C_λ d^{−γ}`, `u_ε ≤ λ^bV ≤ C_λ d`.

## 5. Sobre "crescimento exponencial crítico" (prioridade baixa)

Na §4 o crescimento exponencial é controlado porque a bola de Galerkin tem raio
`ρ ≤ ρ₀` com `8aρ₀^{N'} ≤ α_N` (13)–(14), ou seja, tudo acontece abaixo do
nível de Trudinger–Moser e não há perda de compacidade. Na §5 o mesmo papel é
feito por `u ≤ k₁ = Mλ^b` pequeno. Isso é correto, mas um referee pode
observar que a criticalidade não é realmente enfrentada. Sugiro uma frase no
fim da §1 dizendo explicitamente que os resultados são do tipo "λ pequeno /
norma pequena", e que por isso a condição crítica entra só através da
constante `α_N`. O mesmo vale para "without an Ambrosetti–Rabinowitz
assumption": o método não é variacional, então a ausência dessa hipótese é
natural — vale dizer isso para não parecer que se está vendendo mais do que se
entrega.

## 6. Notação (prioridade baixa)

`η` denota o sinal do acoplamento em (P), mas também é usada como função-teste
na prova do Lema 2.2 (`η ∈ L^{p'}`, `η_j`) e na prova da Prop. 3.3
(`η = (w_λ − u_k)^+`). Sugiro outra letra (`ζ`, `χ`) nesses dois lugares.

## 7. Referências (prioridade baixa)

- [19] Ruiz: "The schrödinger–poisson equation" → "The Schrödinger–Poisson equation".
- Formato inconsistente das páginas de artigos com número de artigo:
  [4] "pp. Paper No. 108535, 8.", [9] "pp. Paper No. 147, 18.",
  [14] "pp. 123647, 21.", [17] "pp. Paper No. 29, 21." — padronizar
  (ex.: "Paper No. 108535, 8 pp.").
- Os dados bibliográficos que conferi de memória (volume, ano, páginas) de
  [1], [2], [6], [13], [15], [16], [18], [20], [21] estão corretos.

## 8. Seção 6 — reprodução numérica (prioridade baixa)

Reproduzi o experimento com o esquema exatamente como descrito (5 pontos,
`n = 100`, `h = 1/101`, `λ = 10⁻²`, `γ = 1/2`, `ω = 0.5`, `u⁰ = 0.1`,
`φ⁰ = 0`, `φ̂^{(k+1)}` calculado com `u^{(k+1)}`, parada em `𝓡_k ≤ 10⁻⁹`):

| Quantidade | Artigo | Reprodução |
|---|---|---|
| Iterações | 39 | 39 |
| `‖R_u‖_∞` | 7.96 × 10⁻¹⁰ | 7.955 × 10⁻¹⁰ |
| `‖R_φ‖_∞` | 1.08 × 10⁻¹⁴ | 1.075 × 10⁻¹⁴ |
| `max u` | 9.15146 × 10⁻³ | 9.151456 × 10⁻³ |
| `max φ` | 3.97852 × 10⁻⁶ | 3.978521 × 10⁻⁶ |
| `e^{u²} − 1` | ≲ 8.4 × 10⁻⁵ | 8.375 × 10⁻⁵ |
| `𝓡_1` (Fig. 3) | ≈ 10³ | 1.020 × 10³ |

Tudo confere. A afirmação de positividade (matriz-M + coeficiente diagonal
não negativo + relaxação convexa) também se confirma: `min u ≈ 1.08 × 10⁻⁴ > 0`,
`min φ ≈ 3.5 × 10⁻⁹ > 0`.

**Sugestão.** A seção "Data availability" poderia apontar para o código (por
exemplo, um repositório público), já que o resultado é totalmente reproduzível.

---

## Apêndice — o que foi conferido e está correto

**Seção 1–2**
- Constante de Trudinger–Moser `α_N = N|S^{N−1}|^{1/(N−1)}` e a cota `C_N|Ω|`.
- Lema 2.2 (Vitali + densidade de `L^∞` em `L^{p'}`).
- Teorema 2.4: expoente `1/(N−1)` compatível com a homogeneidade do `N`-Laplaciano.
- Teorema 2.5: continuidade de `S` via monotonicidade e `|h_n − h|_∞‖u_n − u‖`.
- Lema 2.7 (comparação fraca) e Teorema 2.8 (gradientes de `η₁`, `η₂`, Picone).
- Lema 2.9: Hardy + Hölder com expoentes `N` e `N'`, e o argumento em `Ω_δ`.

**Seção 3**
- Definição (6) de `f_k`: continuidade nos pontos de junção e imparidade.
- Lema 3.1: convergência local uniforme de `f_k` para `f`.
- Lema 3.2: `f_k(s) ≤ L(2s)^r e^{α(2s)^{N'}}` para `s ≥ 1/k` e
  `f_k(s) ≤ L2^r e^{α2^{N'}} k^{1−r}s` para `0 ≤ s ≤ 1/k`; daí (7) e (8).
- Prop. 3.3: escolha de `δ_λ` em (12), subsolução `w_λ = δ_λ z` e comparação
  com `(w_λ − u_k)^+`; cota `(u_k + 1/k)^{−γ} ≤ C_λ d^{−γ}`.

**Seção 4**
- (13)–(14) e as estimativas de Hölder com expoentes `8r`, `8(r−1)`, `1/8`, `4`.
- (17)–(20), em particular `v/(v_+ + 1/k)^γ ≤ v_+^{1−γ}`.
- Expansão de `⟨Φ_{m,k}(ζ,ξ),(ζ,ξ)⟩`, (22), escolha de `ρ`, `k₀` e
  `Λ = ρ^{N−1+γ}/(4C_γ)`, que dá `≥ ρ^N/4`.
- Passagem `m → ∞` com `(S)_+`, sinais de `u` e `φ` via `u_−`, `φ_−`.
- Passagem `k → ∞`: (30)–(34) e o uso do Lema 2.9 na convergência forte.

**Seção 5**
- `(35) ⇒ r > N − 1`, via `(N−1)/t + t − N = (1−t)(N−1−t)/t`.
- Lema 5.1: `−Δ_N F(φ₀) = B − φ₀^{1−γ} + (1−γ)φ₀^{−γ}|∇φ₀|^N`,
  positividade de `H` em `Ω̄` (Hopf na fronteira), e a justificativa fraca
  via funções-teste de suporte compacto.
- Identidades (43)–(45): `k₂ k₁^{r/(N−1)} … = Z k₁^R e^{αN'k₁^{N'}}` usando
  `r + r/(N−1) = rN'` e `1 + 1/(N−1) = N'`; `b(N−1) = 1 − γ` em (50).
- Lema 5.2 (invariância de `𝒜_λ`), Schauder em §5.3, limite `ε → 0` em §5.4.
- §5.5: desigualdade de Young e (57)–(58).

**Seção 6**
- Ver item 8 (reprodução integral).
