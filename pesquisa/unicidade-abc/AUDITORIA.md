# Auditoria — "On the Ambrosetti–Brezis–Cerami conjecture"

**Autores:** Anderson L.A. de Araujo, Luiz F.O. Faria
**Arquivo auditado:** `abc_unicidade.pdf` (26 páginas; removido do repositório — a fonte está em `abc_unicidade_original.tex`)
**Data da auditoria:** 06/11/2026
**Status:** armazenamento temporário (rascunho de trabalho)

---

## 1. Resumo do veredito

A estrutura da prova está **correta** e as contas principais **conferem**. Não encontrei erro que invalide o Teorema 1.1. A estratégia é:

1. Escala `u(x) = λ^{1/(p-q)} U(λ^γ x)` → problema sem parâmetro na bola `B_R`, `R = λ^γ`.
2. Amplitudes pequenas: limite côncavo puro + unicidade tipo Brezis–Oswald (Picone) abaixo de `A* = ((1-q)/(p-1))^{1/(p-q)}`.
3. Amplitudes intermediárias: cota inferior elementar para o primeiro zero `R(A)`.
4. Amplitudes grandes: variáveis de Emden–Fowler, estabilidade exponencial do equilíbrio `(ℓ,0)` (onde entra `p > p_c`) e argumento de aprisionamento robusto.
5. Simetria (moving planes para `f = Lipschitz + não decrescente`) → unicidade entre todas as soluções clássicas positivas.

Há, porém, **duas lacunas pequenas de rigor** (seção 3), **frases condicionais que enfraquecem o texto** (seção 4) e vários **problemas editoriais** (seções 5–6) que convém corrigir antes da submissão.

---

## 2. Verificações feitas (todas conferem)

| Item | Onde | Resultado |
|---|---|---|
| Escala (12)–(14): expoentes `1/(p-q) + 2γ = p/(p-q)`, `1 + q/(p-q) = p/(p-q)` | Seção 2 | ✅ |
| Monotonia (17) dos perfis de tiro | Seção 2 | ✅ |
| Derivada de `s^{q-1}+s^{p-1}` e o limiar `A*` | Lema 3.1 | ✅ |
| Identidade de Picone com `v+ε`, Fatou, soma (21)+(22), conclusão `u ≡ v` | Lema 3.2 | ✅ (argumento padrão de Brezis–Oswald) |
| Cota `ρ_q ≤ √(2N/(1-q))` e transversalidade `V'(ρ_q) < 0` | Lema 3.4 | ✅ |
| Reescala `U_A(r) = A V_A(A^{-(1-q)/2} r)` e equação (26) | Prop. 3.5 | ✅ |
| `liminf ρ_A ≥ ρ_q` (dependência contínua) e `limsup ρ_A ≤ ρ_q` (monotonia de `s^{N-1}V_A'`) | Prop. 3.5 | ✅ |
| Transformação de Emden–Fowler: `A = N-2-2a`, `B = a(N-2-a)`, `2 - ap = -a` | (28)–(29) | ✅ |
| Termo perturbativo `ε e^{mt} x^q` com `m = 2 + a(1-q) = 2(p-q)/(p-1)` | (35)–(36) | ✅ |
| `p > p_c ⇔ a < (N-2)/2 ⇔ A > 0` | Seção 4 / Obs. 8.1 | ✅ |
| Energia negativa, limitação da órbita, LaSalle, exclusão de `(0,0)` pela energia | Lema 4.1 | ✅ |
| Linearização `z'' + Az' + (p-1)Bz = 0` (usa `ℓ^{p-1} = B`) | Seção 4 | ✅ |
| Estimativa de variação de constantes e tempo de saída (`sup < 2η/3`) | Lema 4.2 | ✅ |
| `T_ε = -(1/m) log ε + τ_0`, `ε e^{mT_ε} = e^{mτ_0} < δ/2`, `S(ε) > e^{τ_0} ε^{-(p-1)/(2(p-q))}` | Lema 5.1 | ✅ |
| Reescala `U_A(r) = A v_ε(A^{(p-1)/2} r)`, `ε = A^{q-p}`, cancelamento exato em (37) | Prop. 5.2 / Obs. 8.2 | ✅ |
| `c_{0,1} = (N A_0/(A_1^q + A_1^p))^{1/2}` e contradição `0 > A_0/2` | Lema 5.3 | ✅ |
| Combinação final e escolha `R* < min{R_0, c_{0,1}, c_*}` | Teo. 6.1 | ✅ |
| Assintótica: `(1-q)+(p-1) = p-q`, logo expoente `1/(1-q)`; constante `ρ_q^{-2/(1-q)}` | Prova do Teo. 1.1 | ✅ (confere com o caso côncavo puro `-Δu = λu^q`) |

---

## 3. Lacunas de rigor (corrigir)

### 3.1 Continuidade de `A ↦ R(A)` (Corolário 3.6) — **média**
O texto diz que transversalidade + dependência contínua "implicam persistência e continuidade do primeiro zero". Isso exige comparar `U_{A'}` com `U_A` **um pouco além** de `R(A)`, onde `U < 0` e o termo `U^q` não está definido (nem é Lipschitz em 0).

**Sugestão:** evitar a extensão. Use exatamente os dois argumentos da Prop. 3.5:
- semicontinuidade inferior por dependência contínua em compactos de `[0, R(A))`;
- semicontinuidade superior pela monotonia de `r^{N-1}U'_A` (como no passo `limsup ρ_A ≤ ρ_q`).

Observação: por (17), `U_A'(R(A)) < 0` para **todo** `A` com `R(A) < ∞`, então a continuidade vale em todo o conjunto `{R(A) < ∞}`, não só em `(0, δ)`. Alternativa ainda mais simples para a existência: citar a solução minimal de Ambrosetti–Brezis–Cerami [1].

### 3.2 Classe de regularidade na Prop. 3.3 — **baixa**
O Lema 3.2 exige `u, v ∈ H¹₀(Ω) ∩ C(Ω̄)`, mas a Prop. 3.3 fala de "solução positiva" sem especificar a classe. Uma solução clássica `C²(B_R) ∩ C(B̄_R)` não está automaticamente em `H¹₀`.
**Sugestão:** enunciar a Prop. 3.3 para soluções em `H¹₀ ∩ C(B̄_R)`, ou observar que ela só é aplicada a perfis radiais de tiro, que são `C¹` até a fronteira.

### 3.3 Afirmação adiantada no Teorema 6.1 — **baixa**
"Em particular, pela Prop. 3.5, o parâmetro de tiro tende a zero ao longo do ramo pequeno quando `R → 0+`." Isso só é **provado** depois (na prova do Teo. 1.1, usando a exclusão `A ≥ A_0`). Sugiro remover a frase ou movê-la.

### 3.4 Teorema de simetria (Seção 7) — **verificar referência**
O argumento depende inteiramente de [5] (Moreira dos Santos–Nornberg), que trata **sistemas totalmente não lineares** em solução de viscosidade. Recomendo:
- citar o **enunciado exato** (número do teorema e hipóteses) de [5] e conferir hipótese por hipótese (domínio, classe da solução, decomposição `f₁ + f₂`, sinal exigido em `f₂`);
- verificar se Gidas–Ni–Nirenberg [2] já cobre o caso `f = Lipschitz + não decrescente` — se sim, é a referência mais clássica e natural (e [2] já está na bibliografia, mas não é citado).

---

## 4. Frases condicionais que enfraquecem o texto

Em três lugares o texto soa como se a unicidade fosse uma **hipótese**, quando ela já foi provada na Prop. 3.3:

- Cor. 3.6 (enunciado): *"Proposition 3.3 applies to all positive solutions of (13), then this solution is unique…"* → reescrever: *"By Proposition 3.3, this solution is unique among all positive solutions satisfying ‖U‖ < A*."*
- Cor. 3.6 (prova): *"Finally, if Proposition 3.3 gives uniqueness…"* → *"Finally, Proposition 3.3 gives uniqueness…"*
- Teo. 6.1 (prova): *"Assuming, as in Proposition 3.3, uniqueness in the class…"* → *"By Proposition 3.3, there is at most one…"*

Um revisor pode ler isso como sinal de que a prova está incompleta.

---

## 5. Conflitos de notação

| Símbolo | Usos conflitantes |
|---|---|
| `A` | amplitude `U_A(0) = A` **e** coeficiente de amortecimento `A = N-2-2a` |
| `A_1` | variável local na prova do Cor. 3.6 (`A_1 < A_0`) **vs** limiar de amplitudes grandes na Prop. 5.2 / Teo. 6.1. Só legibilidade, sem conflito lógico (o `A_0` compartilhado entre Cor. 3.6 e Teo. 6.1 é intencional) |
| `λ*` | limiar do Teo. 1.1 **e** parâmetro da solução singular de Zhao–Zhong (introdução) — significados diferentes |
| `ε` | parâmetro da perturbação (Seção 5) **e** variável auxiliar no Lema 3.2 e na Prop. 3.5 |
| `η`, `δ` | Prop. 3.5 **e** Lema 4.2 / Cor. 3.6 |
| `c_0` | Prop. 3.5 **vs** `c_{0,1}` no Lema 5.3 |

**Sugestão:** renomear o amortecimento para `α` (ou `κ`), a amplitude grande para `A_∞`, o limiar do teorema para `λ_0`. (No .tex o termo do Lema 4.2 já é `\mathcal N(Z)`; a confusão com `N` vinha só da extração de texto do PDF.)

---

## 6. Problemas editoriais

- **Referências não citadas:** [2] Gidas–Ni–Nirenberg e [4] Ni–Nussbaum estão na bibliografia mas nunca são citados no texto. Citar ou remover.
- **Conflito com Ma [3]:** o artigo prova unicidade, e [3] afirma multiplicidade radial incluindo o caso supercrítico. Vale dizer explicitamente que o Teo. 1.1 **contradiz** [3] para `λ` pequeno (e, se possível, indicar onde o argumento de [3] falha).
- **Lema 4.1:** o resultado é clássico (comportamento `r^{2/(p-1)} v_0(r) → ℓ` no caso supercrítico). Sugiro citar a literatura clássica (por exemplo Joseph–Lundgren; Budd–Norbury; Merle–Peletier), mantendo a prova por completude.
- **Lema 4.1:** "Consequently," repetido em duas frases seguidas.
- **Marcas de revisão em vermelho:** há trechos em vermelho em várias páginas (por exemplo pp. 22–24: prova do Teo. 6.1, Seção 7, prova do Teo. 1.1). Remover a cor antes de submeter.
- **Teo. 6.1 (prova):** a desigualdade `R* < min{c_{0,1}, c_*}` é repetida logo após `R* < min{R_0, c_{0,1}, c_*}`; a segunda é redundante.
- **p. 24:** espaço solto no início de " Since u_λ is radial…".
- **Lema 4.2:** observar que `K ≥ 1`, então `η_0 < η/(2K)` já implica `η_0 < η` (a frase "We may also assume η_0 < η" fica redundante).
- **Lema 3.2:** o lema é a desigualdade de Brezis–Oswald; vale citar Brezis–Oswald (1986).
- **Título:** ABC formularam uma *questão em aberto*; "conjecture" segue o título de Ma [3], mas convém confirmar a redação usada em [1].
- **"Acknowledgment of scope":** título incomum; sugiro "Concluding remark" ou "Scope of the result".

---

## 7. Prioridades

1. Fechar a continuidade de `R(A)` (3.1).
2. Reescrever as frases condicionais (seção 4).
3. Conferir as hipóteses exatas do teorema de simetria de [5] (3.4).
4. Resolver os conflitos de notação (seção 5).
5. Ajustes editoriais (seção 6).

---

## 8. Revisão implementada no .tex

Arquivos:
- `abc_unicidade_original.tex` — versão recebida, sem alterações.
- `abc_unicidade_revisado.tex` — versão com as correções marcadas (o PDF compilado não fica no repositório).

Marcação (legenda no topo da 1ª página):
- **azul** = texto acrescentado; **vermelho riscado** = texto removido;
- **[Removed: …]** em vermelho = trecho removido que contém fórmulas;
- **verde [Note to authors: …]** = pendências que só os autores podem resolver;
- as antigas marcas vermelhas (`\textcolor{red}`) viraram `\oldrev{}`, sem cor (conteúdo mantido).

O que foi feito:
1. **Novo Lema 3.6** (continuidade de `A ↦ R(A)` no conjunto aberto `{R(A) < ∞}`), provado sem estender a equação além do primeiro zero. A prova do Cor. 3.6 passou a usá-lo. *Atenção: a numeração mudou — o antigo Cor. 3.6 agora é Cor. 3.7.*
2. **Prop. 3.3** enunciada na classe `H¹₀ ∩ C(B̄_R)`, com observação de que os perfis radiais de tiro estão nessa classe.
3. **Frases condicionais** do Cor. 3.6 e do Teo. 6.1 reescritas como afirmações.
4. **Teo. 6.1**: removida a frase adiantada sobre Prop. 3.5 e a desigualdade redundante `R* < min{c_{0,1}, c_*}`.
5. **Prova do Teo. 1.1**: o passo `U_R(0) → 0` agora tem justificativa explícita (Lema 5.3 em `[a, A_1]` + Prop. 5.2, para todo `a > 0`). Antes dizia só "pela exclusão de `A ≥ A_0`", o que dá `A < A_0`, mas não `A → 0`.
6. **Lema 4.2**: "We may also assume η₀ < η" trocado pela observação `K ≥ 1`.
7. **Notação**: amortecimento `A → α`; parâmetro de Zhao–Zhong `λ_* → λ_ZZ`; constantes locais da Seção 3 renomeadas (`θ`, `σ`, `δ₁`, `η₁`, `c̄`, `Ã`).
8. **Citações**: [2] GNN e [4] Ni–Nussbaum agora citados; novas referências Brezis–Oswald (Lema 3.2) e Joseph–Lundgren (Lema 4.1); frase sobre a incompatibilidade com Ma [3].
9. **Edição**: "Consequently" repetido, espaço solto antes de "Since u_λ…", título "Acknowledgment of scope" → "Scope of the result".

Pendências para os autores (notas verdes no PDF):
- número exato do teorema de [5] e conferência das hipóteses; ver se GNN [2] já cobre "Lipschitz + não decrescente";
- confirmar a referência clássica para o Lema 4.1 (e os dados bibliográficos novos);
- se possível, indicar onde falha o argumento de Ma [3].
