# O paradoxo dos pesticidas com difusão espacial

Generalizações de Y. C. Li e Y. Yang, *On the paradox of pesticides*,
Commun. Nonlinear Sci. Numer. Simul. 29 (2015) 179–187.

## Versão planar com refúgio (principal)
- `paradoxo_pesticidas_planar.tex` / `.pdf`: Ω ⊂ ℝ², pesticida aplicado só na zona ω
  (termo α ρ_ε(t) χ_ω(x) u), refúgio Ω \ ω. Inclui existência/positividade, limite impulsivo
  (ordem √ε em L¹, ótima, com camada limite e mortalidade de borda), identidade de energia,
  critério nas médias em ω, cota pós-aplicação e zona ótima (princípio da banheira).
- `codigo/pest2d.py`: solver de volumes finitos 2D (Neumann, χ_ω por fração de área, BDF),
  limite impulsivo e constante de borda c*(α, φ).
- `codigo/figuras2d.py`: figuras `figuras/p3_*.pdf` a `figuras/p7_*.pdf`
  (`python3 figuras2d.py f3`, etc.).

## Versão anterior (Ω ⊂ ℝⁿ, perfil espacial suave)
- `paradoxo_pesticidas_espacial.tex` / `.pdf`, `codigo/pest.py`, `codigo/figuras.py`
  (figuras `figuras/fig1_*.pdf` a `figuras/fig6_*.pdf`; `fig1` e `fig2` também são usadas na versão planar).

Requisitos: Python 3 com numpy, scipy e matplotlib; pdflatex para compilar os textos.
