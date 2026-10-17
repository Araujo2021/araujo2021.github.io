# O paradoxo dos pesticidas com difusão espacial e aplicação de suporte compacto

Generalização de Y. C. Li e Y. Yang, *On the paradox of pesticides*,
Commun. Nonlinear Sci. Numer. Simul. 29 (2015) 179–187.

- `paradoxo_pesticidas_espacial.tex` / `.pdf`: texto com modelo, teoremas, provas e simulações.
- `codigo/pest.py`: solver (diferenças finitas + BDF, Neumann, pulso de suporte compacto e limite impulsivo).
- `codigo/figuras.py`: gera as figuras em `figuras/` (`python3 figuras.py` ou `python3 figuras.py fig3`).

Requisitos: Python 3 com numpy, scipy e matplotlib; pdflatex para compilar o texto.
