# Generalizações do critério de centro de Alves–Lemos–Araujo (conjectura de Lijun–Yun)

Notas e artigo sobre a equação de Abel `x' = f(t) x^3 + t^(n-1) x^2`, `t ∈ [-1,1]`, a partir de

> A. M. Alves, A. Lemos, A. L. A. de Araujo, *A partial positive answer to a Lijun–Yun conjecture*,
> Boll. Unione Mat. Ital. 13 (2020) 49–59.

- `generalizacoes_lijun_yun.tex` / `.pdf` — notas de trabalho (português).
- `artigo/lijun_yun_higher_order.tex` / `.pdf` — versão em artigo (inglês).
- `codigo/centro_abel.py` — verificações em SymPy (identidades, determinante, exemplo de Giné, Gröbner).
- `codigo/center_check.sing` — implementação independente da recorrência e do teste de Gröbner em Singular.
- `codigo/center_check_reduced.py` — gera e roda o teste em Singular após eliminar as condições lineares
  `m_0 = m_1 = m_2 = 0` (mais rápido; usado na Tabela 1 do artigo).

Requer `sympy`, `mpmath` e `Singular`. Exemplos:

```
python3 codigo/centro_abel.py identidades
python3 codigo/centro_abel.py gine
python3 codigo/centro_abel.py determinante
python3 codigo/centro_abel.py hiperplano 2
python3 codigo/centro_abel.py groebner 2 3 3 9        # SymPy: n=2, deg f_e<=6, deg f_o<=7, K=9
python3 codigo/center_check_reduced.py 2 3 3 9        # Singular, mesmo caso
Singular -q -c 'int n=2; int E=3; int R=3; int K=9;' codigo/center_check.sing
```
