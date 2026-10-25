# Generalizações do critério de centro de Alves–Lemos–Araujo (conjectura de Lijun–Yun)

Notas de trabalho que organizam as generalizações de

> A. M. Alves, A. Lemos, A. L. A. de Araujo, *A partial positive answer to a Lijun–Yun conjecture*,
> Boll. Unione Mat. Ital. 13 (2020) 49–59.

para a equação de Abel `x' = f(t) x^3 + t^(n-1) x^2`, `t ∈ [-1,1]`.

- `generalizacoes_lijun_yun.tex` / `.pdf` — o texto.
- `codigo/centro_abel.py` — verificações (identidades, determinante, exemplo de Giné, bases de Gröbner).

Requer `sympy` e `mpmath`. Exemplos:

```
python3 codigo/centro_abel.py identidades
python3 codigo/centro_abel.py gine
python3 codigo/centro_abel.py determinante
python3 codigo/centro_abel.py hiperplano 2
python3 codigo/centro_abel.py groebner 2 3 3 9    # n=2, deg f_e<=6, deg f_o<=7, K=9
```
