"""
Generate and run a Singular script checking

    u_1(1) = ... = u_K(1) = 0   ==>   f_e = 0

for x' = f(t) x^3 + t^(n-1) x^2 on [-1, 1], after eliminating the linear
conditions m_0 = m_1 = m_2 = 0 (Theorem 1.7 of Alves-Lemos-Araujo).

We write f_e = sum_i c(i) e_i(t), where e_i is a basis of the even polynomials
of degree <= 2E with int_{-1}^{1} e_i G^k dt = 0 (k = 0, 1, 2), G = (t^n-1)/n,
and f_o = sum_{j<=R} b(j) t^(2j+1). The script checks that u_1(1), u_2(1),
u_3(1) vanish identically on this family and that, for each i, the ideal
<u_4(1), ..., u_K(1), 1 - z c(i)> is the unit ideal (Rabinowitsch trick),
i.e. c(i) lies in the radical of the ideal of center conditions.

Usage: python3 center_check_reduced.py n E R K
Requires: sympy, Singular.
"""
import os
import subprocess
import sys
import tempfile

import sympy as sp


def kernel_basis(n, E):
    t = sp.Symbol('t')
    G = (t**n - 1) / n
    a = sp.symbols('a0:%d' % (E + 1))
    fe = sum(a[j] * t**(2 * j) for j in range(E + 1))
    eqs = [sp.integrate(fe * G**k, (t, -1, 1)) for k in range(3)]
    sol = sp.solve(eqs, a[:3], dict=True)[0]
    basis = []
    for i in range(3, E + 1):
        e = sp.expand(fe.subs(sol).subs({a[j]: (1 if j == i else 0) for j in range(3, E + 1)}))
        den = sp.ilcm(*[sp.fraction(c)[1] for c in sp.Poly(e, t).coeffs()])
        basis.append(sp.expand(e * den))
    return basis


def singular_script(n, E, R, K):
    basis = kernel_basis(n, E)
    d = len(basis)
    fe = ' + '.join('c(%d)*(%s)' % (i, str(basis[i]).replace('**', '^')) for i in range(d))
    fo = ' + '.join('b(%d)*t^%d' % (j, 2 * j + 1) for j in range(R + 1))
    return f"""
ring r = 0,(z,t,c(0..{d-1}),b(0..{R})),dp;
option(redSB);
proc prim(poly p)
{{
  poly q = 0; int i; matrix M = coeffs(p, t);
  for (i = 1; i <= nrows(M); i++) {{ q = q + M[i,1]*t^i/i; }}
  return(q - subst(q, t, -1));
}}
poly f = {fe} + {fo};
list u; list S; u[1] = prim(t^({n-1})); S[1] = 1;
int k; int i; poly s; ideal I;
for (k = 1; k <= {K}; k++)
{{
  u[k+1] = prim(f*S[k]);
  s = 0; for (i = 0; i <= k-1; i++) {{ s = s + u[i+1]*S[k-i]; }}
  S[k+1] = s;
  if (k >= 4) {{ I = I + subst(u[k+1], t, 1); }}
  else {{ if (subst(u[k+1], t, 1) != 0) {{ print("ERROR: u_k(1) not identically zero for k<=3"); }} }}
}}
int ok = 1; ideal J; int j;
for (j = 1; j <= {d}; j++)
{{
  J = std(I + ideal(1 - z*var(2+j)));
  if (!(size(J) == 1 && J[1] == 1)) {{ ok = 0; }}
}}
print("n={n} deg fe<={2*E} deg fo<={2*R+1} K={K} : fe=0 forced? " + string(ok));
quit;
"""


if __name__ == '__main__':
    n, E, R, K = map(int, sys.argv[1:5])
    with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False) as fh:
        fh.write(singular_script(n, E, R, K))
        path = fh.name
    out = subprocess.run(['Singular', '-q', '--no-warn', path], capture_output=True, text=True)
    os.unlink(path)
    print(out.stdout.strip(), out.stderr.strip())
