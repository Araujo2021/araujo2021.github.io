"""
Verificacoes computacionais para a nota
"Generalizacoes do criterio de centro para x' = f(t) x^3 + t^(n-1) x^2".

Uso:
  python3 centro_abel.py identidades        # confere u_4(1), u_5(1), u_6(1)
  python3 centro_abel.py gine               # exemplo de Gine et al. (m_3 = pi/2)
  python3 centro_abel.py determinante       # formula de Cauchy para det M_N
  python3 centro_abel.py hiperplano 2       # condicoes u_4, u_5 para deg f_e <= 6 (n dado)
  python3 centro_abel.py groebner n E R K   # centro => f impar? (deg f_e<=2E, deg f_o<=2R+1,
                                            #  usando u_1(1)=...=u_K(1)=0)

Requer: sympy, mpmath.

Esquema: x = rho/(1 - rho u), u(-1)=0, u' = g + rho f/(1-rho u).
Escrevendo u = sum u_k rho^k e S = 1/(1-rho u) = sum S_k rho^k:
  u_0 = G,  u_k = int_{-1}^t f S_{k-1},  S_0 = 1,  S_k = sum_{i<k} u_i S_{k-1-i}.
O centro em x=0 equivale a u_k(1) = 0 para todo k >= 0.
"""
import sys, time
import sympy as sp

t = sp.Symbol('t')


def P(e, dom):
    return sp.Poly(e, t, domain=dom)


def prim(p):
    """primitiva que se anula em t = -1"""
    q = p.integrate()
    return q - q.eval(-1)


def I(p):
    """integral em [-1, 1]"""
    q = p.integrate()
    return q.eval(1) - q.eval(-1)


def ucoeffs(f, g, K):
    u = [prim(g)]
    S = [f.one]
    for k in range(1, K + 1):
        u.append(prim(f * S[k - 1]))
        S.append(sum((u[i] * S[k - 1 - i] for i in range(k)), f.zero))
    return u


def identidades():
    for n in (2, 4, 6):
        c = sp.symbols('c0:8')
        dom = sp.QQ[c]
        f = P(sum(c[i] * t**i for i in range(8)), dom)
        g = P(t**(n - 1), dom)
        G, F = prim(g), prim(f)
        H = prim(f * G)
        m = [I(f * G**k) for k in range(6)]
        u = ucoeffs(f, g, 6)
        sol = sp.solve([m[0], m[1], m[2]], [c[0], c[2], c[4]], dict=True)[0]
        r = sp.Rational
        d4 = (u[4].eval(1) - (m[3] - r(1, 2) * I(g * F**2))).subs(sol)
        d5 = (u[5].eval(1) - (m[4] - 2 * I(g * G * F**2))).subs(sol)
        d6 = (u[6].eval(1) - (m[5] - r(1, 2) * I(g * F**3)
                               - r(9, 2) * I(g * G**2 * F**2) + I(f * G**2 * H))).subs(sol)
        print(f'n={n}:', [sp.expand(d) == 0 for d in (d4, d5, d6)])


def gine():
    import mpmath as mp
    mp.mp.dps = 30
    th = sp.Symbol('th')
    f = sp.sin(th) - sp.sin(2 * th) + sp.sin(3 * th)
    g = sp.cos(th) + 2 * sp.cos(2 * th)
    G = sp.integrate(g, (th, 0, th))
    F = sp.integrate(f, (th, 0, th))
    fl, gl, Gl, Fl = (sp.lambdify(th, e, 'mpmath') for e in (f, g, G, F))
    m = [mp.quad(lambda s: fl(s) * Gl(s)**k, [0, 2 * mp.pi]) for k in range(4)]
    print('m_0..m_3 =', [mp.nstr(x, 12) for x in m])
    print('(1/2) int g F^2 =', mp.nstr(mp.quad(lambda s: gl(s) * Fl(s)**2, [0, 2 * mp.pi]) / 2, 20),
          ' pi/2 =', mp.nstr(mp.pi / 2, 20))


def determinante():
    def M(n, N):
        Gx = (t**n - 1) / n
        return sp.Matrix(N + 1, N + 1, lambda k, j: sp.integrate(t**(2 * j) * Gx**k, (t, -1, 1)))

    def formula(n, N):
        num = 2**(N + 1) * sp.Rational(1, n)**(N * (N + 1) // 2)
        for i in range(N + 1):
            for k in range(i + 1, N + 1):
                num *= 2 * n * (k - i)**2
        den = 1
        for i in range(N + 1):
            for j in range(N + 1):
                den *= n * i + 2 * j + 1
        return num / den
    for n in (2, 4, 6, 8):
        print(f'n={n}:', [sp.simplify(M(n, N).det() - formula(n, N)) == 0 for N in range(1, 5)])


def hiperplano(n):
    e = sp.symbols('e0:4')
    g = P(t**(n - 1), sp.QQ)
    G = prim(g)
    fe = sum(e[j] * t**(2 * j) for j in range(4))
    ms = [sp.integrate(fe * G.as_expr()**k, (t, -1, 1)) for k in range(3)]
    sol = sp.solve(ms, e[:3], dict=True)[0]
    eta = sp.expand(fe.subs(sol).subs(e[3], 1))
    print('e_n(t) =', sp.factor(eta))
    a = sp.Symbol('a')
    b = sp.symbols('b0:4')
    dom = sp.QQ[(a,) + b]
    f = P(a * eta + sum(b[j] * t**(2 * j + 1) for j in range(4)), dom)
    u = ucoeffs(f, P(t**(n - 1), dom), 5)
    for k in range(1, 6):
        print(f'u_{k}(1) =', sp.factor(u[k].eval(1)))


def groebner(n, E, R, K):
    a = sp.symbols('a0:%d' % (E + 1))
    b = sp.symbols('b0:%d' % (R + 1))
    syms = a + b
    dom = sp.QQ[syms]
    f = P(sum(a[j] * t**(2 * j) for j in range(E + 1)) + sum(b[j] * t**(2 * j + 1) for j in range(R + 1)), dom)
    t0 = time.time()
    u = ucoeffs(f, P(t**(n - 1), dom), K)
    eqs = [sp.expand(u[k].eval(1)) for k in range(1, K + 1)]
    z = sp.Symbol('z')
    ok = all(list(sp.groebner(eqs + [a[j] * z - 1], z, *syms, order='grevlex')) == [1] for j in range(E + 1))
    print(f'n={n}, deg f_e<={2*E}, deg f_o<={2*R+1}, K={K}: '
          f'u_1(1)=...=u_K(1)=0 => f_e=0 ? {ok}   ({time.time()-t0:.1f}s)')


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'identidades'
    if cmd == 'identidades':
        identidades()
    elif cmd == 'gine':
        gine()
    elif cmd == 'determinante':
        determinante()
    elif cmd == 'hiperplano':
        hiperplano(int(sys.argv[2]))
    elif cmd == 'groebner':
        groebner(*map(int, sys.argv[2:6]))
