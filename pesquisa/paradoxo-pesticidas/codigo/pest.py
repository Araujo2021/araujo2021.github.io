"""Modelo de Lotka-Volterra com difusao e aplicacao de pesticida
com funcoes de suporte compacto no tempo e no espaco.

  H_t = D_H H_xx + H(a - bP) - alpha_H rho_eps(t-T) psi(x) g(H)
  P_t = D_P P_xx + P(bcH - d) - alpha_P rho_eps(t-T) psi(x) g(P)

em (0, L), condicao de Neumann homogenea.
g(u) = u   (mortalidade proporcional, 'prop')
g(u) = 1   (remocao aditiva, como em Li & Yang 2015, 'add')
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.sparse import diags, bmat

# ---------------------------------------------------------------- perfis
def _f(z):
    z = np.asarray(z, dtype=float)
    out = np.zeros_like(z)
    m = z > 0
    out[m] = np.exp(-1.0 / z[m])
    return out

_BUMP_Z = quad(lambda s: np.exp(-1.0 / (s * (1 - s))), 0, 1)[0]

def phi_bump(tau):
    """Bump C^infty com suporte [0,1] e integral 1."""
    tau = np.asarray(tau, dtype=float)
    out = np.zeros_like(tau)
    m = (tau > 0) & (tau < 1)
    out[m] = np.exp(-1.0 / (tau[m] * (1 - tau[m]))) / _BUMP_Z
    return out

def phi_rect(tau):
    """Retangulo do artigo original: 1 em [0,1]."""
    tau = np.asarray(tau, dtype=float)
    return ((tau >= 0) & (tau <= 1)).astype(float)

def smooth_step(z):
    """C^infty, 0 para z<=0, 1 para z>=1."""
    return _f(z) / (_f(z) + _f(1 - z))

def psi_plateau(x, x0, w, delta):
    """Indicadora suavizada (C^infty) do talhao [x0-w/2, x0+w/2]:
    vale 1 no interior, 0 fora, transicoes de largura delta."""
    a, b = x0 - w / 2, x0 + w / 2
    return smooth_step((x - a) / delta) * smooth_step((b - x) / delta)

# ---------------------------------------------------------------- modelo
class Model:
    def __init__(self, a=1., b=1., c=.5, d=1., DH=.1, DP=.05, L=10., N=201,
                 kill='prop'):
        self.a, self.b, self.c, self.d = a, b, c, d
        self.DH, self.DP, self.L, self.N = DH, DP, L, N
        self.x = np.linspace(0, L, N)
        self.dx = self.x[1] - self.x[0]
        self.kill = kill
        self.Hs, self.Ps = d / (b * c), a / b
        # Laplaciano com Neumann (pontos fantasmas)
        e = np.ones(N)
        Lap = diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]).tolil()
        Lap[0, 1] = 2.0
        Lap[N - 1, N - 2] = 2.0
        self.Lap = (Lap / self.dx ** 2).tocsr()
        # pesos de quadratura (trapezio)
        self.w = np.full(N, self.dx); self.w[[0, -1]] = self.dx / 2

    def integral(self, u):
        return self.w @ u

    def energy(self, H, P):
        """E(t) = int [ c(H - H* - H* ln(H/H*)) + (P - P* - P* ln(P/P*)) ] dx >= 0"""
        Hs, Ps, c = self.Hs, self.Ps, self.c
        eH = c * (H - Hs - Hs * np.log(H / Hs))
        eP = (P - Ps - Ps * np.log(P / Ps))
        return self.integral(eH + eP)

    def rhs(self, t, y, apps, psiH, psiP, phi):
        N = self.N
        H, P = y[:N], y[N:]
        dH = self.DH * (self.Lap @ H) + H * (self.a - self.b * P)
        dP = self.DP * (self.Lap @ P) + P * (self.b * self.c * H - self.d)
        for (T, eps, aH, aP) in apps:
            if eps > 0 and T <= t <= T + eps:
                r = phi((t - T) / eps) / eps
                gH = H if self.kill == 'prop' else 1.0
                gP = P if self.kill == 'prop' else 1.0
                dH -= aH * r * psiH * gH
                dP -= aP * r * psiP * gP
        return np.concatenate([dH, dP])

    def jump(self, H, P, aH, aP, psiH, psiP):
        """Limite impulsivo (eps -> 0)."""
        if self.kill == 'prop':
            return H * np.exp(-aH * psiH), P * np.exp(-aP * psiP)
        return H - aH * psiH, P - aP * psiP

    def simulate(self, H0, P0, tf, apps=(), psiH=None, psiP=None,
                 phi=phi_bump, dt_out=0.05, rtol=1e-8, atol=1e-10):
        """apps: lista de (T, eps, alphaH, alphaP). eps=0 -> impulsivo."""
        N = self.N
        psiH = np.ones(N) if psiH is None else psiH
        psiP = psiH if psiP is None else psiP
        # pontos de quebra
        brk = {0.0, float(tf)}
        for (T, eps, *_ ) in apps:
            brk.add(T); brk.add(T + eps)
        brk = sorted(b for b in brk if 0 <= b <= tf)
        y = np.concatenate([H0, P0]).astype(float)
        ts, Ys = [0.0], [y.copy()]
        jac = bmat([[self.Lap, diags(np.ones(N))], [diags(np.ones(N)), self.Lap]])
        jac_sp = (jac != 0).astype(int)
        for t0, t1 in zip(brk[:-1], brk[1:]):
            # saltos impulsivos em t0
            for (T, eps, aH, aP) in apps:
                if eps == 0 and abs(T - t0) < 1e-14:
                    Hn, Pn = self.jump(y[:N], y[N:], aH, aP, psiH, psiP)
                    y = np.concatenate([Hn, Pn])
                    ts.append(t0); Ys.append(y.copy())
            if t1 <= t0:
                continue
            npts = max(2, int(np.ceil((t1 - t0) / dt_out)) + 1)
            te = np.linspace(t0, t1, npts)
            in_pulse = any(eps > 0 and T <= t0 < T + eps for (T, eps, *_ ) in apps)
            ms = (t1 - t0) / 50 if in_pulse else np.inf
            sol = solve_ivp(self.rhs, (t0, t1), y, t_eval=te, method='BDF',
                            jac_sparsity=jac_sp, rtol=rtol, atol=atol,
                            max_step=ms, args=(apps, psiH, psiP, phi))
            if not sol.success:
                raise RuntimeError(sol.message)
            y = sol.y[:, -1].copy()
            ts.extend(sol.t[1:]); Ys.extend(sol.y[:, 1:].T)
        Y = np.array(Ys)
        return np.array(ts), Y[:, :N], Y[:, N:]


def H_max_from_level(e, c=.5, Hs=2.):
    """Maior raiz h>H* de c(h - H* - H* ln(h/H*)) = e (cota do ciclo)."""
    from scipy.optimize import brentq
    if e <= 0:
        return Hs
    f = lambda h: c * (h - Hs - Hs * np.log(h / Hs)) - e
    hi = Hs * 2
    while f(hi) < 0:
        hi *= 2
    return brentq(f, Hs, hi)
