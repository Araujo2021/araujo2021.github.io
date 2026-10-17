"""Modelo planar de Lotka-Volterra com difusao e pesticida aplicado
apenas numa zona omega do retangulo Omega = (0,Lx) x (0,Ly):

  H_t = D_H Lap H + H(a - bP) - alpha_H rho_eps(t-T) chi_omega(x) H
  P_t = D_P Lap P + P(bcH - d) - alpha_P rho_eps(t-T) chi_omega(x) P

Neumann homogeneo em todo o bordo. Volumes finitos centrados na celula;
chi_omega e aproximada pela fracao de area de cada celula contida em omega.
Com Ny = 1 o problema reduz-se exatamente ao caso de dados invariantes em y.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, identity, kron, bmat
from pest import phi_bump, phi_rect  # perfis temporais de suporte compacto


def lap1d(n, h):
    if n == 1:
        return diags([np.zeros(1)], [0])
    e = np.ones(n)
    L = diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]).tolil()
    L[0, 0] = -1.0          # Neumann (reflexao) em celulas centradas
    L[n - 1, n - 1] = -1.0
    return (L / h ** 2).tocsr()


class Planar:
    def __init__(self, a=1., b=1., c=.5, d=1., DH=.1, DP=.05, Lx=10., Ly=10.,
                 Nx=64, Ny=64):
        self.a, self.b, self.c, self.d = a, b, c, d
        self.DH, self.DP = DH, DP
        self.Lx, self.Ly, self.Nx, self.Ny = Lx, Ly, Nx, Ny
        self.hx, self.hy = Lx / Nx, Ly / Ny
        self.x = (np.arange(Nx) + .5) * self.hx
        self.y = (np.arange(Ny) + .5) * self.hy
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')
        self.Lap = (kron(lap1d(Nx, self.hx), identity(Ny)) +
                    kron(identity(Nx), lap1d(Ny, self.hy))).tocsr()
        self.cell = self.hx * self.hy
        self.area = Lx * Ly
        self.Hs, self.Ps = d / (b * c), a / b
        self.n = Nx * Ny

    # ------------------------------------------------------------ zonas
    def chi(self, inside, sub=8):
        """Fracao de area de cada celula dentro de omega = {inside(x,y)}."""
        s = (np.arange(sub) + .5) / sub - .5
        acc = np.zeros((self.Nx, self.Ny))
        for sx in s:
            for sy in s:
                acc += inside(self.X + sx * self.hx, self.Y + sy * self.hy)
        return acc / sub ** 2

    def integral(self, u):
        return u.sum() * self.cell

    def energy(self, H, P):
        Hs, Ps, c = self.Hs, self.Ps, self.c
        return self.integral(c * (H - Hs - Hs * np.log(H / Hs)) +
                             (P - Ps - Ps * np.log(P / Ps)))

    # ------------------------------------------------------------ dinamica
    def rhs(self, t, y, apps, chiH, chiP, phi):
        n = self.n
        H, P = y[:n], y[n:]
        dH = self.DH * (self.Lap @ H) + H * (self.a - self.b * P)
        dP = self.DP * (self.Lap @ P) + P * (self.b * self.c * H - self.d)
        for (T, eps, aH, aP) in apps:
            if eps > 0 and T <= t <= T + eps:
                r = phi((t - T) / eps) / eps
                dH -= aH * r * chiH * H
                dP -= aP * r * chiP * P
        return np.concatenate([dH, dP])

    def simulate(self, H0, P0, tf, apps=(), chiH=None, chiP=None, phi=phi_bump,
                 t_out=None, dt_out=0.1, rtol=1e-7, atol=1e-9, store=True):
        """apps: lista de (T, eps, alphaH, alphaP); eps = 0 => salto impulsivo.
        Retorna (t, H, P) com H, P de forma (nt, Nx, Ny)."""
        n = self.n
        chiH = np.ones(n) if chiH is None else np.ravel(chiH)
        chiP = chiH if chiP is None else np.ravel(chiP)
        brk = {0.0, float(tf)}
        for (T, eps, *_) in apps:
            brk.add(T); brk.add(T + eps)
        brk = sorted(b for b in brk if 0 <= b <= tf)
        if t_out is not None:
            t_out = np.asarray(t_out, float)
        y = np.concatenate([np.ravel(H0), np.ravel(P0)]).astype(float)
        ts, Ys = [0.0], [y.copy()]
        I = identity(n)
        jac_sp = (bmat([[self.Lap + I, I], [I, self.Lap + I]]) != 0).astype(int)
        for t0, t1 in zip(brk[:-1], brk[1:]):
            for (T, eps, aH, aP) in apps:
                if eps == 0 and abs(T - t0) < 1e-14:
                    y = np.concatenate([y[:n] * np.exp(-aH * chiH),
                                        y[n:] * np.exp(-aP * chiP)])
                    if store:
                        ts.append(t0); Ys.append(y.copy())
            if t1 <= t0:
                continue
            if t_out is None:
                npts = max(2, int(np.ceil((t1 - t0) / dt_out)) + 1)
                te = np.linspace(t0, t1, npts)
            else:
                te = np.unique(np.concatenate([[t0, t1], t_out[(t_out > t0) & (t_out < t1)]]))
            in_pulse = any(eps > 0 and T <= t0 < T + eps for (T, eps, *_) in apps)
            ms = (t1 - t0) / 50 if in_pulse else np.inf
            sol = solve_ivp(self.rhs, (t0, t1), y, t_eval=te if store else [t1],
                            method='BDF', jac_sparsity=jac_sp, rtol=rtol, atol=atol,
                            max_step=ms, args=(apps, chiH, chiP, phi))
            if not sol.success:
                raise RuntimeError(sol.message)
            y = sol.y[:, -1].copy()
            if store:
                ts.extend(sol.t[1:]); Ys.extend(sol.y[:, 1:].T)
        if not store:
            ts, Ys = [tf], [y]
        Y = np.array(Ys)
        sh = (len(ts), self.Nx, self.Ny)
        return np.array(ts), Y[:, :n].reshape(sh), Y[:, n:].reshape(sh)


def c_star(alpha, phi=phi_bump, Y=12., N=4801):
    """Constante do deficit de massa na camada limite (problema reescalado):
    u_tau = u_yy - alpha phi(tau) 1_{y>0} u em R, u(0)=1, tau in (0,1);
    c* = int [e^{-alpha 1_{y>0}} - u(1)] dy > 0."""
    from scipy.sparse import diags as _d
    y = np.linspace(-Y, Y, N); h = y[1] - y[0]
    e = np.ones(N)
    L = _d([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]).tolil()
    L[0, 1] = 2; L[-1, -2] = 2
    L = (L / h ** 2).tocsr()
    w = np.full(N, h); w[[0, -1]] = h / 2
    chi = (y > 0) + 0.5 * (y == 0)
    f = lambda t, u: L @ u - alpha * phi(np.array([t]))[0] * chi * u
    s = solve_ivp(f, (0, 1), np.ones(N), method='BDF', rtol=1e-10, atol=1e-12,
                  max_step=0.005)
    return w @ np.exp(-alpha * chi) - w @ s.y[:, -1]
