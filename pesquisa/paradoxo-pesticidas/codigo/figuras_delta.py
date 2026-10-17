"""Figuras da familia com largura de transicao delta. Uso: python3 figuras_delta.py [g1 g2]"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pest2d import Planar, c_star_lambda, k_suave, faixa_delta
from pest import phi_bump

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figuras')
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})
CS = ['#1f5fa8', '#2e8b57', '#d68910', '#8e44ad', '#c0392b']
L = 10.

def g1():
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4))
    # (a) curva c*(lambda) e colapso de simulacoes lineares na faixa
    lams = np.logspace(-2, 1.7, 30)
    for al, col in zip([0.6, 1.2, 2.4], CS):
        c0 = c_star_lambda(al, 0)
        cl = np.array([c_star_lambda(al, l) for l in lams])
        k = k_suave(al)[0]
        ax[0].loglog(lams, cl, color=col, label=f'$\\alpha={al}$')
        ax[0].loglog(lams[lams > 1], k / lams[lams > 1], ':', color=col, lw=.8)
        ax[0].axhline(c0, color=col, ls='--', lw=.6)
        print(f'g1a alpha={al}: c*(0)={c0:.4f} k={k:.4f}', flush=True)
    m = Planar(a=0., b=1e-30, c=.5, d=0., Nx=4000, Ny=1, Lx=L, Ly=1.)
    h0, al, D = 4., 1.2, m.DH
    H0 = np.full((m.Nx, 1), h0); P0 = np.ones((m.Nx, 1))
    pts = []
    for delta in [0.02, 0.1, 0.5]:
        z = faixa_delta(m.X, 3., 7., delta)
        for eps in [0.05, 0.2, 0.8]:
            _, He, _ = m.simulate(H0, P0, eps, [(0., eps, al, 0.)], chiH=z, store=False,
                                  rtol=1e-10, atol=1e-12)
            exc = m.integral(H0 * np.exp(-al * z) - He[-1])
            lam = delta / np.sqrt(D * eps)
            pts.append((lam, exc / (2 * h0 * np.sqrt(D * eps))))
    pts = np.array(pts)
    ax[0].loglog(pts[:, 0], pts[:, 1], 'o', mfc='none', color='k', ms=6,
                 label='simulações 2D (faixa), $\\alpha=1{,}2$')
    print('g1a pontos', pts.tolist(), flush=True)
    ax[0].set_xlabel('$\\lambda=\\delta/\\sqrt{D\\varepsilon}$')
    ax[0].set_ylabel('$c^*(\\alpha,\\varphi,\\lambda)$')
    ax[0].set_title('(a) constante de borda: nítida (- -) e suave $k/\\lambda$ (:)')
    ax[0].legend(fontsize=7, loc='lower left'); ax[0].set_ylim(1e-3, 0.5)
    # (b) modelo completo na faixa: erro L1 em t=T+eps para varios delta
    m = Planar(Nx=4000, Ny=1, Lx=L, Ly=1.)
    X = m.X
    H0 = 1 + .8 * np.cos(2 * np.pi * X / L); P0 = .5 + .4 * np.cos(2 * np.pi * X / L)
    T = 10.5
    _, HT, PT = m.simulate(H0, P0, T, store=False, rtol=1e-10, atol=1e-12)
    HT, PT = HT[-1], PT[-1]
    epss = 0.2 / 2 ** np.arange(0, 7)
    for delta, col in zip([0, 0.03, 0.1, 0.3, 1.0], CS):
        z = faixa_delta(X, 3., 7., delta)
        err = []
        for eps in epss:
            _, He, Pe = m.simulate(HT, PT, eps, [(0., eps, 1., .5)], chiH=z, store=False,
                                   rtol=1e-10, atol=1e-12)
            _, Hr, Pr = m.simulate(HT, PT, eps, [(0., 0., 1., .5)], chiH=z, store=False,
                                   rtol=1e-10, atol=1e-12)
            err.append(m.integral(abs(He[-1] - Hr[-1]) + abs(Pe[-1] - Pr[-1])))
        err = np.array(err)
        ordens = np.log2(err[:-1] / err[1:])
        print(f'g1b delta={delta}:', ['%.2e' % e for e in err], 'ordens', np.round(ordens, 2), flush=True)
        ax[1].loglog(epss, err, 'o-', ms=3, color=col, label=f'$\\delta={delta}$')
        if delta > 0:
            es = delta ** 2 / m.DH
            if epss[-1] < es < epss[0]:
                ax[1].plot(es, np.interp(np.log(es), np.log(epss[::-1]), err[::-1]), 'kx', ms=7)
    ax[1].loglog(epss, 0.9 * np.sqrt(epss / epss[0]), 'k:', lw=.8, label='$\\sqrt\\varepsilon$')
    ax[1].loglog(epss, 0.25 * epss / epss[0], 'k-.', lw=.8, label='$\\varepsilon$')
    ax[1].set_xlabel('$\\varepsilon$'); ax[1].set_ylabel('$\\|u_\\varepsilon-u_0\\|_{L^1}$ em $t=T+\\varepsilon$')
    ax[1].set_title('(b) erro $\\sim\\varepsilon/\\max(\\delta,\\sqrt{D\\varepsilon})$;  x: $\\varepsilon=\\delta^2/D_H$')
    ax[1].legend(fontsize=7, ncol=2)
    # (c) seletividade: razao predador/praga da mortalidade de borda (mesma dose e densidade)
    lam = np.logspace(-2, 1.6, 40)
    al = 1.2
    cache = {}
    def cs(l):
        key = round(float(l), 10)
        if key not in cache:
            cache[key] = c_star_lambda(al, l)
        return cache[key]
    for r, col in [(2., CS[4]), (4., CS[3]), (0.5, CS[0])]:
        R = [np.sqrt(r) * cs(l / np.sqrt(r)) / cs(l) for l in lam]
        ax[2].semilogx(lam, R, color=col, label=f'$D_P/D_H={r:g}$')
        ax[2].axhline(np.sqrt(r), color=col, ls='--', lw=.6)
        ax[2].axhline(r, color=col, ls=':', lw=.6)
        print(f'g1c r={r}: R(lam min)={R[0]:.3f} R(lam max)={R[-1]:.3f}', flush=True)
    ax[2].set_xlabel('$\\delta/\\sqrt{D_H\\varepsilon}$')
    ax[2].set_ylabel('mortalidade de borda predador / praga')
    ax[2].set_title('(c) seletividade: $\\sqrt{D_P/D_H}$ (- -) $\\to D_P/D_H$ (:)')
    ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f'{OUT}/p8_familia_delta.pdf'); plt.close(fig)

if __name__ == '__main__':
    for f in sys.argv[1:] or ['g1']:
        globals()[f]()
