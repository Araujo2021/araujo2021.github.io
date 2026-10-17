"""Figuras do modelo planar com zona de refugio (chi_omega).
Uso: python3 figuras2d.py [f3 f4 f5 f6 f7]"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pest2d import Planar, c_star
from pest import phi_bump, phi_rect, psi_plateau, H_max_from_level

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figuras')
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})
C1, C2, C3, C4, CG = '#1f5fa8', '#c0392b', '#2e8b57', '#8e44ad', '#888888'
L = 10.

def init2d(m):
    X, Y = m.X, m.Y
    return (1 + .8 * np.cos(2 * np.pi * X / L) * np.cos(np.pi * Y / L),
            .5 + .4 * np.cos(2 * np.pi * X / L))

def init_strip(m):
    X = m.X
    return 1 + .8 * np.cos(2 * np.pi * X / L), .5 + .4 * np.cos(2 * np.pi * X / L)

def at(t, arr, tq):
    """valor em t=tq (ultimo indice com t<=tq+tol, i.e. apos eventual salto)"""
    k = np.where(np.abs(t - tq) < 1e-9)[0]
    return arr[k[-1]]

# ---------------------------------------------------------------- f3
def f3():
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    # (a) problema linear (sem reacao) na faixa: deficit de massa
    m = Planar(a=0., b=1e-30, c=.5, d=0., Nx=4000, Ny=1, Lx=L, Ly=1.)
    chi = m.chi(lambda X, Y: (X > 3) & (X < 7))
    h0, T, D = 4., 1., m.DH
    H0 = np.full((m.Nx, 1), h0); P0 = np.ones((m.Nx, 1))
    cs = c_star(1.0)
    epss = 0.1 / 4 ** np.arange(0, 5)
    tr, Hr, _ = m.simulate(H0, P0, T + 5, [(T, 0., 1., 0.)], chiH=chi, t_out=[T + 5],
                           rtol=1e-11, atol=1e-13)
    dm = []
    for eps in epss:
        te, He, _ = m.simulate(H0, P0, T + 5, [(T, eps, 1., 0.)], chiH=chi, t_out=[T + 5],
                               rtol=1e-11, atol=1e-13)
        dm.append(m.integral(Hr[-1] - He[-1]))
    pred = 2 * cs * np.sqrt(D * epss) * h0
    ax[0].loglog(epss, dm, 'o', color=C1, label='numérico (faixa, sem reação)')
    ax[0].loglog(epss, pred, 'k-', lw=.8, label='$c^*\\sqrt{D\\varepsilon}\\int_{\\partial\\omega}H\\,d\\sigma$')
    ax[0].set_xlabel('$\\varepsilon$'); ax[0].set_ylabel('massa morta a mais')
    ax[0].set_title('(a) mortalidade de borda $\\propto\\sqrt{\\varepsilon}$'); ax[0].legend(fontsize=8)
    print('f3a', list(zip(epss, dm, pred)))
    # (b) modelo completo na faixa: erros em t = T + eps
    m = Planar(Nx=4000, Ny=1, Lx=L, Ly=1.)
    H0, P0 = init_strip(m)
    T = 10.5
    zonas = {'$\\chi_\\omega$': m.chi(lambda X, Y: (X > 3) & (X < 7)),
             '$\\psi$ suave': psi_plateau(m.X, 5, 4, 1.)}
    epss = 0.2 / 2 ** np.arange(0, 6)
    res = {}
    for nm, z in zonas.items():
        e1, ei = [], []
        for eps in epss:
            to = [T + eps]
            te, He, Pe = m.simulate(H0, P0, T + eps, [(T, eps, 1., .5)], chiH=z, t_out=to,
                                    rtol=1e-10, atol=1e-12)
            tr, Hr, Pr = m.simulate(H0, P0, T + eps, [(T, 0., 1., .5)], chiH=z, t_out=to,
                                    rtol=1e-10, atol=1e-12)
            dH, dP = He[-1] - Hr[-1], Pe[-1] - Pr[-1]
            e1.append(m.integral(abs(dH) + abs(dP))); ei.append(max(abs(dH).max(), abs(dP).max()))
        res[nm] = (e1, ei)
        print('f3b', nm, ['%.2e' % v for v in e1], ['%.2e' % v for v in ei])
    for (nm, (e1, ei)), col in zip(res.items(), [C2, C1]):
        ax[1].loglog(epss, e1, 'o-', color=col, ms=4, label=nm + ', $L^1$')
        ax[1].loglog(epss, ei, 's--', color=col, ms=4, label=nm + ', $L^\\infty$')
    e = res['$\\chi_\\omega$'][0]
    ax[1].loglog(epss, e[-1] * np.sqrt(epss / epss[-1]), 'k:', lw=.8, label='$\\sqrt{\\varepsilon}$')
    e = res['$\\psi$ suave'][1]
    ax[1].loglog(epss, e[-1] * epss / epss[-1], 'k-.', lw=.8, label='$\\varepsilon$')
    ax[1].set_xlabel('$\\varepsilon$'); ax[1].set_title('(b) erro em $t=T+\\varepsilon$')
    ax[1].legend(fontsize=7, ncol=2)
    # (c) c*(alpha, phi)
    al = np.linspace(.1, 3, 15)
    for phi, nm, col in [(phi_bump, 'bump $C^\\infty_c$', C1), (phi_rect, 'retângulo', C2)]:
        cv = [c_star(a, phi) for a in al]
        ax[2].plot(al, cv, 'o-', ms=3, color=col, label=nm)
        print('f3c', nm, ['%.4f' % v for v in cv])
    ax[2].set_xlabel('dose $\\alpha$'); ax[2].set_ylabel('$c^*(\\alpha,\\varphi)$')
    ax[2].set_title('(c) constante de borda'); ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f'{OUT}/p3_convergencia_refugio.pdf'); plt.close(fig)

# ---------------------------------------------------------------- f4
def f4():
    m = Planar(Nx=64, Ny=64)
    H0, P0 = init2d(m)
    xr = 2.5
    chi = m.chi(lambda X, Y: X > xr)
    aH, aP, eps, tf = 1.2, .6, .2, 40
    cases = [(9.5, 'T=9.5 ($\\bar H$ alto)', C1), (6.5, 'T=6.5 ($\\bar H$ baixo)', C2)]
    snaps = lambda T: [T, T + eps, T + 2.5, T + 5.5]
    out = {}
    tn, Hn, Pn = m.simulate(H0, P0, tf, dt_out=.1)
    out['sem'] = (tn, Hn, Pn, None)
    for T, lab, col in cases:
        to = np.unique(np.concatenate([np.arange(0, tf + 1e-9, .1), snaps(T)]))
        t, H, P = m.simulate(H0, P0, tf, [(T, eps, aH, aP)], chiH=chi, t_out=to)
        out[lab] = (t, H, P, T)
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(3, 4)
    for r, (T, lab, col) in enumerate(cases):
        t, H, P, _ = out[lab]
        for k, ts in enumerate(snaps(T)):
            a = fig.add_subplot(gs[r, k])
            i = np.argmin(abs(t - ts)) if k else np.where(t <= T)[0][-1]
            im = a.imshow(H[i].T, origin='lower', extent=[0, L, 0, L], vmin=0, vmax=8,
                          cmap='viridis', interpolation='nearest')
            a.axvline(xr, color='w', ls='--', lw=.8)
            a.set_title(f'{lab.split(" ")[0]}: H em t={t[i]:.1f}', fontsize=9)
            a.set_xticks([0, 5, 10]); a.set_yticks([0, 5, 10])
            if k == 3:
                fig.colorbar(im, ax=a, fraction=.046)
    aE = fig.add_subplot(gs[2, 0:2]); aM = fig.add_subplot(gs[2, 2:4])
    for key, col, lab in [('sem', CG, 'sem pesticida')] + [(c[1], c[2], c[1]) for c in cases]:
        t, H, P, T = out[key]
        E = np.array([m.energy(h, p) for h, p in zip(H, P)])
        aE.plot(t, E, color=col, label=lab)
        aM.plot(t, H.max(axis=(1, 2)), color=col, lw=1, label=lab)
        aM.plot(t, H.mean(axis=(1, 2)), color=col, lw=1, ls=':')
        if T is not None:
            k1 = np.argmin(abs(t - (T + eps))); k0 = np.where(t < T)[0][-1]
            hb = H_max_from_level(E[k1] / m.area)
            aM.axhline(hb, color=col, ls='--', lw=.7)
            print('f4', lab, 'E: %.2f -> %.2f' % (E[k0], E[k1]), 'pico pos %.2f' % H[t > T].max(),
                  'media max pos %.2f' % H[t > T].mean(axis=(1, 2)).max(), 'cota %.2f' % hb)
        else:
            print('f4 sem: pico t>6 %.2f' % H[t > 6].max())
    aE.set_title('energia $\\mathcal{E}(t)$'); aE.set_xlabel('t'); aE.legend(fontsize=8)
    aM.set_title('$\\max H$ (cheio), $\\bar H$ (pontilhado), cota de Jensen (tracejado)')
    aM.set_xlabel('t')
    fig.tight_layout(); fig.savefig(f'{OUT}/p4_planar_cenario.pdf', dpi=150); plt.close(fig)

# ---------------------------------------------------------------- f5
def f5():
    m = Planar(Nx=40, Ny=40)
    H0, P0 = init2d(m)
    aH, aP, eps, hor = 1.2, .6, .2, 25
    Ts = np.linspace(4, 12, 17)
    xrs = np.arange(0, 9)          # largura da faixa de refugio
    tn, Hn, Pn = m.simulate(H0, P0, 12 + hor, dt_out=.1)
    peak = np.zeros((len(xrs), len(Ts))); dE = np.zeros_like(peak); pred = np.zeros_like(peak)
    for i, xr in enumerate(xrs):
        chi = m.chi(lambda X, Y: X > xr)
        for j, T in enumerate(Ts):
            to = np.unique(np.concatenate([np.arange(0, T + hor, .1), [T, T + eps]]))
            t, H, P = m.simulate(H0, P0, T + hor, [(T, eps, aH, aP)], chiH=chi, t_out=to)
            k0 = np.where(t <= T)[0][-1]; k1 = np.argmin(abs(t - (T + eps)))
            ref = Hn[(tn > T) & (tn <= T + hor)].max()
            peak[i, j] = H[t > T].max() / ref
            dE[i, j] = m.energy(H[k1], P[k1]) - m.energy(H[k0], P[k0])
            # salto previsto (impulsivo) via medias na zona
            Hm, Pm = H[k0], P[k0]
            jx = m.c * (m.Hs * aH - (1 - np.exp(-aH)) * Hm) + (m.Ps * aP - (1 - np.exp(-aP)) * Pm)
            pred[i, j] = m.integral(chi * jx)
        print('f5 xr', xr, 'feito')
    np.savez(f'{OUT}/p5_dados.npz', Ts=Ts, xrs=xrs, peak=peak, dE=dE, pred=pred)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].plot(tn, Hn.mean(axis=(1, 2)), color=C1, label='$\\bar H$')
    ax[0].plot(tn, Hn.max(axis=(1, 2)), color=C3, label='$\\max H$')
    ax[0].axhline(2, color=CG, ls='--', lw=.8); ax[0].set_xlim(4, 12)
    ax[0].set_xlabel('t'); ax[0].set_title('(a) dinâmica sem pesticida'); ax[0].legend(fontsize=8)
    frac = 1 - xrs / L
    im = ax[1].pcolormesh(Ts, frac, peak, cmap='RdBu_r', vmin=0, vmax=2, shading='auto', rasterized=True)
    ax[1].set_xlabel('instante de aplicação T'); ax[1].set_ylabel('fração pulverizada $|\\omega|/|\\Omega|$')
    ax[1].set_title('(b) pico pós-aplicação / pico sem pesticida'); fig.colorbar(im, ax=ax[1])
    vm = np.abs(dE).max()
    im = ax[2].pcolormesh(Ts, frac, dE, cmap='RdBu_r', vmin=-vm, vmax=vm, shading='auto', rasterized=True)
    ax[2].contour(Ts, frac, dE, [0], colors='k', linewidths=.8)
    ax[2].contour(Ts, frac, pred, [0], colors='y', linewidths=.8, linestyles='--')
    ax[2].set_xlabel('instante de aplicação T'); ax[2].set_title('(c) $[\\![\\mathcal{E}]\\!]$; tracejado: zero previsto')
    fig.colorbar(im, ax=ax[2])
    fig.tight_layout(); fig.savefig(f'{OUT}/p5_refugio_instante.pdf'); plt.close(fig)
    print('f5 peak', peak.min(), peak.max(), 'corr dE/pred', np.corrcoef(dE.ravel(), pred.ravel())[0, 1])

# ---------------------------------------------------------------- f6
def zona_otima(m, H, P, aH, aP, area):
    """Principio da banheira: omega* = {j < -t}, |omega*| = area (ou {j<0})."""
    j = m.c * (m.Hs * aH - (1 - np.exp(-aH)) * H) + (m.Ps * aP - (1 - np.exp(-aP)) * P)
    order = np.argsort(j, axis=None)
    ncell = int(round(area / m.cell))
    sel = np.zeros(j.size, bool)
    neg = np.sort(j, axis=None) < 0
    k = min(ncell, neg.sum())
    sel[order[:k]] = True
    return sel.reshape(j.shape).astype(float), j

def f6():
    m = Planar(Nx=64, Ny=64)
    H0, P0 = init2d(m)
    aH, aP, eps, hor = 1.2, .6, .2, 25
    xr = 2.5
    chi_faixa = m.chi(lambda X, Y: X > xr)
    area = m.integral(chi_faixa)
    tn, Hn, Pn = m.simulate(H0, P0, 12 + hor, t_out=[6.5, 9.5])
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.6))
    for r, T in enumerate([9.5, 6.5]):
        k = np.argmin(abs(tn - T))
        HT, PT = Hn[k], Pn[k]
        zopt, j = zona_otima(m, HT, PT, aH, aP, area)
        zhalf, _ = zona_otima(m, HT, PT, aH, aP, area / 3)
        a = ax[r, 0]
        vm = np.abs(j).max()
        im = a.imshow(j.T, origin='lower', extent=[0, L, 0, L], cmap='RdBu_r', vmin=-vm, vmax=vm)
        a.contour(m.x, m.y, zopt.T, [.5], colors='k', linewidths=1)
        a.axvline(xr, color=C3, ls='--', lw=1)
        a.set_title(f'T={T}: densidade de salto $j(x)$;\npreto: $\\omega^*$, verde: borda da faixa', fontsize=9)
        fig.colorbar(im, ax=a, fraction=.046)
        print(f'f6 T={T}: |omega*|={m.integral(zopt):.1f} (orcamento {area:.1f}), |{{j<0}}|={m.integral(j<0):.1f}')
        strat = [('sem pesticida', None, CG), ('faixa $x>2.5$', chi_faixa, C2),
                 ('$\\omega^*$ (mesma área máx.)', zopt, C1), ('$\\omega^*$ com 1/3 da área', zhalf, C3)]
        for lab, z, col in strat:
            to = np.unique(np.concatenate([np.arange(T - 2, T + hor, .1), [T, T + eps]]))
            to = to[to >= 0]
            apps = [] if z is None else [(T, eps, aH, aP)]
            t, H, P = m.simulate(H0, P0, T + hor, apps, chiH=z, t_out=to)
            mm = t >= T - 2
            ax[r, 1].plot(t[mm], H[mm].max(axis=(1, 2)), color=col, lw=1, label=lab)
            ax[r, 2].plot(t[mm], H[mm].mean(axis=(1, 2)), color=col, lw=1, label=lab)
            k0 = np.where(t <= T)[0][-1]; k1 = np.argmin(abs(t - (T + eps)))
            dE = m.energy(H[k1], P[k1]) - m.energy(H[k0], P[k0])
            ar = 0 if z is None else m.integral(z)
            print(f'   {lab:32s} area={ar:5.1f} salto E={dE:+8.2f} pico={H[t > T].max():.2f} media max={H[t > T].mean(axis=(1, 2)).max():.2f}')
        ax[r, 1].set_title(f'T={T}: $\\max H$'); ax[r, 2].set_title(f'T={T}: $\\bar H$')
        for c in (1, 2):
            ax[r, c].axvline(T, color='k', ls=':', lw=.6); ax[r, c].set_xlabel('t')
        ax[r, 1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{OUT}/p6_zona_otima.pdf'); plt.close(fig)

# ---------------------------------------------------------------- f7
def f7():
    """Mortalidade de borda no plano: problema linear (so difusao + abate),
    densidades uniformes h0, p0, mesma area tratada (36) com perimetros diferentes."""
    m = Planar(a=0., b=1e-30, c=.5, d=0., Nx=160, Ny=160)
    h0, p0, al = 4., 1., 1.2
    H0 = np.full((m.Nx, m.Ny), h0); P0 = np.full((m.Nx, m.Ny), p0)
    cs = c_star(al)
    def blocos(k):
        s, g = 6. / k, 10. / k
        return lambda X, Y: ((np.abs(np.mod(X, g) - g / 2) < s / 2) &
                             (np.abs(np.mod(Y, g) - g / 2) < s / 2))
    r = np.sqrt(36 / np.pi)
    zonas = [(f'{k*k} quadrado' + ('s' if k > 1 else ''), blocos(k), 24. * k) for k in (1, 2, 3, 4)]
    zonas.insert(1, ('disco', lambda X, Y: (X - 5) ** 2 + (Y - 5) ** 2 < r ** 2, 2 * np.pi * r))
    epss = [0.2, 0.8]
    res = []
    for nm, f, per in zonas:
        z = m.chi(f)
        for eps in epss:
            _, He, Pe = m.simulate(H0, P0, eps, [(0., eps, al, al)], chiH=z, store=False,
                                   rtol=1e-9, atol=1e-11)
            Hi, Pi = H0 * np.exp(-al * z), P0 * np.exp(-al * z)
            eH = m.integral(Hi - He[-1]); eP = m.integral(Pi - Pe[-1])
            pH = cs * np.sqrt(m.DH * eps) * h0 * per; pP = cs * np.sqrt(m.DP * eps) * p0 * per
            res.append((nm, per, eps, eH, pH, eP, pP))
            print(f'f7 {nm:12s} per={per:5.1f} area={m.integral(z):5.2f} eps={eps}: '
                  f'H {eH:.4f} (prev {pH:.4f})  P {eP:.4f} (prev {pP:.4f})', flush=True)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    for eps, mk in zip(epss, ['o', 's']):
        R = [x for x in res if x[2] == eps]
        per = [x[1] for x in R]
        ax[0].plot(per, [x[3] for x in R], mk, color=C1, label=f'praga, $\\varepsilon={eps}$')
        ax[0].plot(per, [x[5] for x in R], mk, color=C2, mfc='none', label=f'predador, $\\varepsilon={eps}$')
        pp = np.linspace(15, 100, 2)
        ax[0].plot(pp, cs * np.sqrt(m.DH * eps) * h0 * pp, '-', color=C1, lw=.7)
        ax[0].plot(pp, cs * np.sqrt(m.DP * eps) * p0 * pp, '--', color=C2, lw=.7)
    ax[0].set_xlabel('perímetro de $\\omega$ (área fixa $=36$)'); ax[0].set_ylabel('massa morta a mais')
    ax[0].set_title('(a) pontos: simulação; linhas: previsão'); ax[0].legend(fontsize=7)
    allp = [x[4] for x in res] + [x[6] for x in res]; alln = [x[3] for x in res] + [x[5] for x in res]
    ax[1].plot(allp, alln, 'o', ms=4, color=C3)
    lim = [0, max(allp + alln) * 1.05]
    ax[1].plot(lim, lim, 'k', lw=.7)
    for x in res:
        if x[0] == 'disco':
            ax[1].annotate('disco', (x[4], x[3]), fontsize=7, xytext=(4, -8), textcoords='offset points')
    ax[1].set_xlabel('previsão $c^*\\sqrt{D\\varepsilon}\\int_{\\partial\\omega}u\\,d\\sigma$')
    ax[1].set_ylabel('simulação 2D'); ax[1].set_title('(b) previsão vs. simulação')
    fig.tight_layout(); fig.savefig(f'{OUT}/p7_perimetro.pdf'); plt.close(fig)

if __name__ == '__main__':
    for f in sys.argv[1:] or ['f3', 'f4', 'f5', 'f6', 'f7']:
        globals()[f]()
