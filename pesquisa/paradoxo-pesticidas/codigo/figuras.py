"""Gera as figuras do texto. Uso: python3 figuras.py [fig1 fig2 ...]"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pest import Model, psi_plateau, phi_bump, phi_rect, H_max_from_level

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figuras')
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})
C1, C2, C3, CG = '#1f5fa8', '#c0392b', '#2e8b57', '#888888'
a, b, c, d = 1., 1., .5, 1.
Hs, Ps = d / (b * c), a / b

def V(H, P):  # energia pontual (>=0)
    return c * (H - Hs - Hs * np.log(H / Hs)) + (P - Ps - Ps * np.log(P / Ps))

def ode_model():
    return Model(DH=0., DP=0., N=3)

def ode_run(H0, P0, tf, apps, phi=phi_bump, dt=0.01):
    m = ode_model()
    t, H, P = m.simulate(np.full(3, H0), np.full(3, P0), tf, apps, phi=phi,
                         dt_out=dt)
    return t, H[:, 0], P[:, 0]

# ------------------------------------------------------------------ Fig 1
def fig1():
    sH, sP = 1.0, 0.5
    t, H, P = ode_run(1.0, 0.5, 15, [])
    # primeiro periodo: tempos de H maximo e minimo
    i1 = np.argmax(H[t < 8]); Tgood = t[i1]
    j = np.argmin(H[(t > Tgood) & (t < Tgood + 8)]); Tbad = t[t > Tgood][j]
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    ax[0].plot(H, P, color=CG, lw=1, label='sem pesticida')
    for T, col, lab in [(Tgood, C1, f'T={Tgood:.2f} (H alto)'),
                        (Tbad, C2, f'T={Tbad:.2f} (H baixo)')]:
        tt, HH, PP = ode_run(1.0, 0.5, T + 25, [(T, 0.1, sH, sP)])
        m = tt >= T
        ax[0].plot(HH[m], PP[m], color=col, lw=1.2, label=lab)
        ax[0].plot(HH[np.searchsorted(tt, T)], PP[np.searchsorted(tt, T)], 'o', color=col, ms=4)
    ax[0].plot(Hs, Ps, 'k+')
    ax[0].set_xlabel('H (praga)'); ax[0].set_ylabel('P (predador)')
    ax[0].set_title('(a) plano de fase, $\\varepsilon=0.1$'); ax[0].legend(fontsize=8)
    # varredura do tempo de aplicacao (um periodo)
    m = (t >= Tgood) & (t <= Tgood + 8)
    Ts = np.linspace(0.2, 7.5, 120)
    peak, dV, Hat = [], [], []
    for T in Ts:
        tt, HH, PP = ode_run(1.0, 0.5, T + 20, [(T, 0., sH, sP)], dt=0.02)
        k = np.searchsorted(tt, T, side='right')
        peak.append(HH[k:].max())
        Hm, Pm = np.interp(T, t, H), np.interp(T, t, P)
        Hat.append(Hm)
        dV.append(V(Hm * np.exp(-sH), Pm * np.exp(-sP)) - V(Hm, Pm))
    peak0 = H.max()
    ax[1].plot(Ts, peak, color=C1, label='pico de H após aplicação')
    ax[1].plot(Ts, Hat, color=CG, ls='--', label='H(T) no momento da aplicação')
    ax[1].axhline(peak0, color='k', lw=.8, ls=':', label='pico sem pesticida')
    ax[1].set_xlabel('instante de aplicação T'); ax[1].set_title('(b) efeito do instante')
    ax[1].legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5,-0.18), ncol=2)
    ax[2].plot(Ts, dV, color=C2)
    ax[2].axhline(0, color='k', lw=.8)
    ax[2].fill_between(Ts, dV, 0, where=np.array(dV) > 0, color=C2, alpha=.15)
    ax[2].fill_between(Ts, dV, 0, where=np.array(dV) < 0, color=C1, alpha=.15)
    ax[2].set_xlabel('instante de aplicação T')
    ax[2].set_title('(c) salto $[\\![V]\\!]$ no limite impulsivo')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig1_edo_instante.pdf'); plt.close(fig)
    print('fig1: Tgood', Tgood, 'Tbad', Tbad)

# ------------------------------------------------------------------ Fig 2
def fig2():
    s = np.linspace(1e-3, 4, 200)
    thr = Hs * s / (1 - np.exp(-s))
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].plot(s, thr, color='k')
    ax[0].fill_between(s, thr, 9, color=C1, alpha=.15)
    ax[0].fill_between(s, 0, thr, color=C2, alpha=.15)
    ax[0].axhline(Hs, color=CG, ls='--', lw=.8)
    ax[0].text(0.3, 6.5, 'aplicação reduz $V$\n(controle)', color=C1)
    ax[0].text(2.2, 1.0, 'aplicação aumenta $V$\n(ressurgência)', color=C2)
    ax[0].text(3.2, Hs + .15, '$H^*$', color=CG)
    ax[0].set_xlim(0, 4); ax[0].set_ylim(0, 9)
    ax[0].set_xlabel('dose letal $s=\\alpha_H\\psi$'); ax[0].set_ylabel('H no instante da aplicação')
    ax[0].set_title('(a) limiar $H^* s/(1-e^{-s})$')
    # verificacao numerica com eps pequeno (praga apenas)
    rng = np.random.default_rng(1)
    pred, num = [], []
    for _ in range(40):
        H0, P0 = rng.uniform(.3, 6), rng.uniform(.3, 3)
        sH = rng.uniform(.1, 3)
        for eps in (0.01,):
            tt, HH, PP = ode_run(H0, P0, eps, [(0., eps, sH, 0.)], dt=eps / 20)
            num.append(V(HH[-1], PP[-1]) - V(H0, P0))
            pred.append(c * (Hs * sH - H0 * (1 - np.exp(-sH))))
    ax[1].plot(pred, num, 'o', ms=3, color=C1)
    lim = [min(pred + num), max(pred + num)]
    ax[1].plot(lim, lim, 'k', lw=.8)
    ax[1].set_xlabel('salto previsto (limite impulsivo)')
    ax[1].set_ylabel('variação numérica, $\\varepsilon=0.01$')
    ax[1].set_title('(b) verificação da fórmula do salto')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig2_limiar_dose.pdf'); plt.close(fig)
    print('fig2: max abs err', np.max(np.abs(np.array(pred) - np.array(num))))

# ------------------------------------------------------------------ Fig 3
def fig3():
    m = Model(N=201)
    x = m.x
    H0 = 1 + .8 * np.cos(2 * np.pi * x / m.L); P0 = .5 + .4 * np.cos(2 * np.pi * x / m.L)
    psi = psi_plateau(x, 5, 6, 1.5)
    T, tf = 10.5, 15.5
    epss = 0.4 / 2 ** np.arange(0, 7)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    res = {}
    for kill, aH, aP in [('prop', 1.0, .5), ('add', 1.0, .1)]:
        m.kill = kill
        _, Href, Pref = m.simulate(H0, P0, tf, [(T, 0., aH, aP)], psiH=psi, dt_out=1)
        for phi, nm in [(phi_bump, 'bump $C^\\infty_c$'), (phi_rect, 'retângulo (artigo)')]:
            errs = []
            for eps in epss:
                _, He, Pe = m.simulate(H0, P0, tf, [(T, eps, aH, aP)], psiH=psi,
                                       phi=phi, dt_out=1, rtol=1e-10, atol=1e-12)
                errs.append(max(abs(He[-1] - Href[-1]).max(), abs(Pe[-1] - Pref[-1]).max()))
            res[(kill, nm)] = errs
    for k, (key, errs) in enumerate(res.items()):
        kill, nm = key
        ax[0 if kill == 'prop' else 1].loglog(epss, errs, 'o-', label=nm, ms=4,
                                             color=[C1, C2][k % 2])
    for i, ttl in enumerate(['(a) mortalidade proporcional', '(b) remoção aditiva (artigo)']):
        e0 = res[list(res)[2 * i]][0]
        ax[i].loglog(epss, e0 * epss / epss[0], 'k:', lw=.8, label='$O(\\varepsilon)$')
        ax[i].set_xlabel('$\\varepsilon$'); ax[i].set_ylabel('$\\|u_\\varepsilon-u_0\\|_\\infty$ em $t=T+5$')
        ax[i].set_title(ttl); ax[i].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f'{OUT}/fig3_convergencia.pdf'); plt.close(fig)
    for k, v in res.items():
        r = np.log2(np.array(v[:-1]) / np.array(v[1:]))
        print('fig3', k, ['%.2e' % e for e in v], 'ordens', np.round(r, 2))

# ------------------------------------------------------------------ Fig 4
def fig4():
    m = Model(N=201)
    x = m.x
    H0 = 1 + .8 * np.cos(2 * np.pi * x / m.L); P0 = .5 + .4 * np.cos(2 * np.pi * x / m.L)
    psi = psi_plateau(x, 5, 8, 1.)
    aH, aP, eps, tf = 1.2, .6, .2, 45
    scen = [(None, 'sem pesticida', CG), (10.5, 'T=10.5 ($\\bar H$ alto)', C1),
            (6.0, 'T=6.0 ($\\bar H$ baixo)', C2)]
    out = {}
    for T, lab, col in scen:
        apps = [] if T is None else [(T, eps, aH, aP)]
        t, H, P = m.simulate(H0, P0, tf, apps, psiH=psi, dt_out=.05)
        E = np.array([m.energy(h, p) for h, p in zip(H, P)])
        out[lab] = (t, H, P, E, col, T)
    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 3)
    for k, lab in enumerate(list(out)[1:]):
        t, H, P, E, col, T = out[lab]
        axh = fig.add_subplot(gs[0, k])
        im = axh.pcolormesh(x, t, H, shading='auto', cmap='viridis', vmin=0, vmax=10, rasterized=True)
        axh.axhline(T, color='w', lw=.8, ls='--')
        axh.set_xlabel('x'); axh.set_ylabel('t'); axh.set_title(f'H(x,t), {lab}')
        fig.colorbar(im, ax=axh)
    axp = fig.add_subplot(gs[0, 2])
    axp.plot(x, psi, color='k'); axp.set_title('perfil espacial $\\psi(x)$ (suporte compacto)')
    axp.set_xlabel('x')
    axE = fig.add_subplot(gs[1, 0]); axM = fig.add_subplot(gs[1, 1]); axF = fig.add_subplot(gs[1, 2])
    for lab, (t, H, P, E, col, T) in out.items():
        axE.plot(t, E, color=col, label=lab)
        axM.plot(t, H.max(1), color=col, label=lab, lw=1)
        Hb = H @ m.w / m.L; Pb = P @ m.w / m.L
        T0 = T if T else 0
        mm = (t >= T0) & (t <= T0 + 20)
        axF.plot(Hb[mm], Pb[mm], color=col, lw=1, label=lab)
        if T:
            k = np.searchsorted(t, T + eps)
            hb = H_max_from_level(E[k] / m.L)
            axF.axvline(hb, color=col, ls=':', lw=.8)
    axE.set_title('(d) energia $\\mathcal{E}(t)$'); axE.set_xlabel('t'); axE.legend(fontsize=8)
    axM.axhline(6, color='k', ls='--', lw=.8); axM.text(1, 6.2, 'limiar econômico', fontsize=8)
    axM.set_title('(e) pico local $\\max_x H(x,t)$'); axM.set_xlabel('t')
    axF.set_title('(f) médias espaciais e cota de Jensen'); axF.set_xlabel('$\\bar H$'); axF.set_ylabel('$\\bar P$')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig4_espacial.pdf'); plt.close(fig)
    for lab, (t, H, P, E, col, T) in out.items():
        if T:
            k0 = np.searchsorted(t, T) - 1; k1 = np.searchsorted(t, T + eps)
            print('fig4', lab, 'E antes %.3f depois %.3f' % (E[k0], E[k1]),
                  'pico pos %.2f' % H[t > T].max())
        else:
            print('fig4', lab, 'pico (t>5) %.2f' % H[t > 5].max())

# ------------------------------------------------------------------ Fig 5
def fig5():
    m = Model(N=101)
    x = m.x
    H0 = 1 + .8 * np.cos(2 * np.pi * x / m.L); P0 = .5 + .4 * np.cos(2 * np.pi * x / m.L)
    aH, aP, eps, hor = 1.2, .6, .2, 25
    Ts = np.linspace(4, 12, 33)
    ws = np.array([2, 3, 4, 5, 6, 7, 8, 9])
    t0, H0r, P0r = m.simulate(H0, P0, 12 + hor, [], dt_out=.05)
    peak = np.zeros((len(ws), len(Ts))); dE = np.zeros_like(peak)
    for i, w in enumerate(ws):
        psi = psi_plateau(x, 5, w, .8)
        for j, T in enumerate(Ts):
            t, H, P = m.simulate(H0, P0, T + hor, [(T, eps, aH, aP)], psiH=psi, dt_out=.05)
            k = t > T
            k0 = np.searchsorted(t, T) - 1; k1 = np.searchsorted(t, T + eps)
            ref = H0r[(t0 > T) & (t0 <= T + hor)].max()
            peak[i, j] = H[k].max() / ref
            dE[i, j] = m.energy(H[k1], P[k1]) - m.energy(H[k0], P[k0])
    Hb = H0r @ m.w / m.L
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].plot(t0, Hb, color=C1, label='$\\bar H$'); ax[0].plot(t0, H0r[:, 50], color=C3, label='H(5,t) (centro)')
    ax[0].axhline(Hs, color=CG, ls='--', lw=.8); ax[0].set_xlim(4, 12)
    ax[0].set_title('(a) dinâmica sem pesticida'); ax[0].set_xlabel('t'); ax[0].legend(fontsize=8)
    im = ax[1].pcolormesh(Ts, ws, peak, cmap='RdBu_r', vmin=0, vmax=2, shading='auto', rasterized=True)
    ax[1].set_xlabel('instante de aplicação T'); ax[1].set_ylabel('largura do talhão pulverizado w')
    ax[1].set_title('(b) pico pós-aplicação / pico sem pesticida'); fig.colorbar(im, ax=ax[1])
    vm = np.abs(dE).max()
    im = ax[2].pcolormesh(Ts, ws, dE, cmap='RdBu_r', vmin=-vm, vmax=vm, shading='auto', rasterized=True)
    ax[2].contour(Ts, ws, dE, [0], colors='k', linewidths=.8)
    ax[2].set_xlabel('instante de aplicação T'); ax[2].set_title('(c) salto de energia $[\\![\\mathcal{E}]\\!]$')
    fig.colorbar(im, ax=ax[2])
    fig.tight_layout(); fig.savefig(f'{OUT}/fig5_talhao.pdf'); plt.close(fig)
    np.savez(f'{OUT}/fig5_dados.npz', Ts=Ts, ws=ws, peak=peak, dE=dE)
    print('fig5 peak range', peak.min(), peak.max())

# ------------------------------------------------------------------ Fig 6
def dose_otima(H, P, kappa, smax, lam=0.):
    """Minimiza pontualmente j(s)=c[H*s-H(1-e^-s)]+[P*k s-P(1-e^-ks)]+lam*s
    em s in [0,smax] (regra da Proposicao de dose otima)."""
    from scipy.optimize import brentq
    s = np.zeros_like(H)
    for i, (h, p) in enumerate(zip(H, P)):
        dj = lambda z: c * Hs - c * h * np.exp(-z) + kappa * (Ps - p * np.exp(-kappa * z)) + lam
        if dj(0) >= 0:
            s[i] = 0.
        elif dj(smax) <= 0:
            s[i] = smax
        else:
            s[i] = brentq(dj, 0, smax)
    return s

def fig6():
    m = Model(N=201)
    x = m.x
    H0 = 1 + .8 * np.cos(2 * np.pi * x / m.L); P0 = .5 + .4 * np.cos(2 * np.pi * x / m.L)
    kappa, smax, eps, hor = .5, 3., .2, 30
    unif = 1.2 * psi_plateau(x, 5, 8, 1.)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    t0, Hr, Pr = m.simulate(H0, P0, 6 + hor + 5, [], dt_out=.05)
    for col_ax, T in [(1, 10.5), (2, 6.0)]:
        k = np.searchsorted(t0, T)
        HT, PT = Hr[k], Pr[k]
        sopt = dose_otima(HT, PT, kappa, smax)
        M = 0.5 * m.integral(sopt)
        lo, hi = 0., 5.
        for _ in range(50 if M > 0 else 0):
            lam = (lo + hi) / 2
            if m.integral(dose_otima(HT, PT, kappa, smax, lam)) > M: lo = lam
            else: hi = lam
        sbud = dose_otima(HT, PT, kappa, smax, hi)
        strat = [('sem pesticida', None, CG), ('uniforme', unif, C2),
                 ('ótima (sem orçamento)', sopt, C1), ('ótima, metade da quantidade', sbud, C3)]
        if T == 10.5:
            ax[0].plot(x, HT, 'k', lw=1, label='H(T,x)')
            ax[0].plot(x, PT, 'k--', lw=1, label='P(T,x)')
            for lab, sx, colr in strat[1:]:
                ax[0].plot(x, sx, color=colr, label='dose ' + lab)
            ax[0].axhline(Hs, color=CG, lw=.6, ls=':')
            ax[0].set_title('(a) estado em T=10.5 e doses $s_H(x)$'); ax[0].set_xlabel('x')
            ax[0].legend(fontsize=7, loc='upper right', bbox_to_anchor=(1.0, 0.72))
        print(f'fig6 T={T}')
        for lab, sx, colr in strat:
            apps = [] if sx is None else [(T, eps, 1., 1.)]
            psiH = np.zeros_like(x) if sx is None else sx
            t, H, P = m.simulate(H0, P0, T + hor, apps, psiH=psiH, psiP=kappa * psiH, dt_out=.05)
            mm = t > T
            ax[col_ax].plot(t[t > T - 3], H[t > T - 3].max(1), color=colr, lw=1.1, label=lab)
            k0 = np.searchsorted(t, T) - 1; k1 = np.searchsorted(t, T + eps)
            dE = m.energy(H[k1], P[k1]) - m.energy(H[k0], P[k0])
            print(f'   {lab:28s} quantidade={m.integral(psiH):5.2f}  salto E={dE:+7.3f}  pico pos={H[mm].max():5.2f}')
        ax[col_ax].axvline(T, color='k', lw=.6, ls=':')
        ax[col_ax].set_xlabel('t'); ax[col_ax].set_ylabel('$\\max_x H$')
        ax[col_ax].set_title(f'({"bc"[col_ax-1]}) aplicação em T={T}')
        ax[col_ax].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{OUT}/fig6_dose_otima.pdf'); plt.close(fig)

if __name__ == '__main__':
    which = sys.argv[1:] or ['fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6']
    for f in which:
        globals()[f]()
