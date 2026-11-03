#!/usr/bin/env python3
"""Gera as figuras ilustrativas do material (pasta figuras/).

As figuras que são saída direta de uma célula do texto (marcadas com "% figura:" no .tex)
são geradas por verifica_codigos.py; este script gera as figuras "explicativas".

Uso:  python3 scripts/gera_figuras.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle, Polygon  # noqa: E402

DEST = Path(__file__).resolve().parent.parent / "figuras"
DEST.mkdir(exist_ok=True)

AZUL, VERM, VERDE, LARANJA, ROXO, CINZA = "#1F4E9A", "#B22222", "#1B7F3B", "#D2691E", "#6A3D9A", "#777777"
plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.3, "savefig.bbox": "tight", "figure.dpi": 120,
})


def salva(fig, nome):
    fig.savefig(DEST / f"{nome}.pdf")
    plt.close(fig)
    print("figura:", nome)


def eixos_centrais(ax):
    ax.spines["left"].set_position("zero")
    ax.spines["bottom"].set_position("zero")


# ---------------------------------------------------------------- Aula 3
def aula03():
    # piso e teto
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for ax, func, nome, cor in [(axs[0], np.floor, r"piso $\lfloor x\rfloor$", AZUL),
                                (axs[1], np.ceil, r"teto $\lceil x\rceil$", VERM)]:
        for k in range(-3, 3):
            x = np.array([k, k + 1])
            y = func(k + 0.5)
            ax.plot(x, [y, y], color=cor, lw=2.2)
            fechado = k if func is np.floor else k + 1
            aberto = k + 1 if func is np.floor else k
            ax.plot(fechado, y, "o", color=cor, ms=6)
            ax.plot(aberto, y, "o", mfc="white", color=cor, ms=6)
        ax.set_title(nome)
        ax.set_xlim(-3.3, 3.3); ax.set_ylim(-3.5, 3.5)
        eixos_centrais(ax)
        ax.set_xticks(range(-3, 4)); ax.set_yticks(range(-3, 4))
    salva(fig, "aula03-piso-teto")

    # método babilônico
    x, y = 2.0, 1.0
    ys = [y]
    for _ in range(6):
        y = (y + x / y) / 2
        ys.append(y)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.9), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axs[0]
    cores = [CINZA, AZUL, VERDE, VERM]
    for k, cor in zip(range(4), cores):
        a = ys[k]
        ax.add_patch(Rectangle((0, 0), a, x / a, fill=False, lw=2, ec=cor,
                               label=fr"$y_{k}={a:.4f}$,  $2/y_{k}={x/a:.4f}$"))
    r = np.sqrt(2)
    ax.add_patch(Rectangle((0, 0), r, r, fill=True, fc="gold", alpha=0.25, ec="none", label=r"quadrado de lado $\sqrt{2}$"))
    ax.set_xlim(0, 2.2); ax.set_ylim(0, 2.2); ax.set_aspect("equal")
    ax.set_title("retângulos de área 2"); ax.legend(fontsize=7.5, loc="upper right")
    ax = axs[1]
    erros = [abs(v - np.sqrt(2)) for v in ys]
    k = np.arange(len(ys))
    erros_plot = [max(e, 1e-16) for e in erros]
    ax.semilogy(k, erros_plot, "o-", color=AZUL)
    for kk, e in zip(k, erros):
        ax.annotate(f"{e:.1e}", (kk, max(e, 1e-16)), textcoords="offset points", xytext=(4, 6), fontsize=7.5)
    ax.set_xlabel("iteração $k$"); ax.set_ylabel(r"erro $|y_k-\sqrt{2}|$")
    ax.set_title("convergência")
    salva(fig, "aula03-babilonico")

    # bissecção
    f = lambda t: t**3 - t - 2
    fig, ax = plt.subplots(figsize=(8, 3.6))
    t = np.linspace(0.9, 2.1, 300)
    ax.plot(t, f(t), color=AZUL, lw=2, label=r"$f(x)=x^3-x-2$")
    ax.axhline(0, color="black", lw=0.8)
    a, b = 1.0, 2.0
    for k in range(5):
        yk = -2.6 - 0.55 * k
        ax.plot([a, b], [yk, yk], "|-", color=VERM, lw=2, ms=10)
        ax.text(b + 0.03, yk, f"$[a_{k},b_{k}]$", va="center", fontsize=8, color=VERM)
        m = (a + b) / 2
        ax.plot(m, yk, "o", color=VERDE, ms=4)
        if f(a) * f(m) < 0:
            b = m
        else:
            a = m
    raiz = 1.5213797068
    ax.plot(raiz, 0, "o", color="black")
    ax.annotate("raiz", (raiz, 0), xytext=(raiz - 0.25, 2.5), arrowprops=dict(arrowstyle="->"))
    ax.set_ylim(-5.3, 7); ax.set_xlim(0.9, 2.25)
    ax.legend(loc="upper left")
    salva(fig, "aula03-bisseccao")


# ---------------------------------------------------------------- Aula 4
def aula04():
    f = lambda t: np.exp(-t**2)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.3), sharey=True)
    t = np.linspace(-2.2, 0.3, 300)
    for ax, n in zip(axs, [4, 8]):
        x = np.linspace(-2, 0, n + 1)
        for a, b in zip(x[:-1], x[1:]):
            ax.add_patch(Polygon([(a, 0), (a, f(a)), (b, f(b)), (b, 0)], closed=True,
                                 fc=LARANJA, alpha=0.3, ec=LARANJA))
        ax.plot(t, f(t), color=AZUL, lw=2, label=r"$e^{-x^2}$")
        ax.plot(x, f(x), "o", color=VERM, ms=4)
        h = 2 / n
        I = h * (f(x).sum() - (f(x[0]) + f(x[-1])) / 2)
        ax.set_title(f"$n={n}$:  soma dos trapézios = {I:.5f}")
        ax.set_xlim(-2.2, 0.3); ax.set_ylim(0, 1.1)
        ax.legend(loc="upper left")
    fig.suptitle(r"valor exato: $\int_{-2}^0 e^{-x^2}dx \approx 0.88208$", y=1.02, fontsize=10)
    salva(fig, "aula04-trapezios")

    rng = np.random.default_rng(1)
    x, y = rng.random(2000), rng.random(2000)
    dentro = x**2 + y**2 <= 1
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(x[dentro], y[dentro], ".", ms=2.5, color=AZUL)
    ax.plot(x[~dentro], y[~dentro], ".", ms=2.5, color=VERM)
    s = np.linspace(0, np.pi / 2, 100)
    ax.plot(np.cos(s), np.sin(s), color="black", lw=1.5)
    ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title(f"4 x (azuis)/(total) = {4 * dentro.mean():.3f}")
    salva(fig, "aula04-montecarlo")


# ---------------------------------------------------------------- Aula 5
def aula05():
    fig, ax = plt.subplots(figsize=(5, 4))
    t = np.linspace(-3, 2, 100)
    ax.plot(t, (1 - t) / 2, color=AZUL, lw=2, label=r"$x_0+2x_1=1$")
    ax.plot(t, (2 - 3 * t) / 5, color=VERDE, lw=2, label=r"$3x_0+5x_1=2$")
    ax.plot(-1, 1, "o", color=VERM, ms=8)
    ax.annotate("solução $(-1,\\,1)$", (-1, 1), xytext=(-0.6, 1.6), arrowprops=dict(arrowstyle="->"))
    ax.set_xlabel("$x_0$"); ax.set_ylabel("$x_1$"); ax.legend(loc="lower left")
    ax.set_xlim(-3, 2); ax.set_ylim(-1, 2.2)
    salva(fig, "aula05-retas")

    A = np.array([[2, 1], [1, 2]])
    s_ = np.linspace(0, 2 * np.pi, 200)
    C = np.vstack([np.cos(s_), np.sin(s_)])
    E = A @ C
    lam, V = np.linalg.eigh(A)
    fig, axs = plt.subplots(1, 2, figsize=(9, 4.2))
    for ax, P, tit, esc in [(axs[0], C, r"círculo unitário", [1, 1]), (axs[1], E, r"imagem por $A$", lam)]:
        ax.plot(P[0], P[1], color=AZUL, lw=2)
        for k in range(2):
            v = V[:, k] * esc[k]
            ax.annotate("", xy=v, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=VERM, lw=2))
            ax.annotate("", xy=-v, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=VERM, lw=2))
        ax.set_aspect("equal"); ax.set_xlim(-3.3, 3.3); ax.set_ylim(-3.3, 3.3)
        ax.set_title(tit); eixos_centrais(ax)
    axs[1].text(1.9, 2.5, r"$\lambda=3$", color=VERM)
    axs[1].text(-1.6, 1.0, r"$\lambda=1$", color=VERM)
    salva(fig, "aula05-transformacao")


# ---------------------------------------------------------------- Aula 9
def aula09():
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(1, 2, 1)
    t = np.linspace(0, 4, 200)
    ax.plot(t, np.sqrt(t), color=AZUL, lw=2, label=r"$y=\sqrt{x}$")
    ax.fill_between(t, 0, np.sqrt(t), color=AZUL, alpha=0.2)
    x0 = 2.5
    ax.plot([x0, x0], [0, np.sqrt(x0)], color=VERM, lw=2)
    ax.text(x0 + 0.08, 0.5, "raio $f(x)$", color=VERM)
    ax.set_xlim(0, 4.3); ax.set_ylim(0, 2.3); ax.legend(loc="upper left")
    ax.set_title("região")
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    u = np.linspace(0, 4, 40)
    v = np.linspace(0, 2 * np.pi, 50)
    U, Vv = np.meshgrid(u, v)
    R = np.sqrt(U)
    ax.plot_surface(U, R * np.cos(Vv), R * np.sin(Vv), cmap="Blues", alpha=0.75, linewidth=0)
    ax.plot(x0 * np.ones_like(v), np.sqrt(x0) * np.cos(v), np.sqrt(x0) * np.sin(v), color=VERM, lw=2)
    ax.set_title("sólido de revolução (disco em vermelho)")
    ax.set_xlabel("x"); ax.view_init(elev=20, azim=-60)
    salva(fig, "aula09-revolucao")


if __name__ == "__main__":
    aula03()
    aula04()
    aula05()
    aula09()
