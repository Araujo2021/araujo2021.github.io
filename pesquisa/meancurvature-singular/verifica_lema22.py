"""Verificação numérica do Lema 2.2 e de (25) para a função a(v) definida em (4).

Uso: python3 verifica_lema22.py
"""
import numpy as np

b = 1 / (8 * np.sqrt(2))
c = 7 * b


def a(v):
    return np.where(v < 1, 1 / np.sqrt(1 + v), np.where(v < 2, b * (v - 2) ** 2 + c, c))


def da(v):
    return np.where(v < 1, -0.5 / (1 + v) ** 1.5, np.where(v < 2, 2 * b * (v - 2), 0.0))


s = np.linspace(1e-6, 4, 2_000_001)
v = s**2
num = a(v) + 2 * v * da(v)          # (s a(s^2))' = autovalor radial de D[a(|xi|^2) xi]
ratio = num / a(v)                  # quociente em (11)

i = ratio.argmin()
print(f"(11) min de (s a(s^2))'/a(s^2) = {ratio[i]:.5f} em s^2 = {v[i]:.5f}  (artigo afirma >= 1/2)")
print(f"     valor exato em s^2 = (11-sqrt(77))/2 = {(11-np.sqrt(77))/2:.5f}")
print(f"(14) artigo afirma > 5/8 = 0.625 em [1, sqrt 2)")
j = num.argmin()
print(f"autovalor radial minimo nu = {num[j]:.5f} em s^2 = {v[j]:.5f}; 19/(40 sqrt 2) = {19/(40*np.sqrt(2)):.5f}")
print(f"cota usada no artigo p.5-6 para o numerador: (7 - sqrt 2)/(8 sqrt 2) = {(7-np.sqrt(2))/(8*np.sqrt(2)):.5f}")

# (25): t^2 phi(t) / Phi(t), phi(t) = a(t^2), Phi(t) = int_0^t s a(s^2) ds
t = np.linspace(0, 30, 3_000_001)
g = t * a(t**2)
Phi = np.concatenate([[0], np.cumsum((g[1:] + g[:-1]) / 2 * np.diff(t))])
r = t[10:] ** 2 * a(t[10:] ** 2) / Phi[10:]
print(f"(25) inf t^2 phi/Phi = {r.min():.5f} em t = {t[10 + r.argmin()]:.4f}; sup = {r.max():.5f}  (artigo: [3/2, 2])")
