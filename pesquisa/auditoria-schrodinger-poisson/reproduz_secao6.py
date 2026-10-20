"""Reprodução independente da Seção 6 (ilustração numérica planar).

Esquema descrito no artigo: diferenças finitas de 5 pontos em (0,1)^2,
n = 100, h = 1/101, lambda = 1e-2, gamma = 1/2, relaxação omega = 0.5,
u0 = 0.1, phi0 = 0, parada quando o resíduo não escalado <= 1e-9.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

n, lam, gam, om, tol = 100, 1e-2, 0.5, 0.5, 1e-9
h = 1.0 / (n + 1)
I = sp.identity(n, format="csr")
T = sp.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1])
L = ((sp.kron(I, T) + sp.kron(T, I)) / h**2).tocsc()  # -Delta_h
L_lu = spla.splu(L)

u = np.full(n * n, 0.1)
phi = np.zeros(n * n)


def residuals(u, phi):
    Ru = L @ u + phi * u * np.exp(u**2) - lam * u**(-gam)
    Rp = L @ phi - u**2 * np.exp(u**2)
    return np.abs(Ru).max(), np.abs(Rp).max()


hist = []
for k in range(1, 1001):
    A = (L + sp.diags(phi * np.exp(u**2))).tocsc()
    uhat = spla.spsolve(A, lam * u**(-gam))
    u = (1 - om) * u + om * uhat
    phihat = L_lu.solve(u**2 * np.exp(u**2))
    phi = (1 - om) * phi + om * phihat
    ru, rp = residuals(u, phi)
    hist.append(max(ru, rp))
    if k <= 3 or max(ru, rp) <= tol:
        print(f"it {k:3d}  R_u={ru:.3e}  R_phi={rp:.3e}  R={max(ru, rp):.3e}")
    if max(ru, rp) <= tol:
        break

print("iterações:", k)
print("max u   =", f"{u.max():.6e}")
print("max phi =", f"{phi.max():.6e}")
print("e^{u^2}-1 max =", f"{np.expm1(u.max()**2):.3e}")
print("min u   =", f"{u.min():.3e}", " min phi =", f"{phi.min():.3e}")
