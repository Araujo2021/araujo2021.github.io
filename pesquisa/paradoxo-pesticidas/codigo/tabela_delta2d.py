"""Tabela 2 do texto planar: efeito da largura de transicao delta no campo 2D
(predador mais movel, aplicacao longa eps=1). Uso: python3 tabela_delta2d.py"""
import numpy as np, sys
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from pest2d import Planar, Psi
m=Planar(Nx=80,Ny=80,DH=.05,DP=.2); X,Y=m.X,m.Y; L=10
H0=1+.8*np.cos(2*np.pi*X/L)*np.cos(np.pi*Y/L); P0=.5+.4*np.cos(2*np.pi*X/L)
T,eps,aH,aP,hor=9.0,1.0,1.2,1.2,25
_,HT,PT=m.simulate(H0,P0,T,store=False); HT,PT=HT[-1],PT[-1]
tn,Hn,Pn=m.simulate(HT,PT,hor,dt_out=.1)
print('sem: pico local %.2f  pico media %.2f'%(Hn.max(),Hn.mean(axis=(1,2)).max()))
def box(a,b,c,d,dl):
    if dl==0: return ((X>a)&(X<b)&(Y>c)&(Y<d)).astype(float)
    return Psi((X-a)/dl)*Psi((b-X)/dl)*Psi((Y-c)/dl)*Psi((d-Y)/dl)
def faixa(dl):  # refugio x<2.5
    return (X>2.5).astype(float) if dl==0 else Psi((X-2.5)/dl)
def refugios16(dl):
    r=np.zeros_like(X)
    for i in range(4):
        for j in range(4):
            cx,cy=1.25+2.5*i,1.25+2.5*j
            r+=box(cx-.625,cx+.625,cy-.625,cy+.625,dl)
    return 1-r
for nome,fz in [('faixa',faixa),('16 refugios',refugios16)]:
    for dl in [0,0.5,1,2]:
        z=fz(dl)
        t,H,P=m.simulate(HT,PT,hor,[(0.,eps,aH,aP)],chiH=z,dt_out=.1)
        k1=np.argmin(abs(t-eps))
        Hi,Pi=HT*np.exp(-aH*z),PT*np.exp(-aP*z)
        _,Hr,Pr=m.simulate(Hi,Pi,eps,store=False)
        E0=m.energy(HT,PT); E1=m.energy(H[k1],P[k1])
        print(f'{nome:12s} delta={dl:3}: area tratada={m.integral(z):5.1f} E {E0:.1f}->{E1:.1f}  morte extra (vs impulsivo) H={m.integral(Hr[-1]-H[k1]):+.3f} P={m.integral(Pr[-1]-P[k1]):+.3f}  pico local {H[t>eps].max():.2f}  pico media {H[t>eps].mean(axis=(1,2)).max():.2f}', flush=True)
