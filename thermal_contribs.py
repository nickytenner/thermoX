from .constants import *
import numpy as np

def translationalSE(*args,**kwargs):
    kw = {}
    kw.update(kwargs)
    P = kw.get('P',None)
    M = kw.get('M',None) 
    assert P is not None and M is not None and P > 0.0
    M = M/NA
    Trange = kw.get('Trange',None)
    T = kw.get('T',298.15)
    if Trange is None:
        Trange = np.array([T])
    else:
        Trange = np.array(Trange)
    Strans = []
    Etrans = []
    for T in Trange:
        k1 = KB*T/P
        k2 = 2*np.pi*M*KB*T/H**2
        Strans.append(R*np.log(k1*k2**1.5)+1.5)
        Etrans.append(1.5*R*T)
    return (np.array(Strans),np.array(Etrans))

def vibrationalSE (*args,**kwargs):
    kw = {}
    kw.update(kwargs)
    Trange = kw.get('Trange',None)
    T = kw.get('T',298.15)
    freq = np.array(kw.get('freq',\
           [139.670,149.331,209.610,401.178,815.054,3910.239]))*C
    if Trange is None:
        Trange = np.array([T])
    else:
        Trange = np.array(Trange)
    Svib = []
    Evib = []
    for T in Trange:
        S = 0.0
        E = 0.0
        for nu in freq:
            thetavT = H*nu/(KB*T)
            S += thetavT/(np.exp(thetavT)-1.0)-\
                 np.log(1-np.exp(-1.0*thetavT))
            E += T*thetavT/(np.exp(thetavT)-1.0)
        Svib.append(S)
        Evib.append(E)
    return(R*np.array(Svib),R*np.array(Evib))

def rotaionalSE (*args,**kwargs):
    kw = {}
    kw.update(kwargs)
    P = kw.get('P',None)
    M = kw.get('M',None) 
    assert P is not None and M is not None and P > 0.0
    #M = M*1e-3/NA
    Trange = kw.get('Trange',None)
    T = kw.get('T',298.15)
    if Trange is None:
        Trange = np.array([T])
    else:
        Trange = np.array(Trange)
    Strans = []
    Etrans = []


            
    

