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
            thetav = H*nu/KB
            S += thetavT/(np.exp(thetavT)-1.0)-\
                 np.log(1-np.exp(-1.0*thetavT))
            E += T*thetavT/(np.exp(thetavT)-1.0)+thetav/2.0 #Evib +E_ZPE
        Svib.append(S)
        Evib.append(E)
    return(R*np.array(Svib),R*np.array(Evib))

def rotationalSE (*args,**kwargs):
    kw = {}
    kw.update(kwargs) 
    Trange = kw.get('Trange',None)
    T = kw.get ('T',298.15)
    sig = kw.get('sig',None)
    lin = kw.get('lin',None)
    thetar = kw.get ('thetar',None)
    thetaA = kw.get ('thetaA',None)
    thetaB = kw.get ('thetaB',None)
    thetaC = kw.get ('thetaC',None)
    assert sig is not None and lin is not None
    #M = M*1e-3/NA
    Trange = kw.get('Trange',None)
    T = kw.get('T',298.15)
    if Trange is None:
        Trange = np.array([T])
    else:
        Trange = np.array(Trange)
    Srot = []
    Erot = []
    for T in Trange:
        if lin:
            buf = T/(sig*thetar)
            Srot.append(np.log(buf)+1.0)
            Erot.append(T)
        else:
            buf1 = np.sqrt(np.pi)/sig
            buf2 = T**3/(thetaA*thetaB*thetaC)
            Srot.append(np.log(buf1*buf2**0.5)+1.5)
            Erot.append(1.5*T)
    return (R*np.array(Srot),R*np.array(Erot))


            
    

