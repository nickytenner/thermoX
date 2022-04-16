import sys
import numpy as np
from .thermal_contribs import vibrationalSE, translationalSE, rotationalSE
from .constants import *
import periodictable as pt
import re

dist = lambda a,b:np.linalg.norm(a-b)

class Compound(object):
    def __init__(self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        self.fname = kw.get('file',None)
        self.sigma = None
        self.linear = False
        assert self.fname is not None
        with open(self.fname) as fp:
            lines = fp.readlines()
            buf = ""
            vbuf =""
        for line in lines:
            if search_strings[0] in line:
                buf = ""
            elif search_strings[1] in line:
                coord = buf
            else:
                buf=buf+line
            if search_strings[2] in line:
                vbuf =""
            elif search_strings[3] in line:
                vib = vbuf
            else:
                vbuf = vbuf+line
            if search_strings[4] in line:
                self.SPE = float(line.split()[-1])*627.50960803059
            if search_strings[5] in line:
                self.sigma = int(line.split()[-1])
            if search_strings[6] in line:
                self.linear = True
            if search_strings[7] in line:
                self.G298 = float(line.split()[-2])*627.50960803059
        assert self.sigma is not None
        coord = coord.split('\n')[1:-3]
        self.xyz = {}
        for i,axyz in enumerate(coord):
            a,x,y,z=axyz.split()
            self.xyz.update({a+str(i+1):np.array([float(x),float(y),float(z)])})
        vib = vib.split('\n')[4:-4]
        self.freq ={}
        for i,nu in enumerate(vib):
            freq = float(nu.split()[1])
            if freq == 0.0:
                pass
            elif freq > 0.0:
                self.freq.update({i:freq})
            else:
                print("ERROR: Found negative frequency ({}) in {}\nTERMINATING".\
                      format(freq,self.fname))
                sys.exit()
        self.getMolarMass()
        self.getMomOfInertia()
        self.getThetav()
        self.getSEvib()
        self.getSEtrans()
        self.getSErot()
    def getThetav(self):
        self.thetav=H*C*np.array(list(self.freq.values()))/KB
    def getMolarMass(self):
        self.M =0.0
        self.atomicmass = {}
        self.com = np.zeros(3)
        for key,val in self.xyz.items():
            elem = "".join(re.split("[^a-zA-Z]*", key))
            mass = getattr(pt,elem).mass*1E-3  #atomic mass in kg/mol
            self.M += mass
            self.atomicmass.update({key:mass})
            self.com += np.array(val)*mass
        self.com = self.com/self.M
    def getMomOfInertia(self):
        self.I=0.0; self.Ia=0.0; self.Ib=0.0; self.Ic=0.0
        self.thetar=0.0; self.thetaA=0.0;self.thetaB=0.0;self.thetaC=0.0
        for key,val in self.xyz.items():
            mass = self.atomicmass[key]/NA #mass per atom in Kg
            if self.linear:
                self.I += mass*dist(np.array(val),self.com)**2
            else:
                self.Ia += mass*np.sum([(val[i]-self.com[i])**2 for i in [1,2]])
                self.Ib += mass*np.sum([(val[i]-self.com[i])**2 for i in [0,2]])
                self.Ic += mass*np.sum([(val[i]-self.com[i])**2 for i in [0,1]])
        theta = H**2*1E20/(8*np.pi**2*KB) #1E20 is coming from the Anstroem^-2
        if self.linear:
            self.thetar = theta/self.I
        else:
            self.thetaA = theta/self.Ia
            self.thetaB = theta/self.Ib
            self.thetaC = theta/self.Ic
    def getSEtrans (self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        P = kw.get('P',101325)
        T = kw.get('T',298.15)
        Trange = kw.get('Trange',None)
        self.strans,self.etrans = translationalSE(P=P,Trange=Trange,M=self.M)
    def getSEvib (self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        T = kw.get('T',298.15)
        Trange = kw.get('Trange',None)
        self.svib,self.evib= vibrationalSE(Trange=Trange,freq=list(self.freq.values()))
    def getSErot (self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        T = kw.get('T',298.15)
        Trange = kw.get ('Trange',None)
        self.srot,self.erot = rotationalSE (Trange=Trange,sig=self.sigma,\
                             lin=self.linear,thetar=self.thetar,thetaA=self.thetaA,\
                             thetaB=self.thetaB,thetaC=self.thetaC)
    def getG (self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        T = kw.get('T',298.15)
        Trange = kw.get('Trange',None)
        x = kw.get ('x',1.0)
        P = kw.get ('P',101325)
        Px = P*x #Partial pressure
        self.getSEtrans(P=Px,T=T,Trange=Trange)
        self.getSEvib(T=T,Trange=Trange)
        self.getSErot(T=T,Trange=Trange)
        self.U = self.etrans+self.erot+self.evib+self.SPE
        self.S = self.strans+self.srot+self.svib
        if Trange is None:
            Trange = np.array(T)
        else:
            Trange = np.array(Trange)
        self.PV = R*Trange
        self.G = x*(self.U + self.PV - Trange*self.S)
        #return (self.G)
        
            
          
