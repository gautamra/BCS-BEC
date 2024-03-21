#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:25:27 2020

@author: gautam
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import kwant

from matplotlib import pyplot
#%matplotlib notebook
import tinyarray
import scipy.sparse.linalg as sla
import scipy.linalg as lin

from tqdm import tqdm
import pickle

from scipy import stats

import warnings

tau_x = tinyarray.array([[0, 1], [1, 0]])
tau_y = tinyarray.array([[0, -1j], [1j, 0]])
tau_z = tinyarray.array([[1, 0], [0, -1]])
tau_0 = tinyarray.array([[1, 0], [0, 1]])
tau_up = tinyarray.array([[0, 1], [0, 0]])
tau_dn = tinyarray.array([[0, 0], [1, 0]])


def Lorentzian(eex, ee, gam):
    return (gam/np.pi)*(1/((eex-ee)**2 + gam**2))

def Fermi(eps, beta = 'inf'):
    if beta == 'inf':
        #if (eps >= -pow(10,-10)) and (eps <= pow(10,-10)):
        #    return 0.5
        #else:        
        #    return int(eps<0)
        return int(eps<0)
    else:
        return 1/(1+np.exp(beta*eps))


class TBmodel:
    def __init__(self, LL, ts, us, vs, a = 1, beta = 'inf', constrain_phase=lambda Delta: Delta, Delta_init=None, Pot_init=None):
        # Not sure how to use thetas yet to enforce a phase difference between the two sides
        # Instead, I manually enforce it at the opposite ends of the two sides, allowing the rest of the system
        #     to compute self-consistently!
        
        self.LL = LL
        self.a = a
        self.ts, self.us, self.vs = ts, us, vs   # ts: hopping; us: on-site potential; vs: on-site attraction V
        
        if Delta_init is None:
            self.Delta = np.array([1]*self.LL, dtype = "complex")
        else:
            self.Delta = Delta_init
        if Pot_init is None:
            self.Pot = np.array(self.vs)/2
        else:
            self.Pot = Pot_init
        self.beta = beta
        self.make_syst()
        
        if constrain_phase is not None:
            self.constrain_phase=constrain_phase
        # else:
        #     def constrain_phase(Delta):
        #         return Delta
        #     self.constrain_phase=constrain_phase

    def onsite(self, site, Delta, Pot):
        (x,y) = site.tag
        return (self.us[x]+Pot[x])*tau_z - self.vs[x]*(Delta[x]*tau_up + np.conjugate(Delta[x])*tau_dn) 
   
    def hopping(self,site1,site2):
        (x2,y2) = site2.tag
        (x1,y1) = site1.tag
        # return self.ts[x1]*np.exp(1j*self.thetas[x1])*tau_z
        return self.ts[x1]*tau_z
   
    def make_syst(self):
        self.syst = kwant.Builder()
        self.lat = kwant.lattice.square(self.a, norbs = 2)
        self.syst[(self.lat(x,0) for x in range(self.LL))] = self.onsite
        self.syst[((self.lat(x+1,0),self.lat(x,0)) for x in range(self.LL-1))] = self.hopping
        self.syst[((self.lat(0,0), self.lat(self.LL-1,0)))] = self.hopping 
       
        self.fsyst = self.syst.finalized()
        return
    
    def solve(self, H , cc):
        (evals, evecs) = lin.eigh(H)
        uvecs = evecs[::2]
        vvecs = evecs[1::2]
        return (evals[self.LL:],uvecs[:,self.LL:],vvecs[:,self.LL:])
    
    def iterate(self):
        
        def self_cons(H , cc):
            (evals, uvecs, vvecs) = self.solve(H , cc)
            self.evals, self.uvecs, self.vvecs = (evals, uvecs, vvecs)
           
            Delta = np.zeros(self.LL, dtype = "complex")
            for ee, uvec, vvec in zip(evals, uvecs.T, vvecs.T):
                Delta += (1-2*Fermi(ee, beta = self.beta))*(uvec*vvec.conjugate()).reshape(self.LL)
            Delta=self.constrain_phase(Delta)
            
            occupancy = np.zeros(self.LL)
            for ee, uvec, vvec in zip(evals, uvecs.T, vvecs.T):
                    occupancy += (Fermi(ee, beta = self.beta)*np.abs(uvec)**2 + (1-Fermi(ee, beta = self.beta))*np.abs(vvec)**2).reshape(self.LL)     
            self.occupancy = occupancy
            
            Pot = self.vs*occupancy
        
            # Current = np.zeros(self.LL)
            # SourceCurr = np.zeros(self.LL , dtype = 'complex')
            # m = 1.0;
            # a = 1.0;
            # for ee, uvec, vvec in zip(evals, uvecs.T, vvecs.T):
            #         au = np.angle(uvec); au = np.concatenate((au , [au[0]]))
            #         av = np.angle(vvec); av = np.concatenate((av , [av[0]]))
            #         #Current += (2.0/m)*(Fermi(ee, beta = self.beta)*np.abs(uvec)**2 * (au[1:] - au[:-1])+ (1-Fermi(ee, beta = self.beta))*np.abs(vvec)**2 * (av[1:] - av[:-1])).reshape(self.LL)
            #         vvec = np.concatenate((vvec , [0]))
            #         Current += (2.0/m)*np.imag(vvec[:-1]*np.conjugate(vvec[1:])).reshape(self.LL)
            #         #SourceCurr -= 4*a*(Delta*np.conjugate(uvec)*vvec).reshape(self.LL)
                    
            Pot = Pot + 0.00001*np.ones(len(Pot))  ## Adding perturbation for self-consistency convergence
            return (Delta, Pot)

               
        err_Delta = np.ones(self.LL)
        cc = 0
       # definitions for testing the free energy calcualtion
        self.testDeltaF = []
        oldDeltaF = np.array([np.abs(self.Delta.mean()), 0])
        H = self.fsyst.hamiltonian_submatrix(params = dict(Delta = self.Delta, Pot = self.Pot))
        H_ini= H
        while np.array([(abs(Del)>10**(-5)) and ((abs(err)/abs(Del))>0.001) and cc < 200 for err,Del in zip(err_Delta, self.Delta)] ).any():
            
            H = self.fsyst.hamiltonian_submatrix(params = dict(Delta = self.Delta, Pot = self.Pot))
            newDelta, newPot = self_cons(H , cc) 
            newDelta = newDelta*3/4 + self.Delta*1/4
            newPot = newPot*3/4 + self.Pot*1/4
            err_Delta = np.abs(newDelta - self.Delta)
           
            free_energy = self.get_free_energy()
            DeltaF = np.array([abs(self.Delta.mean()), free_energy])
            self.testDeltaF.append(DeltaF)
           
            cc += 1    
            self.Delta, self.Pot = newDelta, newPot
       
        if cc<200:
            print("Convergence took {} iterations".format(cc))  
        else:
            print("Convergence not reached")
            
        self.Delta, self.Pot = self_cons(H , cc)
        
        self.H = H
        
        return self.Delta, self.Pot, self.evals, self.uvecs, self.vvecs, self.H, self.occupancy
       
    def get_DOS(self, gam = None, Num_es = 1000):
        # need to make these limits more general. Also in get_LDOS()
        emax = np.max(np.abs(self.evals))
        emin = -emax

       
        if gam == None:
            gam = 2*emax/self.LL
           
        eex = np.linspace(emin - (emax - emin)/10,emax + (emax - emin)/10, Num_es)
        DOSu = np.zeros(eex.shape)
        DOSv = np.zeros(eex.shape)
       
        for ee, uvec, vvec in zip(self.evals, self.uvecs.T, self.vvecs.T):
            if ee>0:
                DOSu += np.linalg.norm(uvec)**2*Lorentzian(eex,ee,gam)
                DOSv += np.linalg.norm(vvec)**2*Lorentzian(eex,-ee,gam)
               
        self.DOS = (DOSu + DOSv)/self.LL
        return  self.DOS , eex
   
    # Need to rewrite get_LDOS for 2D SSH model.
    def get_LDOS(self, gam = None, Num_es = 1000):
        emax = np.max(np.abs(self.evals))
        emin = -emax

        if gam == None:
            gam = 2*emax/self.LL
           
        eex = np.linspace(emin - (emax - emin)/10,emax + (emax - emin)/10, Num_es)
        DOSu = np.zeros((self.uvecs.shape[0],eex.shape[0]))
        DOSv = np.zeros(DOSu.shape)
       
        for ee, uvec, vvec in zip(self.evals, self.uvecs.T, self.vvecs.T):
            if ee>0:
                DOSu += (np.abs(uvec)**2)[:,np.newaxis]*Lorentzian(eex,ee,gam)
                DOSv += (np.abs(vvec)**2)[:,np.newaxis]*Lorentzian(eex,-ee,gam)      
           
        self.LDOS = (DOSu + DOSv)/self.LL
        uvecshape = (np.abs(uvec)**2).shape
        uvecplusshape = (np.abs(uvec)**2)[:,np.newaxis]
        return  self.LDOS,eex
   
    def get_free_energy(self):
        Energy_g = 0
        for ee, vvec in zip(self.evals, self.vvecs.T):
            Energy_g += -2*ee*np.linalg.norm(vvec)**2
        Energy_g2 = np.linalg.norm(np.abs(self.vs[0])*self.Delta)**2/np.abs([0])

        Energy_g = Energy_g + Energy_g2
       
        Energy_exc = 0
        for ee in self.evals:
            Energy_exc += 2*ee*Fermi(ee, beta = self.beta)
       
       
        Energy_entropy = 0
        for ee in self.evals:
            if self.beta != 'inf':
                term1 = Fermi(ee, beta = self.beta)*np.log(Fermi(ee, beta = self.beta))
                term2 = (1-Fermi(ee, beta = self.beta))*np.log((1-Fermi(ee, beta = self.beta)))
                Energy_entropy += -2/self.beta * (term1 + term2)
                           
        return Energy_g + Energy_exc + Energy_entropy
       
    def get_ham(self,inds):
        if inds == 'full':
            return self.fsyst.hamiltonian_submatrix(params = dict(Delta = self.Delta, Pot = self.Pot))
        else:
            return self.fsyst.hamiltonian(*inds, params = dict(Delta = self.Delta, Pot = self.Pot))
        
    # def get_current(self, ):
    #     # generalize for finite temperature
    #     current = np.zeros(self.LL, dtype='complex')
    #     for vvec in self.vvecs:
    #         current+=2*np.imag(vvec*np.roll(vvec.conjugate(),-1))
       
       
class DOS_SC:
    def __init__(self, eex, DOS):
        self.eex = eex
        self.DOS = DOS
       
        self.zeroindex = np.argmin(abs(self.eex))
       
   
    def find_coherence_peak(self):
        i = self.zeroindex
        while self.DOS[i-1]>self.DOS[i] or self.DOS[i+1]>self.DOS[i]:
            i += 1
        self.rightPeakIndex = i
       
        i = self.zeroindex
        while self.DOS[i-1]>self.DOS[i] or self.DOS[i+1]>self.DOS[i]:
            i += -1
        self.leftPeakIndex = i
       
        self.leftPeakHeight = self.DOS[self.leftPeakIndex]
        self.rightPeakHeight = self.DOS[self.rightPeakIndex]
        self.gapWidth = self.eex[self.rightPeakIndex] - self.eex[self.leftPeakIndex]
