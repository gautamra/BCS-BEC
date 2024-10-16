#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from BCS_BEC_lib import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import fft, ifft
import sys
import pandas as pd
import os , datetime

plt.style.use(['tableau-colorblind10'])
font = {'family' : 'sans serif',
        'weight' : 'regular',
        'size'   : 7}
mpl.rc('font', **font) 
mpl.rcParams['lines.linewidth'] = 0.66
mpl.rcParams['lines.markersize'] = 3
figsize=(3.375,4)
plt.rcParams["figure.figsize"] = figsize
plt.rcParams["figure.dpi"] = 200

def homogeneous_chain(LL, t, v, u0, PBC):
    if PBC:
        ts = [-t]*LL
    else:
        ts = [0]+([-t]*(LL-1))
    vs = [v]*LL
    u0s= [u0]*LL
    
    return ts, u0s, vs

def segmented_chain(LL_list, t_list, v_list, u0_list, t_int_list):
    assert len(LL_list)==len(t_list)==len(v_list)==len(u0_list)==len(t_int_list), "Lengths of supplied parameters not consistent"
    ts, vs, u0s = [], [], []
    for t, t_int, v, u0, LL in zip(t_list, t_int_list, v_list, u0_list, LL_list):
        ts += [-t_int]+[-t]*(LL-1)
        vs += [v]*LL
        u0s += [u0]*LL
    return ts, u0s, vs    
    
def symmetric_SNS_junction(LS, LN, V, t, U0S, U0N, delphi, edge_sites_with_fixed_phase=5):
    LL_list = [LS, LN, LS]
    t_list = [t, t, t]
    v_list = [V, 0, V]
    u0_list = [U0S, U0N, U0S]
    t_int_list = [0, t, t]
    
    ts, u0s, vs = segmented_chain(LL_list, t_list, v_list, u0_list, t_int_list)
    LL = np.sum(LL_list)
        
    def constrain_phase(Delta):
        Delta = Delta*np.exp(-1j*np.angle(Delta[0]))
        Delta[:edge_sites_with_fixed_phase]=np.abs(Delta[:edge_sites_with_fixed_phase])
        Delta[-edge_sites_with_fixed_phase:]=np.abs(Delta[-edge_sites_with_fixed_phase:])*np.exp(1j*delphi)
        return Delta
    
    return TBmodel(LL, ts, u0s, vs, constrain_phase=constrain_phase, Delta_init=Delta_init)

class SNS_junction:
    def __init__(self, LS, LN, V, t, U0S, U0N, edge_sites_with_fixed_phase=5):
        self.LS=LS
        self.LN=LN
        self.V=V
        self.t=t
        self.U0S=U0S
        self.U0N=U0N
        self.edge_sites_with_fixed_phase=edge_sites_with_fixed_phase
        self.LL_list = [LS, LN, LS]
        self.t_list = [t, t, t]
        self.v_list = [V, 0, V]
        self.u0_list = [U0S, U0N, U0S]
        self.t_int_list = [0, t, t]
        self.LL = np.sum(self.LL_list)

    def constrain_phase(self, Delta, delphi):
        Delta = Delta*np.exp(-1j*np.angle(Delta[0]))
        Delta[:self.edge_sites_with_fixed_phase]=np.abs(Delta[:self.edge_sites_with_fixed_phase])
        Delta[-self.edge_sites_with_fixed_phase:]=np.abs(Delta[-self.edge_sites_with_fixed_phase:])*np.exp(1j*delphi)
        return Delta

    def constrain_phase(self, Delta, delphi):
        phi_avg = np.average(np.angle(Delta))
        Delta[:self.edge_sites_with_fixed_phase]=np.abs(Delta[:self.edge_sites_with_fixed_phase])*np.exp(1j*(phi_avg-delphi/2))
        Delta[-self.edge_sites_with_fixed_phase:]=np.abs(Delta[-self.edge_sites_with_fixed_phase:])*np.exp(1j*(phi_avg+delphi/2))
        return Delta
    
    def get_Model_full(self, delphi):
        ts, u0s, vs = segmented_chain(self.LL_list, self.t_list, self.v_list, self.u0_list, self.t_int_list)
        Delta_init = np.array([np.exp(-1j*delphi/2)]*self.LS+[0.]*self.LN+[np.exp(1j*delphi/2)]*self.LS)
        return TBmodel(self.LL, ts, u0s, vs, constrain_phase=lambda Delta: self.constrain_phase(Delta, delphi), Delta_init=Delta_init)

    def get_Model_piece(self, piece):
        if piece=='S':
            ts, u0s, vs = homogeneous_chain(self.LS, self.t, self.V, self.U0S, True)
            return TBmodel(self.LS, ts, u0s, vs)
        elif piece=='N':
            ts, u0s, vs = homogeneous_chain(self.LN, self.t, 0, self.U0N, True)
            return TBmodel(self.LN, ts, u0s, vs)
        else:
            raise ValueError('piece can only by S or N in a symmetric SNS junction')
            return

            
def main():
    V = -2.0
    mus = sys.argv[1]
    mun = sys.argv[2]
    N = sys.argv[4]
    Junction = SNS_junction(LS=101, LN=int(N), V=V, t=1, U0S=-V/2.0 + np.float64(mus), U0N=np.float64(mun), edge_sites_with_fixed_phase=5)

    ModelS = Junction.get_Model_piece('S')
    ModelN = Junction.get_Model_piece('N')
    
    data_num = int(sys.argv[3])
    dphis = np.linspace(-3.0 * np.pi / 2, 3.0 * np.pi / 2, data_num)
    supcurrents = np.zeros(data_num)
    current1 = np.zeros(data_num)
    current2 = np.zeros(data_num)
    current3 = np.zeros(data_num)
    
    i = 0

    for dphi in dphis:
        Model = Junction.get_Model_full(delphi=dphi)
        resS = ModelS.iterate()
        resN = ModelN.iterate()
        Del, Pot, Evals, Uvecs, Vvecs, H, occu_raw = Model.iterate()
        current = np.zeros(Model.LL, dtype='complex')
        current_n = np.zeros(Model.LL, dtype='complex')
        
        
        n = 0
        for vvec in Model.vvecs.T:
            current_n = -2.0 * (np.array(Model.ts) * (np.roll(vvec, 1) - np.roll(vvec, -1)) * (vvec.conjugate())).imag
            if n == 0:
                current1[i] = current_n[int(Model.LL/2)]
            elif n == 1:
                current2[i] = current_n[int(Model.LL/2)]
            elif n == 3:
                current3[i] = current_n[int(Model.LL/2)]
            n += 1

            current += current_n
        
        source_current = np.zeros(Model.LL, dtype='complex')
        #for uvec, vvec in zip(Model.uvecs.T, Model.vvecs.T):
        #    source_current += 4.0 * np.cumsum((np.array(Model.Delta) * (uvec.conjugate() * vvec))).imag

        for uvec, vvec in zip(Uvecs.T, Vvecs.T):
            source_current += 4.0 * np.cumsum((np.array(Model.Delta) * (uvec.conjugate() * vvec))).imag
            
        N_andreev = 0
        for evals in Evals.T:
            if evals < np.max(np.abs(Del)):
                N_andreev += 1
    
        #supcurrents[i] = current[int(Model.LL / 2)] + source_current[int(Model.LL / 2)]
        supcurrents[i] = current[int(Model.LL / 2)]
        i = i + 1

    fig = plt.figure(figsize=(3.5, 4), dpi=200)
    gs = fig.add_gridspec(2, 3, wspace=0.1, hspace=0.6)

    axl = fig.add_subplot(gs[0, 0])
    axc = fig.add_subplot(gs[0, 1])
    axr = fig.add_subplot(gs[0, 2])
    axb = fig.add_subplot(gs[1, :])
    
    ymax=np.max([(axl.get_ylim()), (axc.get_ylim()), (axr.get_ylim())])
    ymin=np.min([(axl.get_ylim()), (axc.get_ylim()), (axr.get_ylim())])
    ylim=(ymin, ymax)
    #axl.set_ylim(ylim)
    axc.set_xlabel('DOS')
    axl.set_ylabel('$E$')
    axc.set_yticks([]), axr.set_yticks([])

    DOS_S, eex = ModelS.get_DOS()
    DOS_N, eex = ModelN.get_DOS()
    axl.plot(DOS_S, eex)
    axc.plot(DOS_N, eex)
    axr.plot(DOS_S, eex)
    axb.plot(dphis, supcurrents , '.')
    #axb.plot(dphis, current1 , 'r.')
    #axb.plot(dphis, current2 , 'm.')
    #axb.plot(dphis, current3 , 'k.')
    axb.set_ylabel('$J$')
    axb.set_xlabel('$\\Delta\\phi$')
    axb.set_title('$\\mu_S=$ ' + f'{mus}' + ' $\\mu_N$ ' + f'{mun} ' + ' $N=$' + f'{N} ' + ' Num of Andreev Lvls = ' + f'{N_andreev}')
    
    # Saving the output plots
    current_date = datetime.datetime.now().strftime('%Y%m%d')
    folder_name = f"Results_{current_date}"
    os.makedirs(folder_name, exist_ok=True)
    output_plot_path = os.path.join(folder_name, f'Current_mus={mus}_mun={mun}_L={N}.png')
    output_excel_path = os.path.join(folder_name, f'Current_mus={mus}_mun={mun}_L={N}.xlsx')
    plt.savefig(output_plot_path)

    df = pd.DataFrame({'dphis': dphis, 'currents': supcurrents})
    df.to_excel(output_excel_path, index=False)
	


if __name__ == "__main__":
    main()
