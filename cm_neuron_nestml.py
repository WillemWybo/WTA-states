
from collections import namedtuple

import nest
import numpy as np
import matplotlib.pyplot as plt
import random
import statistics as stat
import sys
import yaml
from dataclasses import dataclass


#sleep = True

@dataclass
class ReceptorIdxs:
    GABA_soma: int = 0
    AMPA_soma: int = 1
    GABA_dist: int = 2
    AMPA_dist: int = 3
    REFRACT_soma: int = 4
    ADAPT_soma: int = 5
    BETA_dist: int = 6
    NMDA_soma: int = 7
    AMPA_NMDA_soma: int = 8
    ALPHAexc_soma: int = 9
    ALPHAinh_soma: int = 10
    ALPHAexc_dist: int = 11
    ALPHAinh_dist: int = 12
    NMDA_dist: int = 13
    AMPA_NMDA_dist: int = 14
    PSC_context: int = 15
    I_soma: int = 16
    

def init_Ca_AdEx_nestml(default_params, multi_comp=False):
    if multi_comp:
        soma_params = {
            'C_m': default_params['C_m_s'],                   # [pF] Soma capacitance
            'g_L': default_params['g_L_s'],                   # [nS] Soma leak conductance
            'e_L': default_params['e_L_s'],                   # [mV] Soma reversal potential
            'g_spike': default_params['g_L_s'],          # [nS] Adex conductance
            'v_thr': default_params['e_Na_Adex'],         # [mV] Adex threshold
            'delta_thr': default_params['delta_T'],              # [mV] Adex slope factor
            # 'd_BAP': default_params['d_BAP'],             # delay not yet implemented
            'refr_period': default_params['t_ref'] + 0.2,
            'v_reset': default_params['V_reset'],
            'tau_w': default_params['tau_w'],
            'sth_a': default_params['a'],
            'b': default_params['b'],
            'e_adapt': default_params['e_L_s'],
            'g_adapt': default_params['g_w'],
            'g_refr': 1e6,
        }

        distal_params = {
            'C_m': default_params['C_m_d'],                 # [pF] Distal capacitance
            'g_L': default_params['g_L_d'],                 # [nS] Distal leak conductance
            'g_C': default_params['g_C_d'],                 # [nS] Soma-distal coupling conductance
            'e_L': default_params['e_L_d'],                 # [mV] Distal reversal potential
            'gbar_Ca': default_params['gbar_Ca'],           # [nS] Ca maximal conductance
            'gbar_K': default_params['gbar_K_Ca'],       # [nS] K_Ca maximal conductance
            'e_K': default_params['e_K'],                   # [mV] K reversal potential
            'tau_Ca': default_params['tau_decay_Ca'], # [ms] decay of Ca concentration
            'phi_Ca': default_params['phi'],                   # [-] scale factor
            'm_half_Ca': default_params['m_half'],             # [mV] m half-value for Ca
            'h_half_Ca': default_params['h_half'],             # [mV] h half-value for Ca
            'm_slope_Ca': default_params['m_slope'],           # [-] m slope factor for Ca
            'h_slope_Ca': default_params['h_slope'],           # [-] h slope factor for Ca
            'tau_m_Ca': default_params['tau_m'],               # [ms] m tau decay for Ca
            'tau_h_Ca': default_params['tau_h'],               # [ms] h tau decay dor Ca
            'tau_m_K': default_params['tau_m_K_Ca'],     # [ms] m tau decay for K_Ca
            'c_Ca0': default_params['Ca_0'],                 # [mM] Baseline intracellular Ca conc
            'c_Ca_thr': default_params['Ca_th'],               # [mM] Threshold Ca conc for Ca channel opening
            'exp_K_Ca': default_params['exp_K_Ca'],
            'w_bap': default_params['w_BAP'],          # [-] Exponential factor in K_Ca current with Hay dyn
            'g_adapt': 0.,
            'g_spike': 0.,
            'g_refr': 0.,
        }
        other_params = {'V_th': default_params['V_th']}
    else:
        soma_params = {
            'C_m': default_params['C_m'],                   # [pF] Soma capacitance
            'g_L': default_params['g_L'],                   # [nS] Soma leak conductance
            'e_L': default_params['E_L'],                   # [mV] Soma reversal potential
            'g_spike': default_params['g_L'],          # [nS] Adex conductance
            'v_thr': default_params['V_th'],         # [mV] Adex threshold
            # 'v_thr': default_params['V_peak'],         # [mV] Adex threshold
            'delta_thr': default_params['Delta_T'],              # [mV] Adex slope factor
            # 'd_BAP': default_params['d_BAP'],             # delay not yet implemented
            'g_refr': 1e6,
            'refr_period': default_params['t_ref'] + 0.2,
            'v_reset': default_params['V_reset'],
            'tau_w': default_params['tau_w'],
            'sth_a': default_params['a'],
            'b': default_params['b'],
            'e_adapt': default_params['E_L'],
            'g_adapt': 1., # educated guess
            'g_refr': 1e6,
        }
        distal_params = None
        other_params = {'V_th': default_params['V_peak']}

        
    return(soma_params, distal_params, other_params)


def create_receptor_mapping(multi_comp=True):        
    if multi_comp:        
        # receptor mapping
        ri = ReceptorIdxs(
            # receptors get assigned an index which corresponds to the order in which they
            # are added. For clearer bookkeeping, we explicitly define these indices here.
            ALPHAexc_soma=0, 
            ALPHAinh_soma=1, 
            ALPHAexc_dist=2, 
            ALPHAinh_dist=3, 
            AMPA_NMDA_dist=4, 
            PSC_context=5, 
            I_soma=6,
            # other receptors types are not used, assign something that will generate an error when used
            GABA_soma=float('nan'), 
            AMPA_soma=float('nan'), 
            GABA_dist=float('nan'), 
            AMPA_dist=float('nan'), 
            REFRACT_soma=float('nan'), 
            ADAPT_soma=float('nan'), 
            BETA_dist=float('nan'), 
            NMDA_soma=float('nan'), 
            AMPA_NMDA_soma=float('nan'), 
            NMDA_dist=float('nan'), 
        )
    else:
        # receptor mapping
        ri = ReceptorIdxs(
            # receptors get assigned an index which corresponds to the order in which they
            # are added. For clearer bookkeeping, we explicitly define these indices here.
            ALPHAexc_soma=0, 
            ALPHAinh_soma=1, 
            AMPA_NMDA_soma=2,
            PSC_context=3, 
            I_soma=4,
            # other receptors types are not used, assign something that will generate an error when used
            GABA_soma=float('nan'), 
            AMPA_soma=float('nan'), 
            GABA_dist=float('nan'), 
            AMPA_dist=float('nan'), 
            REFRACT_soma=float('nan'), 
            ADAPT_soma=float('nan'), 
            BETA_dist=float('nan'), 
            NMDA_soma=float('nan'), 
            ALPHAexc_dist=float('nan'), 
            ALPHAinh_dist=float('nan'), 
            AMPA_NMDA_dist=float('nan'), 
            NMDA_dist=float('nan'), 
        )
    return ri


def create_nestml_neuron(n = 1, params = {}, multi_comp = True, neuron_model="ca_adex_2expsyn_nestml"):
    soma_params, distal_params, other_params = init_Ca_AdEx_nestml(params, multi_comp=multi_comp)
    # initialize the model with two compartments
    cm = nest.Create(neuron_model, n)
    cm.V_th = other_params['V_th']

    print(f"? multi_comp: {multi_comp}")

    if multi_comp:
        cm.compartments = [
            {"parent_idx": -1, "params": soma_params},
            {"parent_idx":  0, "params": distal_params}
        ]
        # receptors
        receptors = [
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": .173, "tau_d_syn": .227, "e_syn": 0.}}, # ALPHAexc_soma
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": 1.73, "tau_d_syn": 2.27, "e_syn": -85.}}, # ALPHAinh_soma
            {"comp_idx": 1, "receptor_type": "syn_2exp", "params": {"tau_r_syn": .173, "tau_d_syn": .227, "e_syn": 0.}}, # ALPHAexc_dist
            {"comp_idx": 1, "receptor_type": "syn_2exp", "params": {"tau_r_syn": 1.73, "tau_d_syn": 2.27, "e_syn": -85.}}, # ALPHAinh_dist
            {"comp_idx": 1, "receptor_type": "ampa_nmda"}, # AMPA_NMDA_dist
            {'comp_idx': 1, "receptor_type": "syn_psc", "params": {"tau_r_syn": 0.173, "tau_d_syn": 10.}},
            {"comp_idx": 0, "receptor_type": "inp"}, # I_SOMA
        ]
        cm.receptors = receptors

    else:
        cm.compartments =[
            {"parent_idx": -1, "params": soma_params},
        ]
        # receptors
        receptors = [
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": .173, "tau_d_syn": .227, "e_syn": 0.}}, # ALPHAexc_soma
            {"comp_idx": 0, "receptor_type": "syn_2exp", "params": {"tau_r_syn": 1.73, "tau_d_syn": 2.27, "e_syn": -85.}}, # ALPHAinh_soma
            {"comp_idx": 0, "receptor_type": "ampa_nmda"}, # AMPA_NMDA_soma
            {'comp_idx': 0, "receptor_type": "syn_psc", "params": {"tau_r_syn": 0.173, "tau_d_syn": 10.}},
            {"comp_idx": 0, "receptor_type": "inp"}, # I_SOMA
        ]
        cm.receptors = receptors

    return cm
