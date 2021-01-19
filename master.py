#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 14:13:19 2021

@author: Roxane
"""
# from IPython import get_ipython
# get_ipython().run_line_magic('matplotlib', 'qt')
import numpy as np
import sys,os
from scipy.io import FortranFile
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
# os.chdir('/Users/roxane/Desktop/Roxane')
sys.path.append('/Users/roxane/Desktop/Roxane/PythonFastCycles')
from tools import *
############################
# Defines constant varibles
############################
mu = 30*1e9 # Pa
Dc = 1e-3
sigma_N = -1e8
b = 0.01
a_over_b = 0.75
L_over_Lnuc = 2

Lnuc = (mu * Dc) / (sigma_N * (b - a_over_b * b))

x1 = Lnuc / 2 * L_over_Lnuc 
x2 = - Lnuc / 2 * L_over_Lnuc
    
L = x2 -x1

#%% Run all

V_val_x2 = L/4
V_val_pourc = [0.001, 0.01, 0.1] # Pourcentage of variation

HugeSigma = [1.0, 5.00e-01, 1.00e-01]
SmallSigma = [1.00e-01, 1.00e-02, 1.00e-03]
# list_sigmas = [HugeSigma, SmallSigma]
list_sigmas = [SmallSigma]
# sigmas_names = ['HugeSigma', 'SmallSigma']
sigmas_names = ['SmallSigma']
    
# V_val_x_list = [-L/4, -L/6, 0, L/6, L/4]
V_val_x1_list = [-L/4]
# V_val_legend = ['L_over_4', 'L_over_3', 'L_over_2', '2L_over_3', '3L_over_4']
V_val_legend = ['L_over_2']

for i, sigma in enumerate(list_sigmas):
    for j, V_val_x1 in enumerate(V_val_x1_list):
        # Creates simunames
        abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        simunames = []
        for k in abc:
            simunames.append(sigmas_names[i] + '_'+ V_val_legend[j] +'_' + k)
        
        
        # Run simulations
        prepare_simulations_1fault(simunames, V_val_x1, V_val_x2, V_val_pourc, b, a_over_b, L_over_Lnuc, sigma, sigma, max_it=1000000, stop_criteria=1, final_time=157788000.0)
        run_simulations(simunames)
#%% Plot all
for i, sigma in enumerate(list_sigmas):
    for j, V_val_x in enumerate(V_val_x1_list):
        # Creates simunames
        abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        simunames = []
        for k in abc:
            simunames.append (sigmas_names[i] + '_'+ V_val_legend[j] +'_' + k)        
            
        succed = False
        while not succed:
            try:
                plot_max_vel_comp(simunames, sigmas_names[i] + '_'+ V_val_legend[j], 3, 3)
                succed = True
            except :
                continue    
      
#%% 
def create_config_file2(path, simuname, fric_law = 'RateStateAgeing_R', mu=30*1e9, sigma_N=-1e8, Dc=1e-3, b=0.01, a_over_b=0.75, V_val_x1=-10, V_val_x2=10, V_val_pourc=0.001, sigma11_dot_inf = 0.00e+00, sigma22_dot_inf = 0.00e+00, sigma12_dot_inf = 0.00e+00, sigma31_dot_inf = 0.00e+00, sigma32_dot_inf = 0.00e+00, stop_criteria=1, max_it=10000, final_time=10000000, nf=False) :  
    # Preliminary checks
    if stop_criteria not in [0, 1, 2]:
        raise TypeError('stop_criteria should be 0, 1 or 2')
    
    if V_val_x1 > V_val_x2:
        raise ValueError('V_val_x1 should be inferior to V_val_x2')
    
    # Computes a and Lb
    a = a_over_b * b
    Lb = - (mu * Dc) / (sigma_N * b)
    
    # Get fault number if unspecified 
    if nf == False:
        try:
            with open(path + 'geometry.in') as file:
                content = file.read()
            nf = content.count('/')
        except Exception:  # if geometry.in do not exist
            print('You must specify a fault number or create geometry.in first')
        
    # Content of config.in file 
    content = ['&main\n',
     "simulation_name = '{}'\n".format(simuname),
     '/\n',
     '&friction\n',
     "fric_law = '{}'\n".format(fric_law),
     '/\n',
     '&material_and_loading\n',
     'mu = {} \n'.format(mu),
     'cp = 5.00e+03 \n',
     'cs = 3.50e+03 \n',
     'sigma11_dot_inf = {} \n'.format(sigma11_dot_inf),
     'sigma22_dot_inf = {} \n'.format(sigma22_dot_inf),
     'sigma12_dot_inf = {} \n'.format(sigma12_dot_inf),
     'sigma31_dot_inf = {} \n'.format(sigma31_dot_inf),
     'sigma32_dot_inf = {} \n'.format(sigma32_dot_inf),
     '/\n']

    for i in range(nf):
        content += [ '&fault_friction\n', 
                    "a_distr = 'CST'\n",
                    'a_val = -1.000e+00,-1.000e+00,{},{}\n'.format(a, a),
                    "b_distr = 'CST'\n",
                    'b_val = -1,-1,{},{}\n'.format(b, b),
                    "Dc_distr = 'CST'\n",
                    'Dc_val = -1,-1,{},{}\n'.format(Dc, Dc),
                    "V0_distr = 'CST'\n",
                    'V0_val = -1.000e+00,-1.000e+00,1e-09,1e-09\n',
                    '/\n']
    
    content += ['&fault_initial_condition\n',
     "slip_distr = 'CST'\n",
     'slip_val = -1.000e+00,-1.000e+00,0,0\n',
     "V_distr = 'CST'\n",
     'V_val = {},{},1e-09,{}e-09\n'.format(V_val_x1, V_val_x2, 1 + V_val_pourc),
     "theta_distr = 'CST'\n",
     'theta_val = -1.000e+00,-1.000e+00,1000000.0,1000000.0\n',
     "sigmaN_distr = 'CST'\n",
     'sigmaN_val = -1,-1,{},{}\n'.format(sigma_N, sigma_N),
     'ds = {}\n'.format(Lb/10),
     '/\n']
    
    for i in range(nf - 1):
        content += ['&fault_initial_condition\n',
                    "slip_distr = 'CST'\n",
                    'slip_val = -1.000e+00,-1.000e+00,0,0\n',
                    "V_distr = 'CST'\n",
                    'V_val = -1.000e+00,-1.000e+00,1e-09,1e-09\n',
                    "theta_distr = 'CST'\n",
                    'theta_val = -1.000e+00,-1.000e+00,1000000.0,1000000.0\n',
                    "sigmaN_distr = 'CST'\n",
                    'sigmaN_val = -1,-1,{},{}\n'.format(sigma_N, sigma_N),
                    'ds = {}\n'.format(Lb/10),
                    '/\n']
      
        
    content += ['&simulation_parameter\n',
     "fracture_mode = 'modeIII'\n",
     "static_kernel = 'hmatrix'\n",
     'initial_time = 0.000e+00 \n',
     'time_step_max = 10\n',
     'final_time = {}\n'.format(final_time),
     'stop_criteria = {}\n'.format(stop_criteria),
     'cut_vel = 1.000e-08 \n',
     'max_it = {} \n'.format(max_it),
     'tol_solver = 1.000e-08 \n',
     'tol_interp = 1.000e-08 \n',
     'iprec = 4\n',
     'icheck_interp = 0\n',
     'omp_threads = 4\n',
     '/\n',
     '&output\n',
     'stride_time = 5\n',
     'GPS_stride = 1\n',
     'isave_node = 1\n',
     'freq_writing_file = 1000\n',
     '/']

    # Create config.in
    with open(path + 'config.in', 'w') as file:
        for line in content:
            file.write(line)
            
    return 

#%%
plt.close("all")
x = np.arange(0, 128, 1)

a = 1
b = 100
l = 5
T0 = 50
# y = np.tanh((T0 - x)/2 +1)*A0      
y = 0.5 * ((a + b) + (a - b) * np.tanh((x - T0 + 4*l) / l))
# x = x + 2*l
plt.figure()
plt.plot(x, y) 

#%%
def find_period(time, max_vel):
    # Find peaks in the signal
    peaks, _ = find_peaks(max_vel)
    
    # Remove first peaks 
    peaks = peaks[5:]
    
    # Get time associated to peaks
    time_peaks = time[peaks]
    
    # Compute time difference between two peaks
    diff = time_peaks[1:] - time_peaks[:-1]
    
    # Compute average period in seconds
    period = np.mean(diff)
    
    return period
    
#%%
sdot = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]
list_L_over_Lnuc = [1, 2, 3, 4]

# Creates simunames
simunames = []
for k in range(36):
    simunames.append('simu{}'.format(k+1))


prepare_simulations_1fault(simunames, -20, 20, 0.01, 0.01, 0.75, list_L_over_Lnuc, sdot, sdot, Dc=1e-3, sigma_N=-1e8, mu=30*1e9, stop_criteria=2, max_it=1000000)
run_simulations(simunames)

#%%
path = '/Users/roxane/Desktop/version11/problems/'

period = []
for i, simu in enumerate(simunames):
    
    simulation = ReadBinaryFiles(path + simu +'/')
    simulation.max_vel = compute_max_vel(simulation.velocity)
    period.append(find_period(simulation.time, simulation.max_vel))
    
sdot_simu = sdot * 4
L_simu = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4]

#%%
mask = [~np.isnan(a) for a in period]
period = np.array(period)
period = period / (365.25*3600*24)
sdot_simu = np.array(sdot_simu)
L_simu = np.array(L_simu)

fig, ax = plt.subplots(1, 1)
scatter = ax.scatter(period[mask], sdot_simu[mask], c = L_simu[mask])

handles, labels = scatter.legend_elements(prop="colors", alpha=0.6)
legend2 = ax.legend(handles, labels, loc="best", title=r'$L/L_{nuc}$')
            