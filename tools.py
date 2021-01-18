#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:35:12 2021

@author: roxane
"""
import numpy as np
import sys,os
from scipy.io import FortranFile
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def create_config_file(path, simuname, fric_law = 'RateStateAgeing_R', mu=30*1e9, sigma_N=-1e8, Dc=1e-3, b=0.01, a_over_b=0.75, V_val_x1=-10, V_val_x2=10, V_val_pourc=0.001, sigma11_dot_inf = 0.00e+00, sigma22_dot_inf = 0.00e+00, sigma12_dot_inf = 0.00e+00, sigma31_dot_inf = 0.00e+00, sigma32_dot_inf = 0.00e+00, stop_criteria=1, max_it=10000, final_time=10000000) :  
    # Preliminary checks
    if stop_criteria not in [0, 1, 2]:
        raise TypeError('stop_criteria should be 0, 1 or 2')
    
    if V_val_x1 > V_val_x2:
        raise ValueError('V_val_x1 should be inferior to V_val_x2')
    
    # Computes a and Lb
    a = a_over_b * b
    Lb = - (mu * Dc) / (sigma_N * b)
    
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
     '/\n',
     '&fault_friction\n',
     "a_distr = 'CST'\n",
     'a_val = -1.000e+00,-1.000e+00,{},{}\n'.format(a, a),
     "b_distr = 'CST'\n",
     'b_val = -1,-1,{},{}\n'.format(b, b),
     "Dc_distr = 'CST'\n",
     'Dc_val = -1,-1,{},{}\n'.format(Dc, Dc),
     "V0_distr = 'CST'\n",
     'V0_val = -1.000e+00,-1.000e+00,1e-09,1e-09\n',
     '/\n',
     '&fault_initial_condition\n',
     "slip_distr = 'CST'\n",
     'slip_val = -1.000e+00,-1.000e+00,0,0\n',
     "V_distr = 'CST'\n",
     'V_val = {},{},1e-09,{}e-09\n'.format(V_val_x1, V_val_x2, 1 + V_val_pourc),
     "theta_distr = 'CST'\n",
     'theta_val = -1.000e+00,-1.000e+00,1000000.0,1000000.0\n',
     "sigmaN_distr = 'CST'\n",
     'sigmaN_val = -1,-1,{},{}\n'.format(sigma_N, sigma_N),
     'ds = {}\n'.format(Lb/10),
     '/\n',
     '&simulation_parameter\n',
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


def create_geom_1fault(path, x1, x2):
    content = ['{} 0.000e+00 \n'.format(x1), 
               '{} 0.000e+00 \n'.format(x2), 
               '/\n']
   
    with open(path + 'geometry.in', 'w') as file:
        for line in content:
            file.write(line)
            
    return


def create_geom_2_faults_overlapping(path, Lnuc, D_over_Lnuc, L_over_Lnuc, overlap):
    """
    Create geometry.in file for a 2 faults overlapping geometry (as in Romanet et al. (2018), GRL).
 
         L/Lnuc * overlap
         <--->   
         ==========
          | D/Lnuc  
    ==========
    <-------->
      L/Lnuc

    Parameters
    ----------
    path : string
        Path to the directory.
    L_nuc : real
        Nucleation length.
    D_over_Lnuc : real
        Distance between the two fault express as the ratio D/Lnuc.
    L_over_Lnuc : real
        Length of the faul express as the ratio L/Lnuc.
    overlap : real
        Portion of the fault overlapping. 0 < overlap < 1.

    Returns
    -------
    None.

    """
    content = ['0.000e+00 0.000e+00 \n', 
               '{} 0.000e+00 \n'.format(L_over_Lnuc * Lnuc), 
               '/\n',
               '{} {} \n'.format(L_over_Lnuc * Lnuc * (1 - overlap), D_over_Lnuc * Lnuc),
               '{} {} \n'.format(L_over_Lnuc * Lnuc * (2 - overlap), D_over_Lnuc * Lnuc),
               '/\n']
       
    with open(path + 'geometry.in', 'w') as file:
        for line in content:
            file.write(line)
            
    return


def create_GPS_file(path, x1, x2):
    """
    Creates GPS.in at the location 'path'.

    Parameters
    ----------
    path : string
        Path to the directory.
    x1 : list
        First coordinate of the GPS station.
    x2 : list
        Second cooridnate of the GPS station.

    Returns
    -------
    None.

    """
    if len(x1) != len(x2) :
        print('x1 and x2 are not of the same size !')
    else:    
        content = []
        for i, el in enumerate(x1):
            content.append('{} {} \n'.format(el, x2[i]))  
        content.append('/')
       
        with open(path + 'GPS.in', 'w') as file:
            for line in content:
                file.write(line)
            
    return       


def create_tides_file(path, a1, a2, a3):
    """
    Creates tides.in at the location 'path'.


    """
    if len(a1) != len(a2) != len(a3) :
        print('a1, a2 and a3 are not of the same size !')
    else:    
        content = []
        for i, el in enumerate(a1):
            content.append('{} {} {} \n'.format(el, a2[i], a3[i]))  
        content.append('/')
       
        with open(path + 'tides.in', 'w') as file:
            for line in content:
                file.write(line)
            
    return   


class ReadBinaryFiles:
    def __init__(self, path):
        
        self.path = path
        
        # Reads geometry.out
        with open(path + 'geometry.out', 'r') as file:
            content = file.readlines()
        
        # Extract number of fault from the first line
        s1, s2 = content[0].split('         ')
        nbr_fault, s3 = s2.split('\n')
        self.nbr_fault = int(nbr_fault)
        
        # Extract number of elements per fault 
        self.elements = []  # stores number of elements for each fault 
        for i in range(self.nbr_fault):
            s1, s2, s3 = content[i + 1].split('       ')
            elements, s1 = s3.split('\n')
            self.elements.append(int(elements))
        
        # Total number of elements
        # Each fault is separated by a fake element 
        self.nel = sum(self.elements) + (nbr_fault - 1)
        
        # Reads files in the directory 'path'
        self.files_in_dir = os.listdir(self.path)  # Lists all the files
        
        # Run read_file to read time and velocity files
        self.read_file('time')
        self.read_file('velocity')
        
        
    def read_file(self, file_type):    

        files = []  # list of files of the right type
        # Loops over the files and select those with 'file_type' in their name
        for file_name in self.files_in_dir:
            if file_type in file_name:
                files.append(file_name)
                
        # Sorts selected files in ascending order 
        files.sort()
        
        data = [] # list to store the data 
        # Loops over selected files and read them
        for file in files:
            f = FortranFile(self.path + file,'r')
            print(self.path)
            if file_type == 'velocity':
                content = f.read_reals(dtype='f8')
                data.append(content.reshape(int(len(content)/self.nel), self.nel))
                
            else:
                data.append(f.read_reals(dtype='f8'))
        
        # Concatenate data from all files         
        data = np.concatenate(data)        
        
        # Clear velocity data
        cum = np.cumsum(self.elements) 
        if file_type == 'velocity':
            for i in range(self.nbr_fault - 1):
                data = np.delete(data, cum[i], axis=1)  # Removes false column from data
        
        # Stores data
        if file_type == 'velocity':
            self.velocity = data
        elif file_type == 'time':
            self.time = data
            
            
def compute_max_vel(velocity):
    max_vel = []
    for i, el in enumerate(velocity):
        max_vel.append(np.max(el))
        
    return max_vel    


def plot_max_vel_comp(simunames, figname, nrows=3, ncols=3, savefig=True):
    # Reads time, velocity and plotting info
    times = [] # list to store time of all simulation
    max_velocities = [] # list to store velocity for all simulation
    sigma_plot = [] # list of sigma in the plotting order
    V_val_pourc_plot = [] # list of V_val_pourc in the plotting order
    
    for i, sim in enumerate(simunames):
        path = '/Users/roxane/Desktop/version11/problems/' + sim + '/'
        times.append(ReadBinaryFiles(path).time / (3600*24*362.25)) # in year
        max_velocities.append(compute_max_vel(ReadBinaryFiles(path).velocity))
        
        # Reads plotting info
        with open(path + 'which_params.txt', 'r') as file:
            content = file.readlines()
        s1, sig = content[2].split(' = ')
        sigma_plot.append(float(sig))
        s2, v = content[5].split(' = ')
        V_val_pourc_plot.append(float(v))
    
    ############
    ### Plot ###
    ############
    
    # Create figure 
    fig, axs = plt.subplots(nrows, ncols,figsize=(6, 5))
    
    idx = 0
    for i in range(nrows):
        for j in range(ncols):
            # Plot data
            axs[i, j].plot(times[idx], max_velocities[idx], color='red')
            
            # Set ylim
            axs[i, j].set_ylim(1e-14, 1)
            
            # Plot V_val_pourc legend 
            if i == 0:  # for the first line only
                axs[i, j].set_title('$\Delta_V$ ={}%'.format(V_val_pourc_plot[idx] * 100))  
            
            # Plot sigma legend 
            if j == nrows - 1:  # for the last column only
                # axs[i, j].set(ylabel=r'$\dot \sigma^\infty = {}$'.format(sigma_plot[idx]))
                axs[i, j].set_ylabel(r'$\dot \sigma^\infty = {}$'.format(sigma_plot[idx]), fontsize=11)
                axs[i, j].yaxis.set_label_position("right")
                # axs[i, j].yaxis.tick_right()
        
            # Fill between line and O
            axs[i, j].fill_between(times[idx], 0, max_velocities[idx], color='lightgrey')
            
            # Sets x-axis ticks 
            # if max(times[idx]) >= 3:
            #     axs[i, j].xaxis.set_ticks(np.arange(min(times[idx]), max(times[idx]+1), 1))
            # print('max times ={}'.format(max(times[idx])))    
            # if max(times[idx]) < 1:
            #     axs[i, j].xaxis.set_ticks(np.arange(min(times[idx]), max(times[idx]+0.02), 0.02))    
            
            # Set tick label size
            axs[i, j].yaxis.set_tick_params(labelsize=9)
            axs[i, j].xaxis.set_tick_params(labelsize=9)
    
            # Create grid
            axs[i, j].grid()
            
            idx += 1    
    
    # Common x label 'Time'
    axs[nrows - 1, int(np.floor(ncols / 2))].set(xlabel='Time (year)')
    axs[nrows - 1, int(np.floor(ncols /2))].xaxis.label.set_size(12)
    
    # Set yscale and log scale, avoir white space between start and en of axis
    for ax in axs.flat:
        ax.set_yscale('log')
        ax.set_yticks([1e-14, 1e-10, 1e-6, 1e-2])
        ax.margins(x=0) # sets white space between start and end of x-axis 
    
    # Remove y tick labels for all subplots except left column
    for ax in axs[:, 1:].flat:
        ax.yaxis.set_ticklabels([])

    # Add common ylabel     
    fig.text(0, 0.5, 'Maximum slip rate ($m.s^{-1}$)', va='center', rotation='vertical', fontsize = 12)
    
    # Remove whitespace aroung fig
    plt.tight_layout()
    
    # Saves
    if savefig == True:
        fig.savefig('/Users/roxane/Desktop/Roxane/Figures/' + figname + '.png', dpi = 400)
        
        
def prepare_simulations_1fault(simunames, V_val_x1, V_val_x2, V_val_pourc, b, a_over_b, L_over_Lnuc, sigma12_dot_inf, sigma32_dot_inf, Dc=1e-3, sigma_N=-1e8, mu=30*1e9, stop_criteria=1, max_it=10000, final_time=10000000): 
    # Computes Lnuc, edges of fault (x1 and x2) and L
    Lnuc = (mu * Dc) / (sigma_N * (b - a_over_b * b))
    
    # Create all the files (config.in etc.) for all the simulations 
    simu_idx = 0 # simuname index 
    for i, el in enumerate(L_over_Lnuc):
        x1 = Lnuc / 2 * el
        x2 = - Lnuc / 2 * el
            
        for j, el2 in enumerate(sigma12_dot_inf):
            path = '/Users/roxane/Desktop/version11/problems/' + simunames[simu_idx] + '/'
            # Creates the directory for the simulation
            try : 
                os.mkdir(path)
            except:
                pass   
            
            # Computes final time so that it is 15T
            T = (mu / el2) * (el)**(-5.607)
            final_time = 20*T
            
            # Creates files for the simulation
            create_tides_file(path, [5.0e+4], [86400.0], [0.000])
            create_GPS_file(path, [10], [10])
            create_geom_1fault(path, x1, x2)
            create_config_file(path, simunames[simu_idx], V_val_x1 = V_val_x1, V_val_x2 = V_val_x2, V_val_pourc = V_val_pourc, sigma12_dot_inf = el2, sigma32_dot_inf = el2, stop_criteria=stop_criteria, max_it=max_it, final_time=final_time) 
            
            # Creates a texte file containing key parameters -- Just for a short overview 
            with open(path + 'which_params.txt', 'w') as file:
                file.write('Simuname : {} \n'.format(simunames[simu_idx]))
                file.write('sigma12_dot_inf = {} \n'.format(el2))
                file.write('sigma32_dot_inf = {} \n'.format(el2))
                file.write('V_val_x1 = {} \n'.format(V_val_x1))
                file.write('V_val_x2 = {} \n'.format(V_val_x2))
                file.write('V_val_pourc = {} \n'.format(el))
                file.write('max_it = {} \n'.format(max_it))
                
            simu_idx += 1
            

def run_simulations(simunames):
    os.chdir('/Users/roxane/Desktop/version11/')
    os.system('source ~/.zshrc')
    
    with open('/Users/roxane/Desktop/version11/runproblem', 'w') as file:
        file.write('source /opt/intel/compilers_and_libraries_2020.1.216/mac/bin/compilervars.sh intel64\n')
        for simuname in simunames:
            file.write('./fastcycles ./problems/' + simuname + '/ &')
    
    os.system('/Users/roxane/Desktop/version11/runproblem')