"""
Script of usefull tools.

First version by R. Ferry on January 2021.
"""
import numpy as np
import sys
import os
from scipy.io import FortranFile
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap


def run_simulations(simunames, path_to_makefile='/Users/roxane/Desktop/version11/'):
    os.chdir(path_to_makefile)
    os.system('source ~/.zshrc')

    with open(path_to_makefile + 'runproblem', 'w') as file:
        file.write(
            'source /opt/intel/compilers_and_libraries_2020.1.216/mac/bin/compilervars.sh intel64\n')
        for simuname in simunames:
            file.write('./fastcycles ./problems/' + simuname + '/ &')

    os.system(path_to_makefile + 'runproblem')
    
def nice_colormap(minv=1e-12, ssel=1e-8, eql=1e-3, maxv=1):  
    # Compute proportion of each colormap
    l = abs(np.log(minv) - np.log(maxv))
    l1 = abs(np.log(maxv) - np.log(eql)) / l
    l2 = abs(np.log(eql) - np.log(ssel)) / l
    p1 = int(round(384 * l1))
    p2 = int(round(384 * l2))
    p3 = 384 - p1 - p2
    
    # Define three colormaps
    top = cm.get_cmap('hot_r', 128) # r means reversed version
    top = ListedColormap(top(np.linspace(0.3, 0.8, 256)))
        
    bottom = cm.get_cmap('Blues', 128)
    bottom = ListedColormap(bottom(np.linspace(0, 0.8, 256)))

    middle = cm.get_cmap('viridis', 256)
    middle = ListedColormap(middle(np.linspace(0.3, 1, 256)))
    
    # Combine the colormaps    
    newcolors = np.vstack((bottom(np.linspace(0, 1, p3)), middle(np.linspace(0, 1, p2)),top(np.linspace(0, 1, p1))))
     
    # Create the final colormap
    cmp = ListedColormap(newcolors, name='nicecmp')
    
    return cmp
