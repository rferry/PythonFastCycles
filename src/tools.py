"""
Script of usefull tools.

First version by R. Ferry on January 2021.
"""
import numpy as np
import sys
import os
import re
from scipy.io import FortranFile
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib import cm 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


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
    top = cm.get_cmap('hot_r', 128)  # r means reversed version
    top = ListedColormap(top(np.linspace(0.3, 0.8, 256)))
        
    bottom = cm.get_cmap('Blues', 128)
    bottom = ListedColormap(bottom(np.linspace(0, 0.8, 256)))

    middle = cm.get_cmap('viridis', 256)
    middle = ListedColormap(middle(np.linspace(0.3, 1, 256)))
    
    # Combine the colormaps    
    newcolors = np.vstack((bottom(np.linspace(0, 1, p3)), 
                           middle(np.linspace(0, 1, p2)), 
                           top(np.linspace(0, 1, p1))))
     
    # Create the final colormap
    cmp = ListedColormap(newcolors, name='nicecmp')
    
    return cmp

def clean_simulations(path, simunames):
    """
    Clean simulations' folders by removing extra theta files left from previous
    runs.

    Parameters
    ----------
    path : string
        Path to the simulations's folder (typically 'problems'). Must have '/'
        at the end.
    simunames : list of sting
        Names of the simulations to clean.

    """
    for i, simuname in enumerate(simunames):
        # List all files in simulation directory
        files_in_dir = os.listdir(path + simuname)
        
        # Initialisation
        theta_files = []
        time_files = []
        # Loops over the files and select theta and time files 
        for file_name in files_in_dir:
            if 'theta' in file_name:
                theta_files.append(file_name)
            elif 'time' in file_name:
                time_files.append(file_name)
        
        # Check if both list are of the same size. If yes, there is no problem
        if len(theta_files) == len(time_files):
            pass
        else:
            # Get numbers of theta files
            num_theta = [file[5:] for file in theta_files]
            # Get numbers of time files
            num_time = [file[4:] for file in time_files]
            
            # Find theta files with no corresponding time files
            res = [value for value in num_theta if value not in num_time]
            
            # Add them to the files to delete 
            files_to_delete = ['theta' + value for value in res]
            
            # Delete the files 
            for file in files_to_delete:
                os.remove(path + simuname + '/' + file)

class FixSimulations:
    """
    Class to fix simulations' problems. 
    """
    def __init__(self, path, simunames, verbose=True):
        """
        Initialisation of the class.

        Parameters
        ----------
        path : string
            Path to the folder with simulations (called 'problems' most of the 
            times). Must have '/' at the end.
        simunames : list of strings
            Names of the simulations that need to be check.
        verbose : bool, optional
            If True, will speak to you and ask to replace tol_solver and 
            relaunch simulations. If False, you will need to use the function 
            'replace_tol_solver' and to relaunch simulations by yourself. The 
            default is True.

        Returns
        -------
        None.

        """
        self.path = path
        self.simunames = simunames
        
        # Find simulations that need to be fixed
        self.check_simulations()
        
        if verbose:
            # if some simulations have a problem
            if len(self.simu_affected) != 0:  
                line = 'Simulations with -Infinity values are: \n'
                for i, simu in enumerate(self.simu_affected):
                    line += '- {} \n'.format(simu)
                print(line)
                ans = input('Do you want to replace tol_solver ? (y/n)')
                if ans == 'y':
                    self.replace_tol_solver()
                    ans2 = input('Do you want to relaunch simulations ? (y/n)')
                    if ans2 == 'y':
                        run_simulations(self.simu_affected)
                    elif ans2 == 'n':
                        pass
                    else:
                        print('Your answer is not valid.')
                elif ans == 'n':
                    pass
                else:
                    print('Your answer is not valid.')
            # If no simulations need to be changed         
            else:
                print('No problem detected !')
        
    def check_simulations(self):
        """
        Read MomentRate.out for all simulations in simunames and find those 
        with -Infinity values.
        The names of the simulations that need to be fixed will be stored in
        self.simu_affected.

        """
        # Find simulations with a problem
        simu_affected = []
        for i, simuname in enumerate(self.simunames):
            # Check if MomentRate.out contains -infinity
            with open(self.path + simuname + '/MomentRate.out', 'r') as file:
                content = file.read()
                
            if '-Infinity' in content:
                simu_affected.append(simuname)
        
        # Store results 
        self.simu_affected = simu_affected 
        
    def replace_tol_solver(self):
        """
        Replace tol_solver in config.in by a value of one order of magnitude 
        less for ll 'simu_affected' found by check_simulations.

        """
        for i, simuname in enumerate(self.simu_affected):
            # Read config.in
            with open(self.path + simuname + '/config.in', 'r') as file:
                content = file.readlines()
            
            # Loop over config.in
            for j, line in enumerate(content):
                # Find tol_solver
                if line.startswith('tol_solver'):
                    # Read the value
                    tol = float(line.split()[2])
                    # Replace tol_solver 
                    content[j] = 'tol_solver = {} \n'.format(str(tol/10))
            
            # Write new config.in
            with open(self.path + simuname + '/config.in', 'w') as file:
                for line in content:
                    file.write(line)
                    
def find_simunames(path, regex):
    """
    Find all simulations in the directory "path" wich names correspond to the 
    regular expression "regex". 
    For a guide on regular expressions see:
    https://docs.python.org/3/howto/regex.html#regex-howto

    Parameters
    ----------
    path : string
        Path to the folder with simulations (called 'problems' most of the 
        times). Must have '/' at the end.
    regex : string
        Expression to search expressed in regex format.
        Most common expressions are:
            ^a -> look for "a" at the beginning of the text
            a$ -> look for a at the end of the text
            a|b -> look for "a" OR "b"
            [abc] -> look for "a", "b" or "c"
                ex: "m[abc]s" will look for "mas", "mbs" or "mcs"
            [^abc] -> look for everything except "a", "b" and "c"   
            [a-z] -> look for all lowercase letters from a to z
            [0-9] -> lookf for all digits
            . -> look for any caracter 
            \w -> look for all letters and digits, lower or upper case 
            \d -> equivqlent to [0-9]
            {min,max} -> look for n repetitions, n between min and max 
            {,max} -> look for a maximum number of repetition
            {min,} -> look for a minimum number of repetition
            {number} -> look for an exact number of repetition
                ex: "[0-9]{2,4}" will look for 2 to 4 consecutive digit
            * -> equivalent to {0,}
            + -> equivalent to {1,}
            ? -> equivalent to {,1}
            \ -> remove the special meaning of a caracter (for example, if you
                 want to look for a point, use "\.")            

    Returns
    -------
    simus : list of string
        All the simulations matching the expression.

    """
    # List all files in simulations directory
    files_in_dir = os.listdir(path)
    
    # Initialisation
    simus = []
    
    # Loop over all simulation
    for simu in files_in_dir:
        # Check if simulation name match the expression
        obj = re.match(regex, simu)
        # if there is a match, save the name of the simulations
        if obj != None:
            simus.append(simu)
            
    # Return the list of all matching simulations
    return simus
