"""
Class to read data. 

First version by R.Ferry on January 2021.
"""
import numpy as np
import sys
import os
from scipy.io import FortranFile
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import matplotlib.colors as colors 

from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

from mpl_toolkits.axes_grid1 import ImageGrid

import matplotlib.ticker as ticker


from tools import nice_colormap

class ReadData:
    def __init__(self, path):
        """
        Initialisation of ReadData class. Call functions to read data.

        Parameters
        ----------
        path : string
            Path to simulation directory. 

        Returns
        -------
        None.

        """
        # Store the path
        self.path = path
        
        # Read config.in
        self.read_config()
        
        # Read geometry.out
        self.read_geometry_out()

        # Total number of elements
        # Each fault is separated by a fake element
        self.nel = sum(self.elements) + (self.nbr_fault - 1)

        # Reads files in the directory 'path'
        self.files_in_dir = os.listdir(self.path)  # Lists all the files

        # Run read_file to read time, velocity, tractionN and theta files
        self.read_binary_file('time')
        self.read_binary_file('velocity')
        self.read_binary_file('theta')
        self.read_binary_file('tractionN')

        # Read MomentRate.out
        self.read_moment_rate()

        # Compute max_vel
        self.compute_max_vel()
        
        # Compute L and Lnuc
        self.L = []
        self.Lnuc = []
        for i in range(self.nbr_fault):
            self.L.append(self.elements[i] * self.ds[i])
            self.Lnuc.append((self.mu * self.Dc[i]) / (self.sigmaN[i] * (self.b[i] - self.a[i])))

    def read_binary_file(self, file_type):
        """
        Read Fortran binary files "velocity", "tractionN", "theta" and "time".

        Parameters
        ----------
        file_type : string
            File type to read. It can be "velocity", "tractionN", "theta" or 
            "time".

        Returns
        -------
        None.

        """

        files = []  # list of files of the right type
        # Loops over the files and select those with 'file_type' in their name
        for file_name in self.files_in_dir:
            if file_type in file_name:
                files.append(file_name)

        # Sorts selected files in ascending order
        files.sort()

        data = []  # list to store the data
        data_temp = []  # list to the store temporary data for velocity
        # Loops over selected files and read them
        for file in files:
            f = FortranFile(self.path + file, 'r')
            print(self.path)
            if file_type in ['velocity', 'tractionN', 'theta']:
                content = f.read_reals(dtype='f8')
                data_temp.append(content.reshape(
                    int(len(content)/self.nel), self.nel))

            else:
                data.append(f.read_reals(dtype='f8'))

        # Velocity, tractionN and theta case
        cum = np.cumsum(self.elements)
        if file_type in ['velocity', 'tractionN', 'theta']:
            # Concatenate data from all files
            # for velocity, tractionN and theta
            data_temp = np.concatenate(data_temp)

            # Cleaning data
            for i in range(self.nbr_fault - 1):
                # Removes false column from data
                data_temp = np.delete(data_temp, cum[i], axis=1)

            # Reshape to get data for each fault
            last_el = 0
            for i, el in enumerate(self.elements):
                data.append(data_temp[:, last_el:last_el+el])
                last_el += el

        else:  # if files are time files
            # Concatenate data from all files
            data = np.concatenate(data)

        # Stores data
        if file_type == 'velocity':
            self.velocity = data
        elif file_type == 'time':
            self.time = data
        elif file_type == 'theta':
            self.theta = data
        elif file_type == 'tractionN':
            self.tractionN = data

    def read_geometry_out(self):
        """
        Read "geometry.out" file to get the number of fault, the number of 
        elements per fault and the coordinates of the segments.

        Returns
        -------
        None.

        """
        # Open geometry.out and read content
        with open(self.path + 'geometry.out', 'r') as file:
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
            
        # Extract coordinates of elements
        ex_temp = []  # x coordinate of each elements
        ey_temp = []  # y coordinate of each elements 
        ex = []
        ey = []
        # Loop over content and skip headers
        for i, line in enumerate(content[2+self.nbr_fault:]):
            ex_temp.append(float(line.split()[1]))
            ey_temp.append(float(line.split()[2])) 
        # Reshape to get coordinates for each fault
        last_el = 0
        for i, el in enumerate(self.elements):
            ex.append(ex_temp[last_el:last_el+el])
            ey.append(ey_temp[last_el:last_el+el])
            last_el += el
        self.ex = ex
        self.ey = ey            
            
    def compute_max_vel(self):
        """
        Compute the maximum velocity for all faults.

        Returns
        -------
        None.

        """
        # TODO ! Add check if self.velocity is defined
        max_vel = []
        # Loop over the faults
        for i in range(self.nbr_fault):
            max_vel_temp = []
            # Loop over segments
            for j, el in enumerate(self.velocity[i]):
                max_vel_temp.append(np.max(el))
            max_vel.append(max_vel_temp)  # max_vel for one fault

        # Store results 
        self.max_vel = max_vel
        
    def read_moment_rate(self):
        """
        Read "MomentRate.out" file to get net moment rate and moment rate for 
        all faults.

        Returns
        -------
        None.

        """
        # Open and read MomentRate.out
        with open(self.path + 'MomentRate.out', 'r') as file:
            content = file.readlines()

        netMrate = []  # to store net moment rate
        Mrate = []  # to store moment rate for each fault
        tMrate = []  # to store moment rate time
        # stores Mrate data temporary
        Mrate_temp = np.ones((len(content)-2, self.nbr_fault))
        for i, line in enumerate(content[2:]):  # [2:] to skip headers
            tMrate.append(float(line.split()[1]))
            netMrate.append(float(line.split()[2]))
            for j in range(self.nbr_fault):  # Loop over faults
                Mrate_temp[i, j] = line.split()[3+j]
        for i in range(self.nbr_fault):
            Mrate.append(Mrate_temp[:, i])

        # Stores results
        self.netMrate = netMrate
        self.Mrate = Mrate
        self.tMrate = tMrate

    def read_config(self):
        """
        Read "config.in" file to get simulation parameters

        Returns
        -------
        None.

        """
        # Open config.in and read content
        with open(self.path + 'config.in', 'r') as file:
            content = file.readlines()
            
        # Initialisation
        a = []
        b = []
        Dc = []
        sigmaN = []
        ds = []
        sigma_dot = np.zeros((3, 3)) # sigma_dot components will be zero unless
                                     # specified otherwise 
        
        # Loop over content
        for i, line in enumerate(content):
            if line.startswith('a_val'):
                a.append(float((line.split()[2]).split(',')[2]))
            elif line.startswith('b_val'):
                b.append(float((line.split()[2]).split(',')[2])) 
            elif line.startswith('Dc_val'):
                Dc.append(float((line.split()[2]).split(',')[2])) 
            elif line.startswith('sigmaN_val'):
                sigmaN.append(float((line.split()[2]).split(',')[2]))
            elif line.startswith('ds'):
                ds.append(float(line.split()[2]))
            elif line.startswith('fric_law'):
                self.fric_law = line.split()[2][1:-1]
            elif line.startswith('mu'):
                self.mu = float(line.split()[2])
            elif line.startswith('sigma11_dot_inf'):
                sigma_dot[0, 0] = float(line.split()[2])
            elif line.startswith('sigma22_dot_inf'):
                sigma_dot[1, 1] = float(line.split()[2])
            elif line.startswith('sigma33_dot_inf'):
                sigma_dot[2, 2] = float(line.split()[2])
            elif line.startswith('sigma12_dot_inf') or line.startswith('sigma21_dot_inf'):
                sigma_dot[0, 1] = float(line.split()[2])  
                sigma_dot[1, 0] = float(line.split()[2])
            elif line.startswith('sigma31_dot_inf') or line.startswith('sigma13_dot_inf'):
                sigma_dot[0, 2] = float(line.split()[2])  
                sigma_dot[2, 0] = float(line.split()[2])
            elif line.startswith('sigma32_dot_inf') or line.startswith('sigma23_dot_inf'):
                sigma_dot[2, 1] = float(line.split()[2])  
                sigma_dot[1, 2] = float(line.split()[2])  
        
        # Store values 
        self.a = a
        self.b = b
        self.Dc = Dc
        self.sigmaN = sigmaN
        self.ds = ds
        self.sigma_dot = sigma_dot 
    
    def read_GPS_rate(self):
        # TODO !
        pass

    def read_EQcatalog(self):
        # TODO !
        pass

    def plot_max_vel(self, ssel=1e-8, eql=1e-3, start=0, stop=None, savefig=True):
        """
        Plot mamimum slip for all faults.

        Parameters
        ----------
        ssel : float, optional
            SSE velocity limit to draw horizontal line. The default is 1e-8.
        eql : float, optional
            Earthquake velocity limit to draw horizontal line. The default is 
            1e-3.
        start : int, optional
            Starting index to plot data. The default is 0.
        stop : int, optional
            Last index to plot data. The default is None (i.e. last data)

        Returns
        -------
        None.

        """
        # Determine yaxis lim --> determine min and max velocity
        minv = 1000  # ridiculous value to always be below
        maxv = -9999  # ridiculous value to always be above
        for i in range(self.nbr_fault):
            mint = np.min(self.max_vel[i][start:stop])
            maxt = np.max(self.max_vel[i][start:stop])
            if mint < minv:
                minv = mint
            else:
                pass
            if maxt > maxv:
                maxv = maxt
            else:
                pass

        # Time
        time = self.time[start:stop] / (3600*24*362.25)

        # Create figure
        fig, axs = plt.subplots(self.nbr_fault, 1, squeeze=False)

        for i in range(self.nbr_fault):
            # Plot data
            axs[i, 0].plot(time, self.max_vel[i][start:stop], color='red')

            # Set ylim
            axs[i, 0].set_ylim(minv, maxv*10)

            # Fill between line and O
            axs[i, 0].fill_between(
                time, 0, self.max_vel[i][start:stop], color='lightgrey')

            # Set tick label size
            axs[i, 0].yaxis.set_tick_params(labelsize=9)
            axs[i, 0].xaxis.set_tick_params(labelsize=9)

            # Create grid
            axs[i, 0].grid()

            # Put y axis in log
            axs[i, 0].set_yscale('log')

            # Avoid white space between start and en of axis
            axs[i, 0].margins(x=0)

            # Write subplot title
            axs[i, 0].set_ylabel('Fault {}'.format(i+1), fontsize=11)
            axs[i, 0].yaxis.set_label_position("right")

            # Draw SSE and EQ horizontal lines
            axs[i, 0].axhline(ssel, color='lightgreen',
                              linestyle='--', dashes=(4, 4))
            axs[i, 0].axhline(eql, color='paleturquoise',
                              linestyle='--', dashes=(4, 4))

        # x axis label
        axs[self.nbr_fault - 1, 0].set_xlabel('Time (year)', fontsize=12)

        # Add common ylabel
        fig.text(0.01, 0.5, 'Maximum slip rate ($m.s^{-1}$)',
                 va='center', rotation='vertical', fontsize=12)
        
        # Save if savefig=True
        if savefig:
            fig.savefig(self.path + 'max_vel_evolution.png', dpi=400)
        
    def plot_slip_rate(self, vmask=1e-14, start=0, stop=None, savefig=True):  
        """
        Plot slip rate evolution for all faults.

        Parameters
        ----------
        vmask : float, optional
            Mask velocities below vmask (they will be display in white). The 
            default is 1e-14.
        start : int, optional
            Starting index to plot data. The default is 0.
        stop : int, optional
            Last index to plot data. The default is None (i.e. last data)
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "slip_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # Time in year
        time = self.time[start:stop]/(365.25*24*3600)
        
        # fig, axs = plt.subplots(1, self.nbr_fault + 1, sharey=True, gridspec_kw={'width_ratios': [2, 2, 2]})
        fig, axs = plt.subplots(1, self.nbr_fault + 1, sharey=True)
        
        ################
        # Plot max_vel #
        ################
        
        # Loop over all the fault 
        for i in range(self.nbr_fault):
            axs[0].plot(self.max_vel[i][start:stop], time, label='Fault {}'.format(i+1))
            
        # Put y axis in log
        axs[0].set_xscale('log')
        
        # Invert x-axis
        axs[0].invert_xaxis()
        
        # y axis label
        axs[0].set_ylabel('Time (year)')
        
        # Set y ticks at every power of 10
        axs[0].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))
        
        ##################
        # Plot slip rate #
        ##################
        
        for i in range(self.nbr_fault):
            # Mask (will display in white) velocity values < vmask
            velm = np.ma.masked_where(self.velocity[i] < vmask, self.velocity[i])
            
            # Set the grid
            xx, yy = np.meshgrid(self.ex[i], time)
            
            # Create colormap
            # cmp = nice_colormap(self.velocity[i][0:idxlim].min(), 1e-8, 1e-3, self.velocity[i][0:idxlim].max())
            cmp = nice_colormap()
            
            # Plot
            # norm=colors.LogNorm(vmin=self.velocity[i][0:idxlim].min(), vmax=self.velocity[i][0:idxlim].max())
            cs = axs[i+1].pcolormesh(xx, yy, velm[start:stop], norm=colors.LogNorm(vmin=1e-12, vmax=1), cmap=cmp)
       
            # Remove frame            
            axs[i+1].spines["right"].set_visible(False)
            axs[i+1].spines["left"].set_visible(False)
            axs[i+1].spines["top"].set_visible(False)
            axs[i+1].axes.yaxis.set_visible(False)
            
        
        ############
        # Colorbar #
        ############
        
        # Make some space for the colorbar
        fig.subplots_adjust(right=0.86)

        # Add the colorbar outside of the plot
        box = axs[self.nbr_fault].get_position()
        pad, width = 0.02, 0.04 # set width of the colorbar 
        cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
        fig.colorbar(cs, cax=cax)
        
        # Save figure if savefig=True
        if savefig:
            fig.savefig(self.path + 'slip_rate_evolution.png', dpi=400)

        return
        
    def plot_moment_rate(self):
        # TODO !
        pass