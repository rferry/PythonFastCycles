"""
Class to read data. 

First version by R. Ferry on January 2021.
"""
# External imports
import numpy as np
import os
from scipy.io import FortranFile
import scipy.integrate 
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.ticker as ticker

# Internal import
from tools import nice_colormap

class ReadData:
    def __init__(self, path, tLnuc='Lab'):
        """
        Initialisation of the class. Call functions to read data.

        Parameters
        ----------
        path : string
            Path to simulation directory. 
        tLnuc : string, optional
            Type of Lnuc to use. Can be:
                * Lab: Lnuc = (mu * Dc) / (sigmaN * (b - a))
                * Linf: Lnuc = (2 * mu * Dc) / (pi * sigmaN * (1 - a/b)**2)
            The default is 'Lab'.    
            
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
        
        # Read GPS.in
        self.read_GPS()
        
        # Read GPSRate.out
        self.read_GPS_rate()

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
        
        # Read EQ catalog
        self.read_EQcatalog()
        
        # Compute L, Lnuc, Linf and Lb
        self.L = []
        self.Lab = []  # L_a-b
        self.Linf = []
        self.Lb = []
        for i in range(self.nbr_fault):
            self.L.append(self.elements[i] * self.ds[i])
            self.Lab.append(-(self.mu * self.Dc[i]) / (self.sigmaN[i] * 
                                                    (self.b[i] - self.a[i])))
            self.Linf.append(-(2 * self.mu * self.Dc[i]) / 
                             (np.pi * self.sigmaN[i] * self.b[i] * 
                              (1 - self.a[i] / self.b[i])**2))
            self.Lb.append(-(self.mu * self.Dc[i]) / (self.sigmaN[i] * 
                                                      self.b[i]))
        
        # Choose the correct Lnuc
        if tLnuc == 'Lab':
            self.Lnuc = self.Lab
        elif tLnuc == 'Linf':
            self.Lnuc = self.Linf
        else:
            print("Invalide tLnuc !")
        

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
        # Verbose
        print('Reading {} files...'.format(file_type))
        
        files = []  # list of files of the right type
        # Loops over the files and select those with 'file_type' in their name
        for file_name in self.files_in_dir:
            if file_type in file_name:
                files.append(file_name)

        # Sorts selected files in ascending order
        files.sort()
        
        # Dealing with '*******' files
        if files[0] == file_type + '*******':
            # Move the first element to the last position 
            files.append(files[0])
            del files [0]

        data = []  # list to store the data
        data_temp = []  # list to store temporary data for velocity
        # Loops over selected files and read them
        for file in files:
            f = FortranFile(self.path + file, 'r')
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

        # Verbose
        print("...{} files read. \n".format(file_type))
        
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
        # Initialisation
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
        
        tMrate = []  # to store moment rate time
        netMrate = []  # to store net moment rate
        Mrate = []  # to store moment rate for each fault
        max_vels = []  # to store max slip rate for each fault
        # stores Mrate data temporary
        Mrate_temp = np.ones((len(content)-2, self.nbr_fault))
        max_vels_temp = np.ones((len(content)-2, self.nbr_fault))
        for i, line in enumerate(content[2:]):  # [2:] to skip headers
            tMrate.append(float(line.split()[1]))
            netMrate.append(float(line.split()[2]))
            for j in range(self.nbr_fault):  # Loop over faults
                Mrate_temp[i, j] = line.split()[3+j]
                max_vels_temp[i, j] = line.split()[3+self.nbr_fault+j]
        for i in range(self.nbr_fault):
            Mrate.append(Mrate_temp[:, i])
            max_vels.append(list(max_vels_temp[:, i]))

        # Stores results
        self.netMrate = netMrate
        self.Mrate = Mrate
        self.tMrate = tMrate
        self.max_vels = max_vels

    def read_config(self):
        """
        Read "config.in" file to get simulation parameters.

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
        V0 = []
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
            elif line.startswith('V0_val'):
                V0.append(float((line.split()[2]).split(',')[2])) 
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
            elif line.startswith('sigma12_dot_inf') or \
                    line.startswith('sigma21_dot_inf'):
                sigma_dot[0, 1] = float(line.split()[2])  
                sigma_dot[1, 0] = float(line.split()[2])
            elif line.startswith('sigma31_dot_inf') or \
                    line.startswith('sigma13_dot_inf'):
                sigma_dot[0, 2] = float(line.split()[2])  
                sigma_dot[2, 0] = float(line.split()[2])
            elif line.startswith('sigma32_dot_inf') or \
                    line.startswith('sigma23_dot_inf'):
                sigma_dot[2, 1] = float(line.split()[2])  
                sigma_dot[1, 2] = float(line.split()[2])  
        
        # Store values 
        self.a = a
        self.b = b
        self.Dc = Dc
        self.V0 = V0
        self.sigmaN = sigmaN
        self.ds = ds
        self.sigma_dot = sigma_dot 
    
    def read_tides(self):
        """
        Read "tides.in" to get amplitude, period and phase of the waves.

        Returns
        -------
        None.

        """
        # Open tides.in and read content
        with open(self.path + 'tides.in', 'r') as file:
            content = file.readlines()  
        
        # Initialisation
        Tampli = []
        Tperiod = []
        Tphase = []
        
        # Loop over all waves
        for i, line in enumerate(content[:-1]):
            Tampli.append(float(line.split()[0]))
            Tperiod.append(float(line.split()[1]))
            Tphase.append(float(line.split()[2])) 
        
        # Store data
        self.Tampli = Tampli
        self.Tperiod = Tperiod
        self.Tphase = Tphase 
    
    def read_GPS(self):
        """
        Read "GPS.in" file to get GPS stations coordinates.

        Returns
        -------
        None.

        """
        # Open config.in and read content
        with open(self.path + 'GPS.in', 'r') as file:
            content = file.readlines()
        
        # Initialisation
        GPSx = []
        GPSy = []
        
        # Loop over all stations 
        for i, line in enumerate(content[:-1]):
            GPSx.append(float(line.split()[0]))
            GPSy.append(float(line.split()[1]))
        
        # Store data 
        self.GPSx = GPSx
        self.GPSy = GPSy
        
    def read_GPS_rate(self):
        """
        Read "GPSRate.out" file to get GPS rate time and GPS rate for all GPS
        stations.

        Returns
        -------
        None.

        """
        # Open and read GPSRate.out
        with open(self.path + 'GPSRate.out', 'r') as file:
            content = file.readlines()
            
        tGPSrate = [] # to store GPS rate time
        GPSrate = [] # to store GPS rate for each GPS station
        nGPS = len(self.GPSx) # number of GPS stations
        # stores GPSrate data temporary
        GPSrate_temp = np.ones((len(content)-2, nGPS))
        for i, line in enumerate(content[2:]):  # [2:] to skip headers
            tGPSrate.append(float(line.split()[1]))
            for j in range(nGPS):  # Loop over GPS stations 
                GPSrate_temp[i, j] = line.split()[2+j]
        for i in range(nGPS):
            GPSrate.append(GPSrate_temp[:, i])

        # Stores results
        self.tGPSrate = tGPSrate
        self.GPSrate = GPSrate
        
    def read_EQcatalog(self):
        """
        Read "EQcatalog" file to build a general earthquake catalog 
        (EQgeneral_catalog) and a catalog for each fault (EQcatalog).

        Returns
        -------
        None.

        """
        # Open and read Eqcatalog
        with open(self.path + 'EQcatalog', 'r') as file:
            content = file.readlines()
            
        # Initialisation
        EQcatalog = []  # Catalog per fault
        for i in range(self.nbr_fault):
            EQcatalog.append({})
        EQgeneral_catalog = {} # common catalog for all fault
        string = ['Fault', 'Type', 'Time Beg', 'Nucleation dur', 'Rupture dur', 
                  'After slip dur', 'Time Index Beg', 'Time Index Beg Dyn', 
                  'Time Index End Dyn', 'Time Index End', 'Location X', 
                  'Location Y']
        for i, el in enumerate(string):
            EQgeneral_catalog[el] = []
            for j in range(self.nbr_fault):
                EQcatalog[j][el] = []
        
        for i, line in enumerate(content[2:]):  # [2:] to skip headers
            # Building general catalog
            EQgeneral_catalog['Fault'].append(int(line.split()[1]))
            EQgeneral_catalog['Type'].append(int(line.split()[2]))
            EQgeneral_catalog['Time Beg'].append(float(line.split()[3]))
            EQgeneral_catalog['Nucleation dur'].append(float(line.split()[4]))
            EQgeneral_catalog['Rupture dur'].append(float(line.split()[5]))
            EQgeneral_catalog['After slip dur'].append(float(line.split()[6]))  
            EQgeneral_catalog['Time Index Beg'].append(int(line.split()[7])) 
            EQgeneral_catalog['Time Index Beg Dyn'].append(int(line.split()[8]))  
            EQgeneral_catalog['Time Index End Dyn'].append(int(line.split()[9]))   
            EQgeneral_catalog['Time Index End'].append(int(line.split()[10]))
            EQgeneral_catalog['Location X'].append(float(line.split()[11]))
            EQgeneral_catalog['Location Y'].append(float(line.split()[12])) 
            
            # Add the event to the fault's catalog
            fault = int(line.split()[1]) - 1
            EQcatalog[fault]['Type'].append(int(line.split()[2]))
            EQcatalog[fault]['Time Beg'].append(float(line.split()[3]))
            EQcatalog[fault]['Nucleation dur'].append(float(line.split()[4]))
            EQcatalog[fault]['Rupture dur'].append(float(line.split()[5]))
            EQcatalog[fault]['After slip dur'].append(float(line.split()[6]))  
            EQcatalog[fault]['Time Index Beg'].append(int(line.split()[7])) 
            EQcatalog[fault]['Time Index Beg Dyn'].append(int(line.split()[8]))  
            EQcatalog[fault]['Time Index End Dyn'].append(int(line.split()[9]))   
            EQcatalog[fault]['Time Index End'].append(int(line.split()[10]))
            EQcatalog[fault]['Location X'].append(float(line.split()[11]))
            EQcatalog[fault]['Location Y'].append(float(line.split()[12])) 
             
        # Store data
        self.EQgeneral_catalog = EQgeneral_catalog
        self.EQcatalog = EQcatalog
        
        return

    def compute_EQcatalog(self, eql=1e-3, ssel=1e-8):
        # TODO !
        pass  
    
    def compute_event_moment(self, depth):
        # d = np.zeros(len(self.EQgeneral_catalog['Fault']))
        M0 = np.zeros(len(self.EQgeneral_catalog['Fault']))
        # Loop over all events
        for i, fault in enumerate(self.EQgeneral_catalog['Fault']):
            # Convert string to float 
            fault = int(fault) - 1
            # Index of the beginning 
            beg = self.EQgeneral_catalog['Time Index Beg'][i]
            diff = np.abs(self.time - self.tMrate[beg])
            ibeg = np.where(diff == diff.min())[0][0]
            # Index of the end 
            end = self.EQgeneral_catalog['Time Index End'][i]
            diff = np.abs(self.time - self.tMrate[end])
            iend = np.where(diff == diff.min())[0][0]
            
            # Mean velocities (mean over all duration of the event for each 
            # element of the fault)

            mean_vels = np.mean(self.velocity[fault][ibeg:iend, :], axis=0)
            # Mean displacement 

            d = mean_vels * (self.time[iend] - self.time[ibeg])
            
            threshold = 0.001 * np.max(d)
            print(threshold)
            Lrup = sum(d >= threshold) * self.ds[fault]

            d_mean = d.mean() 
            M0[i] = self.mu * Lrup * depth * d_mean
            
        Mw = (2/3)*np.log10(M0) - 6.01   
            
        return Mw
            
        
    def compute_tractionT(self, f0=0.6):
        """
        Compute shear traction on the faults.
        tau = tractionN * [f0 + aln(V/V0) + bln((theta*V0)/Dc)]

        Parameters
        ----------
        f0 : float, optional
            Steady-state friction at V=V0. The default is 0.6.

        Returns
        -------
        None.

        """
        # Initiaisation
        tractionT = []
        
        # Loop over fault
        for i in range(self.nbr_fault):
            # f = f0 + aln(V/V0) + bln((theta*V0)/Dc)
            f = f0 + self.a[i] * np.log(self.velocity[i]/self.V0[i]) + \
                self.b * np.log((self.theta[i] * self.V0[i]) / self.Dc[i])
            # tau = sigma_N * f    
            tractionT.append(self.tractionN[i] * f)
        
        # Store data 
        self.tractionT = tractionT   
    
    def plot_max_vel(self, eql=1e-3, ssel=1e-8, start=0, stop=None, 
                     vel_tot=False, plot_type='each', savefig=True):
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
        vel_tot : bool, optional
            If True, will use the max_vel from MomentRate.out for each time 
            step. If False will use the max_vel computed from velocity data. 
            The default is False.
        plot_type : string, optional
             Type of plot. Can be:
                * 'each': a subplot for each fault
                * 'all': all faults on the same plot
            The default is 'each'.           
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "slip_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # Set which data to use
        if vel_tot:
            max_vel = self.max_vels
            time = np.array(self.tMrate[start:stop])
        else:
            max_vel = self.max_vel
            time = self.time[start:stop]
        
        # Determine yaxis lim --> determine min and max velocity
        minv = 1000  # ridiculous value to always be below
        maxv = -9999  # ridiculous value to always be above
        for i in range(self.nbr_fault):
            mint = np.min(max_vel[i][start:stop])
            maxt = np.max(max_vel[i][start:stop])
            if mint < minv:
                minv = mint
            else:
                pass
            if maxt > maxv:
                maxv = maxt
            else:
                pass
        
        # Convert time to years 
        time = time / (3600*24*362.25)

        # Create figure        
        if plot_type=='each':
            fig, axs = plt.subplots(self.nbr_fault, 1, squeeze=False)
        elif plot_type=='all':
            fig, axs = plt.subplots(1, 1, squeeze=False)

        for i in range(self.nbr_fault):
            if plot_type=='each':
                # Plot data
                axs[i, 0].plot(time, max_vel[i][start:stop], color='red')
    
                # Fill between line and O
                axs[i, 0].fill_between(
                    time, 0, max_vel[i][start:stop], color='lightgrey')
    
                # Write subplot title if there is more than one fault
                if self.nbr_fault > 1 :
                    axs[i, 0].set_ylabel('Fault {}'.format(i+1), fontsize=11)
                    axs[i, 0].yaxis.set_label_position("right")
                
            elif plot_type=='all':
                axs[0, 0].plot(time, max_vel[i][start:stop], \
                               label='Fault {}'.format(i+1)) 
                    
        for ax in fig.axes:
            # Put y axis in log
            ax.set_yscale('log')
            
            # Set ylim
            ax.set_ylim(minv, maxv + 0.1*(maxv-minv))
            
            # Avoid white space between start and en of axis
            ax.margins(x=0)
            
            # Create grid
            ax.grid()
            
            # Set tick label size
            ax.yaxis.set_tick_params(labelsize=9)
            ax.xaxis.set_tick_params(labelsize=9)
            
            # Draw SSE and EQ horizontal lines
            xmin, xmax = ax.get_xlim()
            ax.axhline(ssel, xmin, xmax, color='lightgreen', \
                     linestyle='--', dashes=(4, 4))
            ax.axhline(eql, xmin, xmax, color='paleturquoise', \
                     linestyle='--', dashes=(4, 4))
            
            # Plot legend if there is a single plot AND more than one fault
            if plot_type=="all" and self.nbr_fault>1:    
                ax.legend(loc='upper left')
        
        # x axis label
        try:  # if "each"
            axs[self.nbr_fault - 1, 0].set_xlabel('Time (year)', fontsize=12)
        except:  # if "all"
            axs[0, 0].set_xlabel('Time (year)', fontsize=12)      
        
        # Add common ylabel
        fig.text(0.01, 0.5, 'Maximum slip rate ($m.s^{-1}$)',
                 va='center', rotation='vertical', fontsize=12)
        
        # Save if savefig=True
        if savefig:
            fig.savefig(self.path + 
                        'max_vel_evolution_{}.png'.format(plot_type), dpi=400)
        
    def plot_slip_rate(self, vmask=1e-14, start=0, stop=None, sharey=True, 
                       savefig=True):  
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
        sharey : bool, optional
            If True, time will be linear. If False, time will depends on 
            timestep. The default is True.
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "slip_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # Time in year
        time = self.time[start:stop]/(365.25*24*3600)
        
        # TODO ? Ajouter gridspec_kw={'width_ratios': [2, 2, 2]}
        fig, axs = plt.subplots(1, self.nbr_fault + 1, sharey=sharey)
        
        ################
        # Plot max_vel #
        ################
        
        # Loop over all the fault 
        for i in range(self.nbr_fault):
            axs[0].plot(self.max_vel[i][start:stop], time, 
                        label='Fault {}'.format(i+1))
            
        # Put y axis in log
        axs[0].set_xscale('log')
        
        # Invert x-axis
        axs[0].invert_xaxis()
        
        # x axis label
        axs[0].set_xlabel('$V \ (m.s^{-1})$', fontsize=12)
        
        # y axis label
        axs[0].set_ylabel('Time (year)', fontsize=12)
        
        # Set y ticks at every power of 10
        axs[0].xaxis.set_major_locator(ticker.LogLocator(base=10.0, 
                                                         numticks=5))
        
        # Remove white at begenning and end of y axis 
        axs[0].margins(y=0)
        
        ##################
        # Plot slip rate #
        ##################
        
        for i in range(self.nbr_fault):
            # Mask (will display in white) velocity values < vmask
            velm = np.ma.masked_where(self.velocity[i] < vmask, 
                                      self.velocity[i])
            
            # Set the grid
            if sharey:
                y = time
            else:
                y = np.arange(len(time))
            xx, yy = np.meshgrid(self.ex[i], y)
            
            # Create colormap
            cmp = nice_colormap()
            
            # Plot
            cs = axs[i+1].pcolormesh(xx, yy, velm[start:stop], 
                                     norm=colors.LogNorm(vmin=1e-12, vmax=1), 
                                         shading='nearest', cmap=cmp)
       
            # Remove frame            
            axs[i+1].spines["right"].set_visible(False)
            axs[i+1].spines["left"].set_visible(False)
            axs[i+1].spines["top"].set_visible(False)
            axs[i+1].axes.yaxis.set_visible(False)
            
            # x axis coordinate in X/Lnuc
            # TODO ! Needs to be improved 
            nticks = 5 # Number of ticks
            L_over_Lnuc = round(self.L[i] / self.Lnuc[i], 1)
            labels = np.linspace(-L_over_Lnuc/2, L_over_Lnuc/2, nticks)
            labels = np.round(labels, 2)
            left, right = axs[i+1].get_xlim()
            tick_loc = np.linspace(left, right, nticks)
            axs[i+1].set_xticks(tick_loc)
            axs[i+1].set_xticklabels(labels)
            
            # x axis label
            axs[i+1].set_xlabel('$X/L_{nuc}$', fontsize=12)
            
        
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

        return fig
        
    def plot_moment_rate(self, start=0, stop=None, savefig=True):
        """
        Plot moment rate evolution for all fault and net moment rate.

        Parameters
        ----------
        start : int, optional
            Starting index to plot data. The default is 0.
        stop : int, optional
            Last index to plot data. The default is None (i.e. last data)
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "moment_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # Determine yaxis lim --> determine min and max velocity
        minv = 1000  # ridiculous value to always be below
        maxv = -9999  # ridiculous value to always be above
        for i in range(self.nbr_fault):
            mint = np.min(self.Mrate[i][start:stop])
            maxt = np.max(self.Mrate[i][start:stop])
            if mint < minv:
                minv = mint
            else:
                pass
            if maxt > maxv:
                maxv = maxt
            else:
                pass

        # Time
        time = np.array(self.tMrate[start:stop]) / (3600*24*365.25)

        # Create figure
        fig, axs = plt.subplots(self.nbr_fault + 1, 1, squeeze=False)

        # Plot moment rate for each fault
        for i in range(self.nbr_fault):
            # Plot data
            axs[i, 0].plot(time, self.Mrate[i][start:stop], color='red')

            # Set ylim
            axs[i, 0].set_ylim(minv, maxv*10)

            # Fill between line and O
            axs[i, 0].fill_between(
                time, 0, self.Mrate[i][start:stop], color='lightgrey')

            # Write subplot title for each fault
            axs[i, 0].set_ylabel('Fault {}'.format(i+1), fontsize=11)
            axs[i, 0].yaxis.set_label_position("right") 
                
        # Plot net moment rate 
        axs[self.nbr_fault, 0].plot(time, self.netMrate[start:stop], color='k')
        
        # Fill between line and O
        axs[self.nbr_fault, 0].fill_between(time, 0, \
                                self.netMrate[start:stop], color='lightgrey')
        
        # Write subplot title for net moment rate
        axs[self.nbr_fault, 0].set_ylabel('Net Mo Rate', fontsize=11)
        axs[self.nbr_fault, 0].yaxis.set_label_position("right")
        
        for ax in fig.axes:
            # Put y axis in log
            ax.set_yscale('log')
            
            # Avoid white space between start and en of axis
            ax.margins(x=0)
            
            # Create grid
            ax.grid()
            
            # Set tick label size
            ax.yaxis.set_tick_params(labelsize=9)
            ax.xaxis.set_tick_params(labelsize=9)
        
        # x axis label
        axs[self.nbr_fault, 0].set_xlabel('Time (year)', fontsize=12)

        # Add common ylabel
        fig.text(0.01, 0.5, 'Moment Rate ($Nm.s^{-1})$',
                 va='center', rotation='vertical', fontsize=12)
        
        # TODO ! Compute a better EQcatalog to have the right locations
        # ymin, ymax = axs[self.nbr_fault, 0].get_ylim()
        
        # ymax = [ymax] * len(self.EQgeneral_catalog['Time Index Beg'])
        # tEQ = np.array(self.EQgeneral_catalog['Time Beg']) / (365.25*24*3600)
        # axs[self.nbr_fault, 0].scatter(tEQ, ymax)
        
        # Save if savefig=True
        if savefig:
            fig.savefig(self.path + 'moment_rate_evolution.png', dpi=400)
        pass
    
    def plot_geometry(self, scale='Lnuc', savefig=True):
        """
        Plot the fault system geometry.

        Parameters
        ----------
        scale : str, optional
            Scale of axis. Can be 'Lnuc' (normalised by Lnuc) or 'X'. The 
            default is 'Lnuc'.
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "slip_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # TODO ! Add GPS station plot 
        # Computing figure limits 
        maxi = -np.inf
        mini = np.inf
        for i in range(self.nbr_fault):
            maxtemp = np.max(self.ey[i])
            mintemp = np.min(self.ey[i])
            if maxtemp > maxi:
                maxi = maxtemp
            if mintemp < mini:
                mini = mintemp
        if scale == 'Lnuc':        
            mini = mini / self.Lnuc[0]
            maxi = maxi / self.Lnuc[0]
        Ly = maxi-mini # extend in y direction of the fault system
        
        # Initialise figure
        fig, ax = plt.subplots(1, 1)
        
        # Plot each fault 
        for i in range(self.nbr_fault):
            if scale == 'Lnuc':
                x = [el/self.Lnuc[0] for el in self.ex[i]]
                y = [el/self.Lnuc[0] for el in self.ey[i]]
            else:
                x = self.ex[i]
                y = self.ey[i]
            
            ax.plot(x, y, 'b')
        
            # Add fault "label"
            xtext = 0.5 *(np.max(x) + np.min(x))
            ytext = 0.5 *(np.max(y) + np.min(y))
            ax.text(xtext, ytext, 'F{}'.format(i+1), bbox=dict(fc='yellow', 
                                   ec='none', pad=1), ha='center', va='center')
            
        # Set axis limits and aspect 
        ax.set_ylim(mini - Ly*0.1, maxi + Ly*0.1)
        ax.set_aspect('equal')
        
        # Axis labels
        if scale == 'Lnuc':
            ax.set_xlabel('$X/L_{nuc}$', fontsize=12)
            ax.set_ylabel('$Y/L_{nuc}$', fontsize=12)
        else:
            ax.set_xlabel('$X$', fontsize=12)
            ax.set_ylabel('$Y$', fontsize=12)
            
        # Change tick size    
        ax.tick_params(axis='both', which='major', labelsize=10)    
        
        if savefig:
            fig.savefig(self.path + 'geometry_out_' + scale + '.png', dpi=400, 
                        bbox_inches='tight')
            
    def plot_GPS_rate(self, start=0, stop=None, savefig=True):
        """
        Plot GPS rate for all stations.

        Parameters
        ----------
        start : int, optional
            Starting index to plot data. The default is 0.
        stop : int, optional
            Last index to plot data. The default is None (i.e. last data)
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "GPS_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # TODO ! See why we have weird negative values
        # Determine yaxis lim --> determine min and max velocity
        minv = 1000  # ridiculous value to always be below
        maxv = -9999  # ridiculous value to always be above
        for i in range(len(self.GPSrate)):
            mint = np.min(self.GPSrate[i][start:stop])
            maxt = np.max(self.GPSrate[i][start:stop])
            if mint < minv:
                minv = mint
            else:
                pass
            if maxt > maxv:
                maxv = maxt
            else:
                pass

        # Time
        time = np.array(self.tGPSrate[start:stop]) / (3600*24*365.25)

        # Create figure
        fig, axs = plt.subplots(len(self.GPSrate), 1, squeeze=False)

        # Plot GPS rate for each fault
        for i in range(len(self.GPSrate)):
            # Plot data
            axs[i, 0].plot(time, self.GPSrate[i][start:stop], color='red')

            # Set ylim
            # axs[i, 0].set_ylim(minv, maxv*10)

            # Fill between line and O
            axs[i, 0].fill_between(
                time, 0, self.GPSrate[i][start:stop], color='lightgrey')

            # Write subplot title for each station
            axs[i, 0].set_ylabel('Station {}'.format(i+1), fontsize=11)
            axs[i, 0].yaxis.set_label_position("right") 
                
        
        for ax in fig.axes:
            # Put y axis in log
            ax.set_yscale('log')
            
            # Avoid white space between start and en of axis
            ax.margins(x=0)
            
            # Create grid
            ax.grid()
            
            # Set tick label size
            ax.yaxis.set_tick_params(labelsize=9)
            ax.xaxis.set_tick_params(labelsize=9)
        
        # x axis label
        axs[len(self.GPSrate)-1, 0].set_xlabel('Time (year)', fontsize=12)

        # Add common ylabel
        fig.text(0.01, 0.5, 'GPS Rate ($m.s^{-1})$',
                 va='center', rotation='vertical', fontsize=12)
        
        # Save if savefig=True
        if savefig:
            fig.savefig(self.path + 'GPS_rate_evolution.png', dpi=400, \
                        bbox_inches='tight')
            
    def plot_GPS_disp(self, start=0, stop=None, plot_type='all', savefig=True):
        """
        Plot GPS cumulative displacement over time.

        Parameters
        ----------
        start : int, optional
            Starting index to plot data. The default is 0.
        stop : int, optional
            Last index to plot data. The default is None (i.e. last data)
        plot_type : string, optional
            Type of plot. Can be:
                * 'all': all stations on the same plot
                * 'each': a subplot for each stations
            The default is 'all'.
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "GPS_displacement.png". The default is True.

        Returns
        -------
        None.

        """
        ########################
        # Compute displacement #
        ########################
        
        disp = []  # to store displacement for each GPS station
        
        # Loop over GPS stations
        for i in range(len(self.GPSrate)):
            disp_temp = np.zeros(len(self.GPSrate[i]))
            disp_temp[:] = scipy.integrate.cumtrapz(self.GPSrate[i], \
                                                    self.tGPSrate, initial=0)
            disp.append(disp_temp)
        
        # Store data
        self.disp = disp
        
        ########
        # Plot #
        ########
        
        # Compute time in year
        time = np.array(self.tGPSrate[start:stop]) / (3600*24*365.25)
        
        # Create figure 
        if plot_type=='each':
            fig, axs = plt.subplots(len(self.GPSrate), 1, squeeze=False)
        elif plot_type=='all':
            fig, axs = plt.subplots(1, 1, squeeze=False)
        
        # Plot each GPS stations
        for i in range(len(self.GPSrate)):
            if plot_type=='each':
                axs[i, 0].plot(time, self.disp[i][start:stop], color='k')
                axs[i, 0].set_ylabel('Station {}'.format(i+1), fontsize=11)
                axs[i, 0].yaxis.set_label_position("right")
            elif plot_type=='all':
                axs[0, 0].plot(time, self.disp[i][start:stop], \
                               label='Station {}'.format(i+1))                

        # EQ time to plot lines for EQ
        tEQ = np.array(self.EQgeneral_catalog['Time Beg']) / (365.25*24*3600)
        # Colors of EQ
        colors=[]
        for i, el in enumerate(self.EQgeneral_catalog['Type']):
            if el==1:
                colors.append('pink') # EQ are in red
            elif el==2:
                colors.append('cyan') # Slow slip are in blue
        
        # Loop over axes 
        for ax in fig.axes:
            # Avoid white space between start and en of axis
            ax.margins(x=0)
            ax.margins(y=0)
            
            # Set tick label size
            ax.yaxis.set_tick_params(labelsize=9)
            ax.xaxis.set_tick_params(labelsize=9)
            
            # Plot lines for EQ
            ymin, ymax = ax.get_ylim()
            ax.vlines(tEQ, ymin, ymax, colors=colors, linestyle='--')
        
        # x axis label
        try:  # if "each"
            axs[len(self.GPSrate)-1, 0].set_xlabel('Time (year)', fontsize=12)
        except:  # if "all"
            axs[0, 0].set_xlabel('Time (year)', fontsize=12)  
        
        # Add common ylabel
        fig.text(0.01, 0.5, 'Displacement $(m)$', va='center', \
                 rotation='vertical', fontsize=12)           

        # Plot legend if there is a single plot AND more than one station
        if plot_type=="all" and len(self.GPSrate) != 1:
            fig.legend(loc='upper left')
            
        # Save if savefig=True
        if savefig:
            fig.savefig(self.path + 'GPS_displacements.png', dpi=400, \
                        bbox_inches='tight')

        return fig, axs