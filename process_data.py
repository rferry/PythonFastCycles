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


class ReadData:
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
        self.max_vel = self.compute_max_vel()

    def read_binary_file(self, file_type):

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

    def compute_max_vel(self):
        # TODO ! Add check if self.velocity is defined
        max_vel = []
        # Loop over the faults
        for i in range(self.nbr_fault):
            max_vel_temp = []
            # Loop over segments
            for j, el in enumerate(self.velocity[i]):
                max_vel_temp.append(np.max(el))
            max_vel.append(max_vel_temp)  # max_vel for one fault

        return max_vel

    def read_moment_rate(self):
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

    def read_GPS_rate(self):
        # TODO !
        pass

    def read_EQcatalog(self):
        # TODO !
        pass

    def plot_max_vel(self, ssel=1e-8, egl=1e-3):
        # Determine yaxis lim --> determine min and max velocity
        minv = 1000  # ridiculous value to always be below
        maxv = -9999  # ridiculous value to always be above
        for i in range(self.nbr_fault):
            mint = np.min(self.max_vel[i])
            maxt = np.max(self.max_vel[i])
            if mint < minv:
                minv = mint
            else:
                pass
            if maxt > maxv:
                maxv = maxt
            else:
                pass

        # Time
        time = self.time / (3600*24*362.25)

        # Create figure
        fig, axs = plt.subplots(self.nbr_fault, 1, squeeze=False)

        for i in range(self.nbr_fault):
            # Plot data
            axs[i, 0].plot(time, self.max_vel[i], color='red')

            # Set ylim
            axs[i, 0].set_ylim(minv, maxv*10)

            # Fill between line and O
            axs[i, 0].fill_between(
                time, 0, self.max_vel[i], color='lightgrey')

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
            axs[i, 0].axhline(egl, color='paleturquoise',
                              linestyle='--', dashes=(4, 4))

        # x axis label
        axs[self.nbr_fault - 1, 0].set_xlabel('Time (year)', fontsize=12)

        # Add common ylabel
        fig.text(0.01, 0.5, 'Maximum slip rate ($m.s^{-1}$)',
                 va='center', rotation='vertical', fontsize=12)
