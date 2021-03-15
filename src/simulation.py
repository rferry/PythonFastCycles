"""
Class to create all files for a simulation.

First version by R. Ferry on January 2021.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

class Simulation:
    def __init__(self, path, simuname, mu, a, b, fric_law='RateStateAgeing_R',
                 frac_mode='ModeIII', sigma_N=-1e8, Dc=1e-3, tLnuc='Lab'):
        """
        Initialise Simulation class. Compute Lnuc. 

        Parameters
        ----------
        path : string
            Path to problems directory. 
        simuname : string
            Name of the simulation.
        mu : float
            Shear stress modulus in Pa.
        a : float
            Frictional empirical parameter.
        b : float
            Frictional empirical parameter.
        fric_law : str, optional
            Friction law used by FastCycles. Can be :
                * RateStateAgeing_R
                * RateStateAgeing
                * RateStateSlip_R
                * RateStateSlip
            '_R' means Regularized version of the Friction Law. The default is 
            'RateStateAgeing_R'.
        frac_mode : str, optional
            Fracture mode. Can be'ModeII' or 'ModeIII'. The default is 
            'ModeIII'.
        sigma_N : float, optional
            Normal stress. The default is -1e8.
        Dc : float, optional
            Critical distance Dc. The default is 1e-3.
        tLnuc : string, optional
            Type of Lnuc to use. Can be:
                * Lab: Lnuc = (mu * Dc) / (sigmaN * (b - a))
                * Linf: Lnuc = (2 * mu * Dc) / (pi * sigmaN * (1 - a/b)**2)
            The default is 'Lab'.   

        Returns
        -------
        None.

        """
        # Store parameters 
        self.path = path + simuname + '/'
        self.simuname = simuname
        self.mu = mu
        self.fric_law = fric_law
        self.sigma_N = sigma_N
        self.Dc = Dc
        self.frac_mode = frac_mode
        self.a = a
        self.b = b

        # Compute Lab and Linf
        self.Lab = - (self.mu * self.Dc) / (self.sigma_N * (self.b - self.a))
        self.Linf = -(2 * self.mu * self.Dc) / (np.pi * self.sigma_N * self.b 
                                                * (1 - self.a / self.b)**2)
        
        # Choose Lnuc
        if tLnuc == 'Lab':
            self.Lnuc = self.Lab
        elif tLnuc == 'Linf':
            self.Lnuc = self.Linf
        else:
            print("Invalide tLnuc !")
        

    def create_all_files(self, sigma_dot, geom_type, L_over_Lnuc=2, show=False,
                         lengths=None, angles=None, xs=None, ys=None, 
                         D_over_Lnuc=0.1, overlap=0.5, GPSx=[10], GPSy=[10], 
                         Tampli=[0.0], Tperiod=[1.0], Tphase=[0.0], 
                         Vval_x1='default', Vval_x2='default', \
                         Vval_pourc=0.001, stop_crit=1, max_it=10000, 
                         final_time=10000000, tol_solver=1.00e-8, nf=False):
        """
        Create all files for a simulation, i.e. "config.in", "geometry.in", 
        "tides.in" and "GPS.in".

        Parameters
        ----------
        sigma_dot : array 3*3 of float
            Loading rate symmetrical matrix.
        geom_type : string
            Geometry of the fault system. Can be 
                * "1fault" for a single fault of length L defines with 
                  L_over_Lnuc 
                * "2faults_overlapping" for a configuration like Romanet et al. 
                  (2018).
                  Space between two faults can be specified with D_over_Lnuc
                  and the overlap with overlap       
                                 L/Lnuc * overlap
                         <--->   
                         ==========
                          | D/Lnuc  
                    ==========
                    <-------->
                      L/Lnuc
                * "multiple" for a multiple faults configuration. Faults are 
                  defined with a length (lengths), an angle (angles) and the 
                  distance in x (xs) and y (ys) of one edge from the first 
                  fault.
                  
                      /
                     /        
                    /
                   /
                  + <-- Point defining the fault (given as x and y distances 
                        from the origin). Here the angle is positive and is 
                        ~60°.
                 
        L_over_Lnuc : real
            Length of the fault express as the ratio L/Lnuc. Default is 2.
        show : bool
            If True, plot the geometry. Default is False.
        lengths : list of float
            Lengths of the faults normalised by Lnuc (L/Lnuc).
        angles : list of float
            Angles of the faults.
        xs : list of float
            Distance in x normalised by Lnuc of the points defining the faults 
            and the origin.
        ys : list of float
            Distance in y normalised by Lnuc of the points defining the faults 
            and the origin.    
        D_over_Lnuc : real
            Distance between the two fault express as the ratio D/Lnuc. Default
            is 0.1.
        overlap : real
            Portion of the fault overlapping. 0 < overlap < 1. Default is 0.5.
            overlap = 0:
                      =======
               ======= 
            overlap = 1:
                =======
                =======
        GPSx : list
            First coordinate of the GPS station. Default is [10].
        GPSy : list
            Second cooridnate of the GPS station. Default is [10].
        Tampli : array of float, optional
            Amplitude of the waves. The default is [0.0]
        Tperiod : array of float, optional
            Period of the waves. The default is [0.0].
        Tphase : array of float, optional
            Phase of the waves. The default is [0.0].
        Vval_x1 and Vval_x2: float, optional
            The initial velocity perturbation is imposed between Vval_x1 and 
            Vval_x2. The default is between 0.4*Lnuc and 0.6*Lnuc.
        Vval_pourc : float, optional
            Amplitude of the initial perturbation as a percentage of V0. The 
            default is 0.001.
        stop_crit : int, optional
            Criteria to stop the simulation. Can be 0, 1 or 2.
                * 0 : simulation will stop after the first event
                * 1 : simulation will stop after max_it iterations
                * 2 : simulation will stop at final_time
            The default is 1.
        max_it : int, optional
            Number of iteration for the simulation. It will only be taken into
            consideration if stop_crit = 1. The default is 10000.
        final_time : float, optional
            Time in seconds at which the simulation will stop. It will only be 
            taken into consideration if stop_crit = 2. The default is 10000000.
        nf : int, optional
            Number of fault. If not specified, the number of fault will be 
            read in the "geometry.in" file. The default is not specified.

        Returns
        -------
        None.

        """
        # Creates the directory for the simulation
        try:
            os.mkdir(self.path)
        except:
            pass

        self.L_over_Lnuc = L_over_Lnuc
        self.create_tides_file(Tampli, Tperiod, Tphase)
        self.create_GPS_file(GPSx, GPSy)

        if geom_type == '1fault':
            x1 = 0
            x2 = self.Lnuc * self.L_over_Lnuc
            self.create_geom_1fault(x1, x2, show=show)
        elif geom_type == '2faults_overlapping':
            self.create_geom_2_faults_overlapping(
                D_over_Lnuc, L_over_Lnuc, overlap, show=show)
        elif geom_type == "multiple":
            self.create_geom_multiple_faults(lengths, angles, xs, ys, 
                                             show=show)
        else:
            raise Exception('geom_type does not exist.')

        self.create_config_file(sigma_dot, Vval_x1, Vval_x2, Vval_pourc, 
                                stop_crit, max_it, final_time, tol_solver, nf)

    def create_GPS_file(self, GPSx=[10], GPSy=[10]):
        """
        Creates GPS.in at the location 'path'.

        Parameters
        ----------
        GPSx : list
            First coordinate of the GPS station. Default is [10].
        GPSy : list
            Second cooridnate of the GPS station. Default is [10].

        Returns
        -------
        None.

        """
        # Store GPS coordinates
        self.GPSx = GPSx
        self.GPSy = GPSy

        if len(self.GPSx) != len(self.GPSy):
            print('x1 and x2 are not of the same size !')
        else:
            content = []
            for i, x in enumerate(self.GPSx):
                content.append('{} {} \n'.format(x, self.GPSy[i]))
            content.append('/')

            with open(self.path + 'GPS.in', 'w') as file:
                for line in content:
                    file.write(line)

        return

    def create_tides_file(self, Tampli=[0.0], Tperiod=[1.0], Tphase=[0.0]):
        """
        Create "tides.in" file in the simulation directory. Default is no 
        tides.

        Parameters
        ----------
        Tampli : array of float, optional
            Amplitude of the waves. The default is [0.0]
        Tperiod : array of float, optional
            Period of the waves. The default is [0.0].
        Tphase : array of float, optional
            Phase of the waves. The default is [0.0].

        Returns
        -------
        None.

        """
        # Store tides informations 
        self.Tampli = Tampli
        self.Tperiod = Tperiod
        self.Tphase = Tphase

        # Check if inputs are consistent
        if len(self.Tampli) != len(self.Tperiod) != len(self.Tphase):
            # If not consistent
            print('a1, a2 and a3 are not of the same size !')
        else:
            # If consistent
            content = []
            # Loop over the waves
            for i, el in enumerate(self.Tampli):
                content.append('{} {} {} \n'.format(
                    el, self.Tperiod[i], self.Tphase[i]))
            content.append('/')

            # Write the file 
            with open(self.path + 'tides.in', 'w') as file:
                for line in content:
                    file.write(line)

        return

    def create_config_file(self, sigma_dot, Vval_x1='default', 
                           Vval_x2='default', Vval_pourc=0.001, stop_crit=1, 
                           max_it=10000, final_time=1e7, tol_solver=1.00e-8, 
                           nf=False):
        """
        Create "config.in" file in the simulation directory.

        Parameters
        ----------
        sigma_dot : array 3*3 of float
            Loading rate symmetrical matrix.
        Vval_x1 and Vval_x2: float, optional
            The initial velocity perturbation is imposed between Vval_x1 and 
            Vval_x2. The default is between 0.4*Lnuc and 0.6*Lnuc.
        Vval_pourc : float, optional
            Amplitude of the initial perturbation as a percentage of V0. The 
            default is 0.001.
        stop_crit : int, optional
            Criteria to stop the simulation. Can be 0, 1 or 2.
                * 0 : simulation will stop after the first event
                * 1 : simulation will stop after max_it iterations
                * 2 : simulation will stop at final_time
            The default is 1.
        max_it : int, optional
            Number of iteration for the simulation. It will only be taken into
            consideration if stop_crit = 1. The default is 10000.
        final_time : float, optional
            Time in seconds at which the simulation will stop. It will only be 
            taken into consideration if stop_crit = 2. The default is 10000000.
        nf : int, optional
            Number of fault. If not specified, the number of fault will be 
            read in the "geometry.in" file. The default is not specified.

        Raises
        ------
        Exception
            Vval_x1 and Vval_x2 must be both specified or both set to default 
            value.            
            You must specify a fault number or create geometry.in first
        TypeError
            Stop_criteria should be 0, 1 or 2..
        ValueError
            Vval_x1 should be inferior to Vval_x2.

        Returns
        -------
        None.

        """
        # TODO ! Implement possibility to specify different a, b, Dc etc. for 
        # all faults
        
        # Store simulation parameters 
        self.sigma_dot = sigma_dot
        self.Vval_x1 = Vval_x1
        self.Vval_x2 = Vval_x2
        self.Vval_pourc = Vval_pourc
        self.stop_crit = stop_crit
        self.max_it = max_it
        self.final_time = final_time
        self.nf = nf

        # Preliminary check for Vval_x1 and Vval_x2
        if Vval_x1 == 'default' and Vval_x2 == 'default':
            self.Vval_x1 = 0.4 * self.Lnuc
            self.Vval_x2 = 0.6 * self.Lnuc
        elif Vval_x1 == 'default' and Vval_x2 != 'default':
            raise Exception(
                'Vval_x1 and Vval_x2 must be both specified or both set to \
                    default value. Here Vval_x1 is default value but Vval_x2 \
                    is not.')
        elif Vval_x1 != 'default' and Vval_x2 == 'default':
            raise Exception(
                'Vval_x1 and Vval_x2 must be both specified or both set to \
                    default value. Here Vval_x1 is default value but Vval_x2 \
                    is not.')
        
        if Vval_x1 > Vval_x2:
            raise ValueError('Vval_x1 should be inferior to Vval_x2')
            
        # Preliminary check for stop_criteria
        if stop_crit not in [0, 1, 2]:
            raise TypeError('stop_criteria should be 0, 1 or 2')

        # Computes Lb
        self.Lb = - (self.mu * self.Dc) / (self.sigma_N * self.b)

        # Get fault number if unspecified
        if self.nf is False:
            try:
                # Read "geometry.in"
                with open(self.path + 'geometry.in') as file:
                    content = file.read()
                self.nf = content.count('/')
            except Exception:  # if geometry.in do not exist
                print('You must specify a fault number or create geometry.in')
        
        # Reference velocity       
        V0_val = 1e-9  # TODO ! Option to change V0_val ?  
        
        # Compute theta_val
        self.theta_val = self.Dc / V0_val
        
        # Content of config.in file
        content = ['&main\n',
                   "simulation_name = '{}'\n".format(self.simuname),
                   '/\n',
                   '&friction\n',
                   "fric_law = '{}'\n".format(self.fric_law),
                   '/\n',
                   '&material_and_loading\n',
                   'mu = {} \n'.format(self.mu),
                   'cp = 5.00e+03 \n',
                   'cs = 3.50e+03 \n',
                   'sigma11_dot_inf = {} \n'.format(self.sigma_dot[0, 0]),
                   'sigma22_dot_inf = {} \n'.format(self.sigma_dot[1, 1]),
                   'sigma12_dot_inf = {} \n'.format(self.sigma_dot[0, 1]),
                   'sigma31_dot_inf = {} \n'.format(self.sigma_dot[0, 2]),
                   'sigma32_dot_inf = {} \n'.format(self.sigma_dot[1, 2]),
                   '/\n']

        # Create one "&fault_friction" section per fault
        for i in range(self.nf):
            content += ['&fault_friction\n',
                        "a_distr = 'CST'\n",
                        'a_val = -1.000e+00,-1.000e+00,{},{}\n'.format(
                            self.a, self.a),
                        "b_distr = 'CST'\n",
                        'b_val = -1,-1,{},{}\n'.format(self.b, self.b),
                        "Dc_distr = 'CST'\n",
                        'Dc_val = -1,-1,{},{}\n'.format(self.Dc, self.Dc),
                        "V0_distr = 'CST'\n",
                        'V0_val = -1.000e+00,-1.000e+00,1e-09,1e-09\n',
                        '/\n']
        
        # Create "&fault_initial_condition" for the first fault, where the 
        # initial perturbation will occur
        content += ['&fault_initial_condition\n',
                    "slip_distr = 'CST'\n",
                    'slip_val = -1.000e+00,-1.000e+00,0,0\n',
                    "V_distr = 'CST'\n",
                    'V_val = {},{},1e-09,{}e-09\n'.format(
                        self.Vval_x1, self.Vval_x2, 1 + self.Vval_pourc),
                    "theta_distr = 'CST'\n",
                    'theta_val = -1.000e+00,-1.000e+00,{},{}\n'.format(
                        self.theta_val, self.theta_val),
                    "sigmaN_distr = 'CST'\n",
                    'sigmaN_val = -1,-1,{},{}\n'.format(
                        self.sigma_N, self.sigma_N),
                    'ds = {}\n'.format(self.Lb/10),
                    '/\n']
        
        # Create "&fault_initial_condition" for all other faults 
        for i in range(self.nf - 1):
            content += ['&fault_initial_condition\n',
                        "slip_distr = 'CST'\n",
                        'slip_val = -1.000e+00,-1.000e+00,0,0\n',
                        "V_distr = 'CST'\n",
                        'V_val = -1.000e+00,-1.000e+00,1e-09,1e-09\n',
                        "theta_distr = 'CST'\n",
                        'theta_val = -1.000e+00,-1.000e+00,{},{}\n'.format(
                            self.theta_val, self.theta_val),
                        "sigmaN_distr = 'CST'\n",
                        'sigmaN_val = -1,-1,{},{}\n'.format(
                            self.sigma_N, self.sigma_N),
                        'ds = {}\n'.format(self.Lb/10),
                        '/\n']

        content += ['&simulation_parameter\n',
                    "fracture_mode = 'modeIII'\n",
                    "static_kernel = 'hmatrix'\n",
                    'initial_time = 0.000e+00 \n',
                    'time_step_max = 10\n',
                    'final_time = {}\n'.format(self.final_time),
                    'stop_criteria = {}\n'.format(self.stop_crit),
                    'cut_vel = 1.000e-08 \n',
                    'max_it = {} \n'.format(self.max_it),
                    'tol_solver = {} \n'.format(tol_solver),
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

        # Write config.in
        with open(self.path + 'config.in', 'w') as file:
            for line in content:
                file.write(line)

        return

    def create_geom_1fault(self, x1, x2, show=False):
        """
        Creates geometry.in file for a simple fault configuration with two 
        points (x1, 0) and (x2, 0).

        Parameters
        ----------
        path : string
            Path to the simulation directory.
        x1 : float
            x coordinate of the first point.
        x2 : float
            x coordinate of the second point.

        Returns
        -------
        None.

        """
        # Checks that x1 < x2
        if x1 > x2:
            raise ValueError('x1 should be inferior to x2')

        # Create file content
        content = ['{} 0.000e+00 \n'.format(x1),
                   '{} 0.000e+00 \n'.format(x2),
                   '/\n']

        # Write geometry.in
        with open(self.path + 'geometry.in', 'w') as file:
            for line in content:
                file.write(line)

        # Plot the geometry
        if show:
            self.plot_geometry()
            
        return

    def create_geom_2_faults_overlapping(self, D_over_Lnuc, L_over_Lnuc, 
                                         overlap, show=False):
        """
        Create geometry.in file for a 2 faults overlapping geometry (as in 
        Romanet et al. (2018), GRL).

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
            Length of the fault express as the ratio L/Lnuc.
        overlap : real
            Portion of the fault overlapping. 0 < overlap < 1.
            overlap = 0:
                      =======
               ======= 
            overlap = 1:
                =======
                =======

        Returns
        -------
        None.

        """
        self.L_over_Lnuc = L_over_Lnuc

        # Preliminary checks
        if (overlap > 1) or (overlap < 0):
            raise ValueError('overlap should be between 0 and 1')

        # Create file content
        content = ['0.000e+00 0.000e+00 \n',
                   '{} 0.000e+00 \n'.format(L_over_Lnuc * self.Lnuc),
                   '/\n',
                   '{} {} \n'.format(L_over_Lnuc * self.Lnuc *
                                     (1 - overlap), D_over_Lnuc * self.Lnuc),
                   '{} {} \n'.format(L_over_Lnuc * self.Lnuc *
                                     (2 - overlap), D_over_Lnuc * self.Lnuc),
                   '/\n']

        # Write geometry.in
        with open(self.path + 'geometry.in', 'w') as file:
            for line in content:
                file.write(line)

        # Plot the geometry
        if show:
            self.plot_geometry()
            
        return

    def create_geom_multiple_faults(self, lengths, angles, xs, ys, show=False):
        """
        Create geometry.in file for a multiple faults geometry.
        
              /
             /        
            /
           /
          + <-- Point defining the fault (given as x and y distances from the 
                origin). Here the angle is positive and is ~60°.
                 
        
        Faults are defined with a length, an angle and the distance in x (xs) 
        and y (ys) of one point from the first fault. 
        Hence len(lengths) = len(angles) = n and len(xs) = len(ys) = n-1, n is 
        the number of faults.
        Length, x and y distances are normalised by Lnuc.
        The angle is positive in the trigonometric direction. 
        The point defining the first fault is (0, 0).

        Parameters
        ----------
        lengths : list of float
            Lengths of the faults normalised by Lnuc (L/Lnuc).
        angles : list of float
            Angles of the faults.
        xs : list of float
            Distance in x normalised by Lnuc between the points defining the 
            faults and the origin.
        ys : list of float
            Distance in y normalised by Lnuc between the points defining the 
            faults and the origin.
        show : bool, optional
            If True, plot the fault system geometry. The default is False.

        Returns
        -------
        None.

        """
        # Preliminary checks
        if len(lengths) != len(angles):
            raise Exception('lengths and angles must have the same size !')
        if len(xs) != len(ys):
            raise Exception('xs and ys must have the same size !')
        if len(xs) != len(lengths) - 1:
            raise Exception('xs and ys must have one element less than \
                            lengths and angles')            
                
        # Initialisation
        xs = [0] + xs
        ys = [0] + ys
        
        # Compute the second edges of the faults 
        x2 = (xs + np.cos(np.deg2rad(angles)) * lengths) * self.Lnuc
        y2 = (ys + np.sin(np.deg2rad(angles)) * lengths) * self.Lnuc
        
        # Assemble first edges (xs and ys) with second edges
        X = [None] * 2 * len(xs)  # Initialisation
        Y = [None] * 2 * len(ys)  # Initialisation
        X[::2] = [el * self.Lnuc for el in xs]
        X[1::2] = [el for el in x2]
        
        Y[::2] = [el * self.Lnuc for el in ys]
        Y[1::2] = [el for el in y2]
        
        # Create the content of the file
        content = []
        for i, el in enumerate(X[::2]):
            content += ['{} {} \n'.format(el, Y[i*2]),
                        '{} {} \n'.format(X[(i*2)+1], Y[(i*2)+1]),
                        '/\n']
        
        # Write geometry.in
        with open(self.path + 'geometry.in', 'w') as file:
            for line in content:
                file.write(line)
        
        # Plot the geometry
        if show:
            self.plot_geometry()

        return

    def plot_geometry(self, scale='Lnuc', savefig=True):
        """
        Plot the fault system geometry from geometry.in.

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
        # Read geometry.im
        with open(self.path + 'geometry.in', 'r') as file:
            content = file.readlines()
            
        # Extract nbr_fault, x amd y coordinates
        # Initialisation
        nbr_fault = 0
        x = []
        y = []
        for line in content:
            if line == '/\n':
                nbr_fault += 1
            else:
                x_temp, y_temp = line.split()
                x.append(float(x_temp))
                y.append(float(y_temp))
        # Store results        
        self.nbr_fault = nbr_fault
        self.x = x
        self.y = y
        
        # Compute figure limits 
        ymaxi = np.max(self.y)
        ymini = np.min(self.y)
        xmaxi = np.max(self.x)
        xmini = np.min(self.x)
        if scale == 'Lnuc':        
            ymini = ymini / self.Lnuc
            ymaxi = ymaxi / self.Lnuc
            xmini = xmini / self.Lnuc
            xmaxi = xmaxi / self.Lnuc
        Ly = ymaxi-ymini # extend in y direction of the fault system
        Lx = xmaxi-xmini # extend in x direction of the fault system

        # Initialise figure
        fig, ax = plt.subplots(1, 1)
        
        # Plot each fault 
        for i in range(nbr_fault):
            if scale == 'Lnuc':
                x = [el / self.Lnuc for el in self.x[i*2:i*2+2]]
                y = [el / self.Lnuc for el in self.y[i*2:i*2+2]]
            else:
                x = self.x[i*2:i*2+2]
                y = self.y[i*2:i*2+2]
            
            ax.plot(x, y, 'b')
        
            # Add fault "label"
            xtext = 0.5 * (np.max(x) + np.min(x))
            ytext = 0.5 * (np.max(y) + np.min(y))
            ax.text(xtext, ytext, 'F{}'.format(i+1), bbox=dict(fc='yellow', 
                    ec='none', pad=1), ha='center', va='center')
            
            # # Set axis limits and aspect 
            if ymaxi == ymini :  # If all faults are aligned horizontally
                ax.set_ylim(ymini -Lx*0.1, ymaxi + Lx*0.1)  
            elif xmaxi == xmaxi :  # If all faults are aligned vertically 
                ax.set_xlim(xmini -Ly*0.1, xmaxi + Ly*0.1) 
            else:
                ax.set_ylim(ymini - Ly*0.1, ymaxi + Ly*0.1)
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
        
        fig.tight_layout()
        
        if savefig:
            fig.savefig(self.path + 'geometry_in.png', dpi=400, \
                        bbox_inches='tight')
    
        return
