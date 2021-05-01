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
                         D_over_Lnuc=0.1, overlap=0.5, W_over_Lnuc=0.1, 
                         Tampli=[0.0], Tperiod=[1.0], Tphase=[0.0], GPSx=[10],
                         GPSy=[10], Vval_x1='default', Vval_x2='default',  
                         Vval_pourc=0.001, stop_crit=None, max_it=10000, 
                         final_time=10000000, tol_solver=1.00e-8, nf=False,
                         times=None, s_amplitudes=None, n_amplitudes=None,
                         version=14, amplitude_decimation_factor=1, 
                         shear_magnitude=0, normal_magnitude=0):
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
                      
                * "2faults_aligned" for a configuration with two faults aligned
                  and separated by a distance W_over_Lnuc.
                  
                    L/Lnuc
                   <------>
                   ========      ========
                           <---->
                           W/Lnuc
                           
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
            Criteria to stop the simulation. Can be (version<14 / version>=14):
                * 0 / 'nucleation' : simulation will stop after the first event
                * 1 / 'time steps' : simulation will stop after max_it 
                  iterations
                * 2 / 'time' : simulation will stop at final_time
            The default is 1 / 'time steps'.
        max_it : int, optional
            Number of iteration for the simulation. It will only be taken into
            consideration if stop_crit = 1. The default is 10000.
        final_time : float, optional
            Time in seconds at which the simulation will stop. It will only be 
            taken into consideration if stop_crit = 2. The default is 10000000.
        nf : int, optional
            Number of fault. If not specified, the number of fault will be 
            read in the "geometry.in" file. The default is not specified.
        times : list, optional
            Times in seconds associated with the change of amplitude. Default 
            is [0,0].
        s_amplitudes : list, optional
            Amplitudes of the shear stress as a percentage of the 
            background loading given in "config.in". Default is [1,1]. 
        n_amplitudes : list, optional
            Amplitudes of the normal stress as a percentage of the normal 
            traction given in "config.in". Default is [1,1]. 
        version : int, optional
            Version of FastCycles. Default is 14.
        amplitude_decimation_factor : float, optional
            Decimation factor of the optimal time step computed for amplitude.
        shear_magnitude : float, optional
            Shear traction is given as background + shear_magnitude*amplitude,
            with amplitude between 0 and 1.
        normal_magnitude : float, optional
            Normal traction is given as background + normal_magnitude*amplitude
            with amplitude between 0 and 1.

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
        self.version = version 
        self.create_GPS_file(GPSx, GPSy)
        if self.version < 14:
            self.create_tides_file(Tampli, Tperiod, Tphase)
            self.create_amplitude_file(times, s_amplitudes)
        else:
            self.create_amplitude_file_version14(times, s_amplitudes, 
                                                 n_amplitudes)   

        if geom_type == '1fault':
            x1 = 0
            x2 = self.Lnuc * self.L_over_Lnuc
            self.create_geom_1fault(x1, x2, show=show)
        elif geom_type == '2faults_overlapping':
            self.create_geom_2_faults_overlapping(
                D_over_Lnuc, L_over_Lnuc, overlap, show=show)
        elif geom_type == '2faults_aligned':
            self.create_geom_2faults_aligned(W_over_Lnuc, L_over_Lnuc, 
                                             show=show)
        elif geom_type == "multiple":
            self.create_geom_multiple_faults(lengths, angles, xs, ys, 
                                             show=show)
        else:
            raise Exception('geom_type does not exist.')

        self.create_config_file(sigma_dot=sigma_dot, Vval_x1=Vval_x1, 
                                Vval_x2=Vval_x2, Vval_pourc=Vval_pourc, nf=nf, 
                                stop_crit=stop_crit, max_it=max_it, 
                                final_time=final_time, tol_solver=tol_solver,
                                amplitude_decimation_factor=amplitude_decimation_factor, 
                                shear_magnitude=shear_magnitude,
                                normal_magnitude=normal_magnitude)

    def create_GPS_file(self, GPSx=[10], GPSy=[10]):
        """
        Create GPS.in at the location 'path'.

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
                           Vval_x2='default', Vval_pourc=0.001, max_it=10000, 
                           stop_crit=None, final_time=1e7, tol_solver=1.00e-8, 
                           nf=False, amplitude_decimation_factor=1, 
                           shear_magnitude=0, normal_magnitude=0):
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
            Criteria to stop the simulation. Can be (version<14 / version>=14):
                * 0 / 'nucleation' : simulation will stop after the first event
                * 1 / 'time steps' : simulation will stop after max_it 
                  iterations
                * 2 / 'time' : simulation will stop at final_time
            The default is 1 / 'time steps'.
        max_it : int, optional
            Number of iteration for the simulation. It will only be taken into
            consideration if stop_crit = 1. The default is 10000.
        final_time : float, optional
            Time in seconds at which the simulation will stop. It will only be 
            taken into consideration if stop_crit = 2. The default is 10000000.
        nf : int, optional
            Number of fault. If not specified, the number of fault will be 
            read in the "geometry.in" file. The default is not specified.
        amplitude_decimation_factor : float, optional
            Decimation factor of the optimal time step computed for amplitude.
        shear_magnitude : float, optional
            Shear traction is given as background + shear_magnitude*amplitude,
            with amplitude between 0 and 1.
        normal_magnitude : float, optional
            Normal traction is given as background + normal_magnitude*amplitude
            with amplitude between 0 and 1.
            
        Raises
        ------
        Exception
            Vval_x1 and Vval_x2 must be both specified or both set to default 
            value.            
            You must specify a fault number or create geometry.in first
        TypeError
            Stop_criteria should be 0, 1 or 2 (version < 14) or 'nucleation', 
            'time steps', 'time (version >= 14).
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
        self.amplitude_decimation_factor = amplitude_decimation_factor
        self.shear_magnitude = shear_magnitude 
        self.normal_magnitude = normal_magnitude 

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
        
            
        # Preliminary check for stop_criteria and set default value
        if self.version >= 14:
            if stop_crit is None:  # if it is the default value
                stop_crit = 'time steps'
            elif stop_crit not in ['nucleation', 'time steps', 'time']:
                raise TypeError("stop_criteria should be 'nucleation', 'time steps', 'time' ")
        else:
            if stop_crit is None:  # if it is the default value
                stop_crit = 1
            elif stop_crit not in [0, 1, 2]:
                raise TypeError("stop_criteria should be 0, 1 or 2")
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
                   'sigma32_dot_inf = {} \n'.format(self.sigma_dot[1, 2])]
        if self.version >= 14:
            content += ['shear_magnitude = {} \n'.format(self.shear_magnitude),
                        'normal_magnitude = {} \n'.format(
                            self.normal_magnitude),
                        '/\n']
        else:
            content += ['/\n']

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
                    'final_time = {}\n'.format(self.final_time)]
        
        if self.version >= 14:
            content += ["stop_criteria = '{}'\n".format(self.stop_crit)]
        else:
            content += ['stop_criteria = {}\n'.format(self.stop_crit)]
            
        content += ['cut_vel = 1.000e-08 \n',
                    'max_it = {} \n'.format(self.max_it),
                    'tol_solver = {} \n'.format(tol_solver),
                    'tol_interp = 1.000e-08 \n',
                    'iprec = 4\n',
                    'icheck_interp = 0\n',
                    'omp_threads = 4\n']
        if self.version >= 14:
            content += ['amplitude_decimation_factor = {}\n'.format(
                self.amplitude_decimation_factor),
                        '/\n']
        else:
            content += ['/\n']
        
        content += ['&output\n',
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
    
    def create_amplitude_file(self, times=None, s_amplitudes=None):
        """
        Create "amplitude.in" file in the simulation directory. Default is no 
        variation of the background loading rate.

        Parameters
        ----------
        times : list, optional
            Times in seconds associated with the change of amplitude. Default 
            is [0, 1].
        s_amplitudes : list, optional
            Amplitudes of the background loading rate as a percentage of the 
            loading given in "config.in". Default
            is [1, 1].
        
        Returns
        -------
        None.

        """            
        # Deal with default value
        if times is None:
            times = [0, 1]
            s_amplitudes = [1, 1]
        else:
            # Check that times and s_amplitudes have the same length
            if len(times) != len(s_amplitudes):
                raise Exception('times and s_amplitudes must have the same length')
            
        # Create content
        content = []
        for i, time in enumerate(times):
            content.append('{} {} \n'.format(time, s_amplitudes[i]))
        content.append('/\n')
        
        # Write amplitude.in
        with open(self.path + 'amplitude.in', 'w') as file:
            for line in content:
                file.write(line)
                
    def create_amplitude_file_version14(self, times=None, s_amplitudes=None, 
                                        n_amplitudes=None):
        """
        Create "amplitude.in" file in the simulation directory for version 14 
        of FastCycles. Default is no variation of the shear and normal 
        amplitudes.

        Parameters
        ----------
        times : list, optional
            Times in seconds associated with the change of amplitude. Default 
            is [0,1].
        s_amplitudes : list, optional
            Amplitudes of the shear stress as a percentage of the 
            background loading given in "config.in". Default is [0,0]. 
        n_amplitudes : list, optional
            Amplitudes of the normal stress as a percentage of the normal 
            traction given in "config.in". Default is [0,0].     

        Returns
        -------
        None.

        """
        # Deal with default values
        if times is None:  # if times is None = if nothing is specify
            times = [0,1]
            s_amplitudes=[0,0]
            n_amplitudes=[0,0]
        # Varying normal stress but constant background loading
        elif s_amplitudes is None:  # time is given then n_amplitudes is also given
            s_amplitudes = [0] * len(n_amplitudes) 
        # Varying background loading but constant normal stress    
        elif n_amplitudes is None:  # time is given then s_amplitudes is also given
            n_amplitudes = [0] * len(s_amplitudes)
            
        # Create content 
        content = []
        for i, time in enumerate(times):
            content.append('{} {} {} \n'.format(time, s_amplitudes[i], 
                                                n_amplitudes[i]))
        content.append('/\n')
        
        # Write amplitude.in
        with open(self.path + 'amplitude.in', 'w') as file:
            for line in content:
                file.write(line)

    def create_geom_1fault(self, x1, x2, show=False):
        """
        Create geometry.in file for a simple fault configuration with two 
        points (x1, 0) and (x2, 0).

        Parameters
        ----------
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

        # Create content
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
    
    def create_geom_2faults_aligned(self, W_over_Lnuc, L_over_Lnuc, 
                                    show=False):
        self.W_over_Lnuc = W_over_Lnuc 
        self.L_over_Lnuc = L_over_Lnuc
        
        # Create file content
        content = ['0.000e+00 0.000e+00 \n',
                   '{} 0.000e+00 \n'.format(L_over_Lnuc * self.Lnuc),
                   '/\n',
                   '{} 0.000e+00 \n'.format((L_over_Lnuc + W_over_Lnuc) * self.Lnuc),
                   '{} 0.000e+00 \n'.format((2 * L_over_Lnuc + W_over_Lnuc) * self.Lnuc),
                   '/\n']
        
        # Write geometry.in
        with open(self.path + 'geometry.in', 'w') as file:
            for line in content:
                file.write(line)

        # Plot the geometry
        if show:
            self.plot_geometry()

    
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

    def plot_geometry(self, scale='Lnuc', GPS=True, savefig=True):
        """
        Plot the fault system geometry from geometry.in.

        Parameters
        ----------
        scale : str, optional
            Scale of axis. Can be 'Lnuc' (normalised by Lnuc) or 'X'. The 
            default is 'Lnuc'.
        GPS : bool, optional
            If GPS is True, plot the GPS stations. Default is True.
        savefig : bool, optional
            If savefig is True, save the figure in the simulation directory 
            under the name "slip_rate_evolution.png". The default is True.

        Returns
        -------
        None.

        """
        # TODO ! Add GPS station plot
        # Read geometry.in
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
        # Faults limits
        ymaxi_faults = np.max(self.y)
        ymini_faults = np.min(self.y)
        xmaxi_faults = np.max(self.x)
        xmini_faults = np.min(self.x)
        if scale == 'Lnuc':        
            ymini_faults = ymini_faults / self.Lnuc
            ymaxi_faults = ymaxi_faults / self.Lnuc
            xmini_faults = xmini_faults / self.Lnuc
            xmaxi_faults = xmaxi_faults / self.Lnuc
        Ly = ymaxi_faults-ymini_faults # fault system extend in y direction 
        Lx = xmaxi_faults-xmini_faults # fault system extend in x direction 
        # GPS limits
        if GPS:
            ymaxi_GPS = np.max(self.GPSy)
            ymini_GPS = np.min(self.GPSy)
            xmaxi_GPS = np.max(self.GPSx)
            xmini_GPS = np.min(self.GPSx)
            
            if scale == 'Lnuc':
                ymaxi_GPS = ymaxi_GPS / self.Lnuc
                ymini_GPS = ymini_GPS / self.Lnuc
                xmaxi_GPS = xmaxi_GPS / self.Lnuc
                xmini_GPS = xmini_GPS / self.Lnuc
             
            ymaxi = max(ymaxi_faults, ymaxi_GPS)
            ymini = min(ymini_faults, ymini_GPS)
            xmaxi = max(xmaxi_faults, xmaxi_GPS)
            xmini = min(xmini_faults, xmini_GPS)
        else :
            ymaxi = ymaxi_faults
            ymini = ymini_faults
            xmaxi = xmaxi_faults
            xmini = xmini_faults
                
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
            elif xmaxi == xmini :  # If all faults are aligned vertically 
                ax.set_xlim(xmini -Ly*0.1, xmaxi + Ly*0.1) 
            else:
                ax.set_ylim(ymini - Ly*0.1, ymaxi + Ly*0.1)
                ax.set_xlim(xmini - Lx*0.1, xmaxi + Lx*0.1)
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
        
        # Plot GPS stations
        if scale == 'Lnuc':
            ax.scatter([el / self.Lnuc for el in self.GPSx], 
                       [el/ self.Lnuc for el in self.GPSy], color='red', s=20)
            # Add number next to the station
            for i, x in enumerate(self.GPSx):
                ax.annotate(i+1, (x / self.Lnuc + 0.03, 
                            self.GPSy[i] / self.Lnuc + 0.03), color='red')
        else :
            ax.scatter(self.GPSx, self.GPSy, color='red', s=20)
            # Add number next to the station
            for i, x in enumerate(self.GPSx):
                ax.annotate(i+1, (x + 7, self.GPSy[i] + 7), color='red')
        
        fig.tight_layout()
        
        if savefig:
            fig.savefig(self.path + 'geometry_in.png', dpi=400, \
                        bbox_inches='tight')
    
        return
