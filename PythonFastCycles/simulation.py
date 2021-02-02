"""
Class to create all files for a simulation.

First version by R. Ferry on January 2021.
"""
import os

class Simulation:
    def __init__(self, path, simuname, mu, a, b, fric_law='RateStateAgeing_R',\
                 frac_mode='ModeIII', sigma_N=-1e8, Dc=1e-3):
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

        # Compute Lnuc
        self.Lnuc = - (self.mu * self.Dc) / (self.sigma_N * (self.b - self.a))

    def create_all_files(self, L_over_Lnuc, sigma_dot, geom_type, \
                         D_over_Lnuc=0.1, overlap=0.5, GPSx=[10], GPSy=[10], \
                         Tampli=[0.0], Tperiod=[0.0], Tphase=[0.0], \
                         Vval_x1='default', Vval_x2='default', \
                         Vval_pourc=0.001, stop_crit=1, max_it=10000, \
                         final_time=10000000, nf=False):
        """
        Create all files for a simulation, i.e. "config.in", "geometry.in", 
        "tides.in" and "GPS.in".

        Parameters
        ----------
        L_over_Lnuc : real
            Length of the fault express as the ratio L/Lnuc.
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

        Raises
        ------
        Exception
            DESCRIPTION.

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
            self.create_geom_1fault(x1, x2)
        elif geom_type == '2faults_overlapping':
            self.create_geom_2_faults_overlapping(
                D_over_Lnuc, L_over_Lnuc, overlap)
        else:
            raise Exception('geom_type does not exist.')

        self.create_config_file(sigma_dot, Vval_x1, Vval_x2, Vval_pourc, \
                                stop_crit, max_it, final_time, nf)

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

    def create_tides_file(self, Tampli=[0.0], Tperiod=[0.0], Tphase=[0.0]):
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

    def create_config_file(self,sigma_dot,Vval_x1='default',Vval_x2='default',\
                           Vval_pourc=0.001, stop_crit=1, max_it=10000, \
                           final_time=10000000, nf=False):
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

        # Write config.in
        with open(self.path + 'config.in', 'w') as file:
            for line in content:
                file.write(line)

        return

    def create_geom_1fault(self, x1, x2):
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

        return

    def create_geom_2_faults_overlapping(self, D_over_Lnuc, L_over_Lnuc, overlap):
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

        return
