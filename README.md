# PythonFastCycles 

This is a Python package to interact with the the Fortran code *FastCycles*.      
*FastCycles* has been developped by Pierre Romanet during his PhD thesis under the direction of R. Madariaga and H. Bhat. See Romanet et al. (2018). [[1]](#1)

## Install

You can install the last version of *PythonFastCycles* using pip.    
Run in a terminal:         
`pip install PythonFastCycles`

## How to use

### Create file for a simulation 

#### Defines path to a simulation folder
`path = 'path/to/fastcycles/problems/' `

#### Defines parameters
```python
fric_law = 'RateStateAgeing_R'
frac_mode = 'ModeIII'

mu = 3e10 
Dc = 1e-3
sigma_N = -1e8
a = 0.0075
b = 0.01
```

The lenght of the fault L is defined with the ratio L/L<sub>nuc</sub> (L<sub>nuc</sub> is computed automatically).      

`L_over_Lnuc = 2`

Defines \\sigma 

```python
s11 = 0.00e+00
s22 = 0.00e+00
s33 = 0.00e+00
s12 = 0.00e+00
s13 = 0.00e+00
s23 = 0.1

sigma_dot = np.array([[s11, s12, s13], [s12, s22, s23], [s13, s23, s33]])
```

## References
<a id="1">[1]</a> 
Romanet et al. (2018). 
ast and slow slip events emerge due to fault geometrical complexity. 
*Geophysical Research Letters*, 45(10), pp.4809-4819
https://doi.org/10.1029/2018GL077579
