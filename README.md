# Waveform Amplitudes for a Spinning Particle (WASP)

Mathematica package for calculating the gravitational wave amplitudes emitted by a spinning test particle in Kerr spacetime. 
The code provides the adiabatic and linear-in-spin corrections calculated in the Teukolsky formalism.
At the moment, only nearly-equatorial orbits are implemented

Mathematica files
- ```WASP.wl```: Mathematica package for the calculation of the numerical solutions to the radial Teukolsky equation.

Installation
The file ```HSCSolverHomogeneousRadialTeukolsky.wl``` must be placed at one of the paths shown in the variables $BaseDirectory (for system-wide installations) and $UserBaseDirectory (for single-user installations). 
Place it in the Applications/ subdirectory (or subfolder) of the paths returned by those variables. Alternatively, you can load the package using the command

```Get["/absolute_path_where_the_package_is_located/HSCSolverHomogeneousRadialTeukolsky.wl"]; ``` 

Required Mathematica packages:
- ```HSCSolverHomogeneousRadialTeukolsky.wl```. Numerical solver for the homogeneous radial Teukolsky equation.
      Available at https://github.com/gabriel-andres-piovano/HSCSolverHomogeneousRadialTeukolsky
- ```AnalyticSpinOrbitsHamiltonJacobiSpectral.wl```. Spectral solver for the spin corrections to the motion.
      Available at https://github.com/gabriel-andres-piovano/Analytic_Solutions_Spinning-Body-Hamilton-Jacobi
- ```SpinWeightedSpheroidalHarmonics``` of the BlackHolePerturbationToolkit https://github.com/BlackHolePerturbationToolkit/SpinWeightedSpheroidalHarmonics
