# Atomistic spin dynamics

_Modelling of magnetisation dynamics at the atomistic scale._

## Summary

The code represents the implementation of the atomistic spin dynamics model[^1] (ASD), which describes the time evolution of magnetic moments of interacting atoms. The lattice positions of the atoms are fixed, and the behaviour of the magnetic moments is governed by the Landau-Lifshitz-Gilbert equation. The latter is solved using the mid-point method[^2]. The atomistic magnetic moments are represented by unit vectors in 3D; the crystallographic lattice can be arbitrary (1D/2D/3D).

The code was originally written for Ref.[^3]. It subsequently evolved over the years, acquiring numerous improvements. The main highlights of the code are:
- the following energy terms are currently implemented: the exchange interaction, the magnetocrystalline anisotropy, the interaction with the external field;
- thermal fluctuations are incorporated via the Langevin dynamics;
- the code has minimalist design.

In the current version of the code excludes the implementation of the long-range dipolar interactions. The latter has been separately implemented in the currently-unsupported branch of the code and can be shared upon request.

The code is written in MATLAB.

## Getting started

The user can run the code by executing the main file and providing the file name where the results should be stored: `mainAtom('tmp1');`

The default example simulates the propagation of a domain wall in a 2D triangular lattice with a notch (three missing atoms).

After the calculations are finished (the default example takes approximately 150 seconds on a laptop with 11th Gen Intel Core i5), the user can run the postprocessing script to create a video file visualising the solution: `animateField;`

For understanding the problem parameters and the setup of the example, the user is referred to Ref.[^3].

## Developer and acknowledgements

The code has been developed by Dr. M. Poluektov during the postdoctoral fellowship at Uppsala University. 

[^1]: B. Skubic _et al._, _J. Phys.: Condens. Matter_ 20(31):315203, 2008, [link](https://iopscience.iop.org/article/10.1088/0953-8984/20/31/315203).
[^2]: M. d'Aquino _et al._, _J. Comput. Phys._ 209(2):730-753, 2005, [link](https://www.sciencedirect.com/science/article/pii/S002199910500197X).
[^3]: M. Poluektov _et al._, _Comput. Methods Appl. Mech. Eng._ 329:219-253, 2018, [link](https://www.sciencedirect.com/science/article/pii/S0045782517302463).