Dispersion Solver for drifting bi-Kappa plasmas


Description of the routines used

- main.f90 : Main program, produce executable dis-k.
- nyquist.f90 : Secundary program, produce executable nyquist.
                Run this code to take the initial seed at k0, and
		include it in the file input.dat.
- params.f90 : Modules to be used.
- load.f90 : Initialize all the variables using the input file.
- plasmadf.f90 : Definition of the plasma dispersion function.
- root_finder.f90 : Computes the real or complex roots of a nunlinear
                    function, using an optimal Muller's method.
- disp.f90 : Compute the elements of the dispersion tensor and return
             the determinant of it.
- solve.f90 : Solve the equation Det{\Lambda}==0 using the root_finder.
              for a fixed value of theta, return the corresponding
	      complex frequency for the entire k vector. 
- polarization.f90 : Compute the polarization and ratios |Ei/E|^2 and 
                     |Bi/B|^2
- solve2d.f90 : Use solve.f90 to find the complex frequencies for the
                complete range in k and theta. Creates a 2d vector
		theta vs. k.
- extras.f90 : Contains multiple definitions for exponentially
               normalized modified Bessel functions.
              
