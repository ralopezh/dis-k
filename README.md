# Dispersion Solver for drifting bi-Kappa plasmas

## Author
Rodrigo A. López

Universidad de Santiago de Chile, Usach

rlopez186@gmail.com

ORCID: [0000-0003-3223-1498](https://orcid.org/0000-0003-3223-1498)

## Credits
If you use this code, please acknowledge the following references:

[1] - López, R., Shaaban, S., & Lazar, M. (2021). General dispersion
properties of magnetized plasmas with drifting bi-Kappa distributions.
DIS-K: Dispersion Solver for Kappa Plasmas. Journal of Plasma
Physics,87(3), 905870310.
doi:[10.1017/S0022377821000593](https://doi.org/10.1017/S0022377821000593)

[2] - López, R.A., Moya, P.S., Shaaban, S.M., Lazar, M., Yoon, P.H.,
Poedts, S. (2021). Advanced Numerical Tools for Studying Waves and
Instabilities in Kappa Distributed Plasmas. In: Lazar, M., Fichtner,
H. (eds) Kappa Distributions. Astrophysics and Space Science Library,
vol 464. Springer, Cham.
doi:[10.1007/978-3-030-82623-9_9](https://doi.org/10.1007/978-3-030-82623-9_9)

This solver is programmed based on the equations derived in the above references.


## Requisites

* Make

* Fortran compiler \ 
gfortran or ifort \
You can select the fortran compiler in the Makefile. You can also choose to include the openMP flag to run in parallel.

* python3 \
  If you want a plot after the program is finished. If you prefer to make the plot yourself, comment lines 98 and 152 in src/main.f90

## Compile
Just run
```
make
```
The executable will be created in the main folder. Then, go to the example folder (or any other folder) and run the code using an input file.

## Routines

Description of the routines used in src/

* main.f90 : The main program produces executable dis-k.

* nyquist.f90 : Secondary program, produce executable nyquist. Run this code to take the initial seed at k0, and include it in the file input.dat.

* params.f90 : Modules to be used.

* load.f90 : Initialize all the variables using the input file.

* plasmadf.f90 : Definition of the plasma dispersion function.

* muller.f90 : Computes the real or complex roots of a nonlinear function, using an optimal Muller's method.

* disp.f90 : Compute the elements of the dispersion tensor and return the determinant of it.
* solve.f90 : Solve the equation $\text{Det}{\Lambda}=0$ using the muller routine. For a fixed value of theta, return the corresponding complex frequency for the entire k vector.

* polarization.f90 : Compute the polarization and ratios $|Ei/E|^2$
and $|Bi/B|^2$

* solve2d.f90 : Use solve.f90 to find the complex frequencies for the complete range in k and theta. Creates a 2d vector theta vs. k.

* extras.f90 : Contains multiple definitions for exponentially normalized modified Bessel functions.
               
## Run
If you do not have an initial guess solution, run first
```
path_to_executable/nyquist input.dat
```
If your input file contains the initial guess, run
```
path_to_executable/dis-k input.dat
```

## How it works

First, you need to set an input file. You can use the examples contained in the folder examples. 
As in any linear dispersion solver, you must provide a complex initial
seed $(w_{0r},w_{0i})$ for the solution at the initial wavenumber, k0, and angle, th0. If you don't have an initial guess, leave those values blank and run the code nyquist using that input file.
The code nyquist plot the contours of  $\text{Re}{\text{Det}(\Lambda)}
= 0$ and $\text{Im}{\text{Det}(\Lambda)}=0$ in a $\text{Re}(w)$ vs.
$\text{Im}(w)$ map.\
All the crossings between those contours correspond to solutions of
the $\text{Det}(\Lamnda)=0. If you don't find intersections, you can adapt the search range at the end of the input file.
Then, you can extract the values of $\text{Re}(w)$ and $\text{Im}(w)$ at those intersections and put them in the input file.
For more details, please read Ref. [2] above.

Now you can run dis-k using the input file. If you are solving for a fixed angle, nth=1, it will create an output.save file and plot using plot_disp.py.
If you solve for a range of angles, nth>1, it will create an output_kth.save and plot using plot_map.py.

