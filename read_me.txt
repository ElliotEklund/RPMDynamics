RPMDynamics

This code is written to produce basic RPMD simulations, primiarly those
referenced in D. E. Manolopoulos and I. R. Craig's 2004 paper.  

(1) Potential Function
The potential on which the ring polymer lives is contained in
"functions.cpp". Specifically, the functions "potential" and "dVdq" need to be
derived by hand and hardcoded.

(2) Equilibrium Population Generation 
Standard Path Integral Monte Carlo (MC) methods are used to generate an 
equilibrium distribution from which dynamic trajectories are sampled. 
ASTRA does not currently support compilers capable of implementing the 
most advanced random number generator (RNG) features available to the C++
language. Thus, RNGs are included from the classic text "Numerical Recipes" 3rd ed. 
These include files must be included in a specific order!

(3) Population Sampling
After equilibrium has been established, phase space points are sampled and
added to coordinate and momentum vectors. These points are sampled at a
frequency determined by the decorrelation length. Additionally, the
decorrelation length is slightly adjusted by a RNG to help decrease the
likelihood two trajectories are artificially correlated.

(4) Dynamics via Velocity Verlet
A Velocity Verlet (VV) integrator is used to propagate classial trajectories.
Correlation functions should be calculated in this section of the code.
Currently, RPMDynamics is equiped to handle auto-correlation functions, but
modifications can readily be made.

(5) Output and Plotting
Two .csv files are produced: "mcCalculations" and "dynCalculations". To plot
csv files, open gnuplot and change the datafile separator to ','. That is, type:

	set datafile separator ','

Any modifications to the correlation function for thermodynamic estimators
will be output in these files.

Good Luck!

Elliot Eklund
