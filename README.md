# AdvPgrm_PM
The PM main repo

 A C code is used to simulate the evolution of dark matter. We use a lagrangian approach through the particle mesh method.

 In order to run it correctly, as of now, you firstly have to run the makefile "makefileInit", then run the Init executable, then run the makefile "makefile" to compile the Main, to finally run Main.
  
examples of bash commands required in order to run:
 make -f makefileInit
 ./Init
 make 
 ./Main

The numbers of step of the leap frog algorithm is still a sore point of this program. Atm you can edit the number of time step taken through the variabile MAX_STEPS inside the Main.c file in line 19.

What to expect?
 As of today this program, if correctly compiled, should produce:
  -  a "DISTRRHOB.dat" binary file containing the distributed particle inside the box.*
  -  a "DISTRRHO" text file containing the same as above*
  -  a "phaseSpace.txt" text file containing the temporal evolution of 5 particles inside the phase space (x - v). This data formatting is meant to be used with gnuplot. When plotting with gnuplot you can use the "u 2n-1:2n" keyword to plot the n-th particle evolution. Eg, to plot the 3rd particle trajectory in the phase space you should type:
 "plot "phaseSpace.txt" u 5:6 w l".

*these files creation is to be set in the param.txt file. Default is to create both .dat and the text file.*

The "makefile_hard" is the future version of the makefile file. It is meant to provide different keyword based on the need of the user, eg run/debug/run with omp.

This is, clearly, a work in progress..
