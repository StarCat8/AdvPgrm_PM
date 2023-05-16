# AdvPgrm_PM
The PM main repo

 A C code is used to simulate the evolution of dark matter. We use a lagrangian approach through the particle mesh method.

 In order to run it just run the RUN.sh bash script.
 
Else you could just type the following in the bash terminal:
 make -f makefileInit
 ./Init
 make 
 ./Main

The numbers of step of the leap frog algorithm is still a sore point of this program. Atm you can edit the number of time step taken through the variabile MAX_STEPS inside the Main.c file in line 19.

What to expect?
 As of today this program, if correctly compiled, should produce:
  -  a "DISTRRHOB.dat" binary file containing the distributed particle inside the box.*
  -  a "DISTRRHO" text file containing the same as above*
  -  5 ("phaseSpace1",..., "phaseSpace5") files containing a screenshot of the phase space of the 4096 particles. Plotting one over the other (using gnuplot) allows us to appreciate the evolution of the system.

**these files creation is to be set in the param.txt file. Default is to create both .dat and the text file.*

The "makefile_hard" is the future version of the makefile file. It is meant to provide different keyword based on the need of the user, eg run/debug/run with omp.

This is, clearly, a work in progress..
