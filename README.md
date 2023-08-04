# AdvPgrm_PM
The PM main repo

 A C code is used to simulate the evolution of dark matter. We use a lagrangian approach through the particle mesh method.

 To run the code first build the file with the "make" command.
Then first run the "Init" executable in order to generate the initial distribution of points. 
 Then you can start the "Main" executable.
 
The Main exe require an argument at the time of call, namely either -NGP -CIC -TSC. This will select the kernel for distributing the points on the grid. Once the program is running the user is required to input either 0 or 1 (though you can replace 1 with any positive integer) to choose whether to launch a brand new integration (0) or to load a previously interrupted integration (1).
 If there is no previously saved session (e.g. is the first time you ever run the code) the execution will stop resulting in a segmentation error.
  
 The output produced will be 6 plots of the phase space evenly distributed in time (0%, 20%..., 100%). You can find the produced outputs inside the "plots" folder. In addition at the same integration time that plots are produced it will be produced the save file to resume the integration in a second moment.
 
 In the save file there will be every point's property such as position, velocity, acceleration, mass. Then the save file will also store the integration time and the time step as calculated at save time. Note that the implementation of the leap-frog method (in the kick kick drift form) is such that it is correct to save the velocity in the i-step and not in the i/2-step.

What to expect?
 As of today this program, if correctly compiled, should produce:
  -  a "DISTRRHOB.dat" binary file containing the distributed particle inside the box.*
  -  a "DISTRRHO" text file containing the same as above*
  -  5 ("phaseSpace1",..., "phaseSpace5") files containing a screenshot of the phase space of the 4096 particles. Plotting one over the other (using gnuplot) allows us to appreciate the evolution of the system.

**these files creation is to be set in the param.txt file. Default is to create both .dat and the text file.*

