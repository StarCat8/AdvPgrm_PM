double NGP(double x, double Cell_Size);
double CIC(double x, double Cell_Size);
double TSC(double x, double Cell_Size);
double densityJthCell(int j, double (*f)(double, double), double BoxLenght, double massa, double partPos, double cellSize);
double timeStep(double cellSize, double minVel);
double halfStepVelocity(double v, double timeStep, double F, double mass);
double posUpdate(double x, double v, double timeStep);

