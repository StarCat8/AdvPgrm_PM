extern double unit_vel;
extern double unit_mass;
extern double unit_lenght;
extern int N_points;
extern int N_grid;
extern double BoxLenght;
extern double A_deltaPar;
extern double H_0;
typedef struct{
    double x;
    double v;
    double a;
    double m;
} particle_data;
extern double *denspot;
extern double *force;
extern double Time, deltaT;
extern double allMass;
extern int outText;
extern int outBin;
extern double unit_Time;
extern double second_in_year;
extern double codeUnit_G;
extern double StopTime;
//extern char *Dist_Name;