//generatore di condizioni iniziali

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "funzioni.h"
#include "allvars.h"

int main(){
    double rho_crit, total_mass, particle_mass;
    long int uu;
    long int *dumm = &uu;
    int NumBin=6;
    char nome_dist[13] = "DISTRRHO";
    char paramFile[] = "param.txt";
    double i_speed = 0;
    //leggi i parametri
    readParamDouble(paramFile);
    readParamInt(paramFile);
    //printf(">%d<, <%d>", outBin, outText);
    rho_crit = 2.775*100*H_0*H_0; //in SunMass / Kpc^3
    total_mass = rho_crit*BoxLenght*BoxLenght*BoxLenght;
    particle_mass = total_mass / N_points;
    rejection(0, rho_crit + rho_crit*A_deltaPar, N_points, BoxLenght, dumm, rho_crit, A_deltaPar, particle_mass, i_speed, nome_dist, outText, outBin);
    makeHist(nome_dist, NumBin);
    unit_Time = unit_lenght/unit_vel/(second_in_year*pow(10, 9));
    printf("Unit of time = %lf in Gy\n", unit_Time);
    double cm1= 1/unit_lenght;
    double g1= 1/unit_mass;
    double sec1= 1/(unit_Time*second_in_year*pow(10, 9));
    codeUnit_G = 6.67430*pow(10, -8) / g1 * cm1* cm1* cm1 / sec1/ sec1;
    printf("Gravitational constant = %lf in code unit", codeUnit_G);
    return 0;
}
