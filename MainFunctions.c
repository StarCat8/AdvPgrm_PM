#include "funzioni.h"
#include "allvars.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double NGP(double x, double Cell_Size){ //Nearest Grid Point
    if(fabs(x)<Cell_Size){
        return 1;
    } else {
        return 0;
    }
}
double CIC(double x, double Cell_Size){ //Cloud In Cell
    double k;
    if(fabs(x)<Cell_Size){
        k=Cell_Size-fabs(x);
        return k/Cell_Size;
    } else {
        return 0;
    }
}
double TSC(double x, double Cell_Size){ //Triangular Shaped Cloud
    double k=0;
    if(fabs(x)<Cell_Size/2.0){
        k = Cell_Size*Cell_Size*0.75 - x*x;
        //printf("u %lf u", k);
        return k/(Cell_Size*Cell_Size);
    }
    if(Cell_Size*0.5<fabs(x) && fabs(x)< 1.5 * Cell_Size){
        k = 0.5*(Cell_Size*1.5-fabs(x))*(Cell_Size*1.5-fabs(x));
        //printf("a %lf a", k);
        return k/(Cell_Size*Cell_Size);
    } else {
        //printf("i");
        return 0;
    }
}

double densityJthCell(int j, double (*f)(double, double), double BoxLenght, double massa, double partPos, double cellSize){ //in realtà è la massa nella j-esima cella. Questa viene distribuita tramite una funzione di distribuzione sia essa NGP CIC TSC
    double retVal = 0;
    double cellPos = 0.5 * cellSize + cellSize*j;
    //printf("%lf", partPos-cellPos);
    retVal = massa*f(partPos-cellPos, cellSize)/cellSize/BoxLenght/BoxLenght;
    return retVal;
}

double timeStep(double cellSize, double maxVel){ //Calcolo del time step. Probabilmente contiene un errore, di fatti non è stata usata (è commentata nel main)
    double w, retVal;
    w = 0.3;
    retVal = w * cellSize/maxVel;
    return retVal;
}

double halfStepVelocity(double v, double timeStep, double F, double mass){ //update della velocità nel leap frog
    double retVal;
    retVal = v + F/mass * 0.5 * timeStep;
    return retVal;
}
double posUpdate(double x, double v, double timeStep){ //Update della posizione nel leap frog
    double retVal;
    retVal = x + v*timeStep;
    return retVal;
}

/*void saveProgress(particle_data *part, int numberOfPoints){
    FILE *save;
    double buff[4];
    save = fopen("save.txt", "wb");
    for(int i = 0; i<numberOfPoints; i++){
        buff[0]=part[i].x;
        buff[1]=part[i].v;
        buff[2]=part[i].a;
        buff[3]=part[i].m;
        fprintf(buff, sizeof(double), 4, save);
    }
}*/