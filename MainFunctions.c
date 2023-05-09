#include "funzioni.h"
#include "allvars.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double NGP(double x, double Cell_Size){ //Nearest Grid Point
    if(fabs(x)<Cell_Size){
        return Cell_Size;
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

double densityJthCell(int j, double N_grid, double N_points, double (*f)(double, double), double BoxLenght, double massa, double partPos, double cellSize){ //in realtà è la massa nella j-esima cella. Questa viene distribuita tramite una funzione di distribuzione sia essa NGP CIC TSC
    double massContr = 0;
    double cellPos = cellSize/2.0+cellSize*j;
    //printf("%lf", partPos-cellPos);
    massContr = massa*f(partPos-cellPos, cellSize);
    return massContr;
}