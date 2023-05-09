//main

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allvars.h"
#include "funzioni.h"
#include "MainFunctions.h"

double (*f)(double, double);

int main(){
    FILE *ptrToIC; //puntatore al file delle condizioni iniziali
    FILE *stampaMassPerGrid; //puntatore a un file dove stampo i valori della massa sulla griglia
    int i, j;
    ptrToIC = fopen("DISTRRHOB.dat", "rb");
    char paramFile[] = "param.txt"; //nome parameter file
    //leggi i parametri
    readParamDouble(paramFile);
    readParamInt(paramFile);
    double cellSize, buff[4], massGrid[N_grid], temp;
    for(i=0; i<N_grid; i++){ //Inizializzo la massa sulla griglia a 0 per sicurtade
        massGrid[i]=0;
    }
    cellSize = BoxLenght/N_grid; //grandezza della cella della griglia
    particle_data particle[N_points]; //inizializzo l'array di struct particle_data
    //riempio le struct con i dati prodotti dal file di condizioni iniziali
    for(i=0; i<N_points; i++){
        fread(buff, sizeof(double), 4, ptrToIC);
        particle[i].x=buff[0];
        particle[i].v=buff[1];
        particle[i].a=buff[2];
        particle[i].m=buff[3];
        //printf("%lf %lf %lf %lf %d %d\n", particle[i].x, particle[i].v, particle[i].a, particle[i].m, i, N_points);
    }
    fclose(ptrToIC);
    f=&TSC; //scelta della funzione di distribuzione da utilizzare, da mettere con un ifdef probabilmente
    /*CALCOLO DELLA MASSA SU GRIGLIA*/
    for(j=0; j<N_grid; j++){
        //printf("\n%d\n", j);
        for(i=0; i<N_points; i++){
            massGrid[j] = massGrid[j] + densityJthCell(j, N_grid, N_points, f, BoxLenght, particle[i].m, particle[i].x, cellSize); 
            if(particle[i].x < 3*cellSize){ //condizioni di periodicità
                temp = particle[i].x+BoxLenght;
                massGrid[j] = massGrid[j] + densityJthCell(j, N_grid, N_points, f, BoxLenght, particle[i].m, temp, cellSize);
            }
            if(particle[i].x > BoxLenght - 3*cellSize){ //condizioni di periodicità
                temp = particle[i].x-BoxLenght;
                massGrid[j] = massGrid[j] + densityJthCell(j, N_grid, N_points, f, BoxLenght, particle[i].m, temp, cellSize);
            }
    }}

    stampaMassPerGrid = fopen("Bugas", "w");
    for(i=0; i<N_grid; i++){
        printf("%lf\n", massGrid[i]);
        fprintf(stampaMassPerGrid, "%lf\n", massGrid[i]);
    }
    fclose(stampaMassPerGrid);
    }
