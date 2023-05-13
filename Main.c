//main

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allvars.h"
#include "funzioni.h"
#include "MainFunctions.h"
#include <fftw3.h>

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
    
    printf("AAAAAAAAAAAAAAAAA");

    fftw_complex *kDensity, *kPot;
    double *Pot, k, norm;
    Pot = (double*)malloc(N_grid * sizeof(double));
    fftw_plan fft_real_fwd, fft_real_bck;
    kDensity = fftw_malloc( N_grid * sizeof(fftw_complex));
    kPot = fftw_malloc( N_grid * sizeof(fftw_complex));
        /*definition of back and fourth FFT*/
    fft_real_fwd = fftw_plan_dft_r2c_1d(N_grid, massGrid, kDensity, FFTW_ESTIMATE );
    fft_real_bck = fftw_plan_dft_c2r_1d(N_grid, kPot, Pot, FFTW_ESTIMATE );
    //forward
    fftw_execute(fft_real_fwd);
    //divido per k^2
    norm = 2*M_PI/BoxLenght;
    for(i = 1; i<N_grid/2+1; i++){
        k = (i*1.0) * norm;
        kPot[i][0] = -kDensity[i][0]/k/k;
        kPot[i][1] = -kDensity[i][1]/k/k;
    }
    //backward
    fftw_execute(fft_real_bck);
    //renorm
    norm = 1.0 / N_grid;
    FILE *mmm;
    mmm = fopen("pot.txt", "w+");
    for(i=0; i<N_grid; i++){
        Pot[i] *= norm;
        fprintf(mmm, "%lf\n", Pot[i]);
    }
    
    }
