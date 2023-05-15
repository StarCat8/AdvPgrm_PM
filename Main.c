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
    int i, j, MAX_STEPS;
    MAX_STEPS=10000; //NUMERO DI VOLTE CHE SI INTEGRA
    double *Pot, k, norm;
    ptrToIC = fopen("DISTRRHOB.dat", "rb");
    char paramFile[] = "param.txt"; //nome parameter file
    //leggi i parametri
    readParamDouble(paramFile); //
    readParamInt(paramFile);
    double cellSize, buff[4], massGrid[N_grid], temp;
    cellSize = BoxLenght/N_grid; //grandezza della cella della griglia
    particle_data particle[N_points]; //inizializzo l'array di struct particle_data
    //riempio le struct con i dati prodotti dal file di condizioni iniziali
    
    f=&TSC; //scelta della funzione di distribuzione da utilizzare, da mettere con un ifdef probabilmente
    /*inizializzazione particelle in struct*/
    for(i=0; i<N_points; i++){
        fread(buff, sizeof(double), 4, ptrToIC);
        particle[i].x=buff[0];
        particle[i].v=buff[1];
        particle[i].a=buff[2];
        particle[i].m=buff[3];
        }
    fclose(ptrToIC);

                  /*FFTW*/
    fftw_complex *kDensity, *kPot;
    Pot = (double*)malloc(N_grid * sizeof(double));
    kDensity = fftw_malloc( N_grid * sizeof(fftw_complex));
    kPot = fftw_malloc( N_grid * sizeof(fftw_complex));
    fftw_plan fft_real_fwd, fft_real_bck;
    fft_real_fwd = fftw_plan_dft_r2c_1d(N_grid, massGrid, kDensity, FFTW_ESTIMATE );
    fft_real_bck = fftw_plan_dft_c2r_1d(N_grid, kPot, Pot, FFTW_ESTIMATE );

                 /*force*/
    double *force, *forceOnParticle;
    force = (double*)malloc(N_grid * sizeof(double));
    forceOnParticle = (double*)malloc(N_grid * sizeof(double));
                    
                    
                    /////* MAIN INTEGRATION LOOP */////
    for(int leapCounter=0; leapCounter<MAX_STEPS; leapCounter++){
                        /*inizializzazione massa*/
        for(i=0; i<N_grid; i++){ //Inizializzo la massa sulla griglia a 0 per sicurtade
            massGrid[i]=0;
        }
        /*CALCOLO DELLA DENSITA' SU GRIGLIA*/
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

                    //FFTW

            /*definition of back and fourth FFT*/
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
        /*force*/
        for(int p = 0; p < N_points ; p++) //loop per calcolare la forza su singola particella p data dalla distribuzione secondo kernel TSC delle forza i sulla griglia.
        {
            force[0] = -(Pot[0] - Pot[N_grid-1]) / (2.0*cellSize);
            forceOnParticle[p] = f( particle[p].x - 0.5*cellSize, cellSize)*force[0];
            for(i=1; i<N_grid-1; i++)
            {
                force[i] = -(Pot[i+1]- Pot[i-1])/ (2.0*cellSize);
                forceOnParticle[p] += f(particle[p].x-i*cellSize - 0.5*cellSize, cellSize)*force[i];
            }
            force[N_grid] = -(Pot[N_grid] - Pot[0]) / (2.0*cellSize);
            forceOnParticle[p] += f( particle[p].x - N_grid*cellSize - 0.5*cellSize, cellSize)*force[N_grid];
            particle[p].a=forceOnParticle[p]/particle[p].m;
            //printf("\n<%lf>", particle[p].a);
        }
                    /*LEAP FROG*/
        double step; 
        double tmp = 0;
        //calcolo timeStep <--- CREDO CI SIANO ERRORI
        int boole = 0;
        //if(leapCounter != 0 ){
        if(boole == 1 ){ // <--- MOTIVO PER CUI IMPONGO UNO STEP FISSO
        for(i = 0; i<N_points; i++){
            if(tmp < fabs(particle[i].v)){
            tmp = fabs(particle[i].v);
            }
        }
        step = timeStep(cellSize, tmp); //STEP FISSO
        } else {
            step = cellSize/1000;
        }
        //update particles' position and velocity
        for(i = 0; i<N_points; i++){
            particle[i].v = halfStepVelocity(particle[i].v, step, forceOnParticle[i], particle[i].m);
            particle[i].x = posUpdate(particle[i].x, particle[i].v, step);
            if(particle[i].x < 0){                            //Boundary cond
                particle[i].x = particle[i].x + BoxLenght;
            }
            if(particle[i].x > BoxLenght){                    //Boundary cond
                particle[i].x = particle[i].x - BoxLenght;
            }
        }

        if(leapCounter==0){ //Stampo su file la traiettoria nello spazio delle fasi (x vs v) di 5 particelle
            FILE *phaseSpace1;
            phaseSpace1 = fopen("phaseSpace1", "w");
            for(i=0; i<N_points; i++){
            fprintf(phaseSpace1, "%lf %lf\n", particle[i].x, particle[i].v);
            }
            fclose(phaseSpace1);
        }
        if(leapCounter==40){
            FILE *phaseSpace2;
            phaseSpace2 = fopen("phaseSpace2", "w");
                    for(i=0; i<N_points; i++){
            fprintf(phaseSpace2, "%lf %lf\n", particle[i].x, particle[i].v);
            }
            fclose(phaseSpace2);
        }
        if(leapCounter==80){
            printf("EHEH");
            FILE *phaseSpace3;
            phaseSpace3 = fopen("phaseSpace3", "w");
                    for(i=0; i<N_points; i++){
            fprintf(phaseSpace3, "%lf %lf\n", particle[i].x, particle[i].v);
            }
            fclose(phaseSpace3);
        }
        if(leapCounter==120){
            FILE *phaseSpace4;
            phaseSpace4 = fopen("phaseSpace4", "w");
            for(i=0; i<N_points; i++){
            fprintf(phaseSpace4, "%lf %lf\n", particle[i].x, particle[i].v);
            }
            fclose(phaseSpace4);
        }
        if(leapCounter==160){
            FILE *phaseSpace5;
            phaseSpace5 = fopen("phaseSpace5", "w");
            for(i=0; i<N_points; i++){
            fprintf(phaseSpace5, "%lf %lf\n", particle[i].x, particle[i].v);
            }
            fclose(phaseSpace5);
        }
        if(leapCounter%100==0)
        printf("%d/100\n", leapCounter/(MAX_STEPS/100));

        //printf("\nleapCount: %lf %lf %d \n", particle[0].x, particle[0].v, leapCounter);
        }
        //free(forceOnParticle);
        //fftw_free(kDensity);
        //fftw_free(kPot);
        //fftw_destroy_plan(fft_real_fwd);
        //fftw_destroy_plan(fft_real_bck);
        //free(force);
        //free(Pot);
        /*fclose(phaseSpace1);
        fclose(phaseSpace2);
        fclose(phaseSpace3);
        fclose(phaseSpace4);
        fclose(phaseSpace5);*/
    }
