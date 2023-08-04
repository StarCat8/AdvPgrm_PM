//main

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allvars.h"
#include "funzioni.h"
#include "MainFunctions.h"
#include <fftw3.h>
#include <sys/time.h>
//#include <omp.h>

double (*f)(double, double);

int main(int argc, char **argv){
    if(strcmp(argv[1],"-NGP")==0)
        {
            f=&NGP;
        }
    if(strcmp(argv[1],"-CIC")==0)
        {
            f=&CIC;
        }
    if(strcmp(argv[1],"-TSC")==0)
        {
            f=&TSC;
        }
    FILE *ptrToIC; //puntatore al file delle condizioni iniziali
    int i, j, MAX_STEPS, leapCounter = 0, boola = 1, boolb = 1, boolc = 1, boold = 1, boolf = 1, boolg = 1; 
    double unit_Time = 0.978462;
    double step;
    MAX_STEPS=5001; //NUMERO DI VOLTE CHE SI INTEGRA
    double *Pot, k, norm, elapsedTime;
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);
    char paramFile[] = "param.txt"; //nome parameter file
    readParamDouble(paramFile); //
    readParamInt(paramFile);
    system("mkdir plots");
    double cellSize, buff[4], massGrid[N_grid], temp, integrationTime, StartTime, buff1[6];

    int newORresume;
    cellSize = BoxLenght/N_grid; //grandezza della cella della griglia
    particle_data particle[N_points]; //inizializzo l'array di struct particle_data
    printf("Si vuole iniziare una nuova simulazione o riprendere una precedentemente interrotta?\n\n (ci si assicuri che esista il file di salvataggio avente nome 'mostRecentSession.dat')\n\n Digitare 0 se si vuole iniziare una nuova integrazione, 1 se si vuole riprendere dall'ultimo salvataggio disponibile:   ");
    scanf("%d", &newORresume);
    if(newORresume==0){ //nuova integrazione o riprendere da una interrotta
    step = StopTime/100;
    ptrToIC = fopen("DISTRRHOB.dat", "rb");
    StartTime = 0;
        for(int i=0; i<N_points; i++)
        {
        fread(buff, sizeof(double), 4, ptrToIC);
        particle[i].x=buff[0];
        particle[i].v=buff[1];
        particle[i].a=buff[2];
        particle[i].m=buff[3];
        }
    //leggi i parametri
    //riempio le struct con i dati prodotti dal file di condizioni iniziali
    
    //scelta della funzione di distribuzione da utilizzare, da mettere con un ifdef probabilmente
    /*inizializzazione particelle in struct*/
    } else {
    ptrToIC = fopen("mostRecentSession.dat", "rb");
    for(int i=0; i<N_points; i++)
        {
        fread(buff1, sizeof(double), 6, ptrToIC);
        particle[i].x=buff1[0];
        particle[i].v=buff1[1];
        particle[i].a=buff1[2];
        particle[i].m=buff1[3];
        StartTime = buff1[4];
        step = buff1[5];
        }
    }
    integrationTime = StartTime;
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
    //#pragma omp parallel for num_threads(4)
    while(integrationTime <= StopTime){
                        /*inizializzazione massa*/
        for(i=0; i<N_grid; i++){ //Inizializzo la massa sulla griglia a 0 per sicurtade
            massGrid[i]=0;
            
        }
        /*CALCOLO DELLA DENSITA' SU GRIGLIA*/
        int temporary;
        double temporar;
        for(i=0; i<N_points; i++){
            temporar = particle[i].x / BoxLenght * N_grid;
            temporary = (int)temporar;
            massGrid[temporary] = massGrid[temporary] + densityJthCell(temporary, f, BoxLenght, particle[i].m, particle[i].x, cellSize); 
            if(particle[i].x < 3*cellSize){ //condizioni di periodicità
                temp = particle[i].x+BoxLenght;
                massGrid[temporary] = massGrid[temporary] + densityJthCell(temporary, f, BoxLenght, particle[i].m, temp, cellSize);
            }
            if(particle[i].x > BoxLenght - 3*cellSize){ //condizioni di periodicità
                temp = particle[i].x-BoxLenght;
                massGrid[temporary] = massGrid[temporary] + densityJthCell(temporary, f, BoxLenght, particle[i].m, temp, cellSize);
            }
        }
        

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
            force[N_grid-1] = -(Pot[N_grid-1] - Pot[0]) / (2.0*cellSize);
            forceOnParticle[p] += f( particle[p].x - N_grid*cellSize - 0.5*cellSize, cellSize)*force[N_grid];
            particle[p].a=forceOnParticle[p]/particle[p].m;
        }
                    /*LEAP FROG*/
        double tmp = 0;
        if(leapCounter > 9 && leapCounter%10 == 0){ // <--- MOTIVO PER CUI IMPONGO UNO STEP FISSO
            for(i = 0; i<N_points; i++){
                if(tmp < fabs(particle[i].v)){
                tmp = fabs(particle[i].v);
                }
            }
        step = timeStep(cellSize, tmp); //STEP FISSO
        }

        //update particles' position and velocity
        for(i = 0; i<N_points; i++){
            particle[i].v = particle[i].v + 0.5*forceOnParticle[i]/particle[i].m * step;
            
            particle[i].v = halfStepVelocity(particle[i].v, step, forceOnParticle[i], particle[i].m);
            particle[i].x = posUpdate(particle[i].x, particle[i].v, step);
            if(particle[i].x < 0){                            //Boundary cond
                particle[i].x = particle[i].x + BoxLenght;
            }
            if(particle[i].x > BoxLenght){                    //Boundary cond
                particle[i].x = particle[i].x - BoxLenght;
            }
        }
        codeUnit_G = 43021.931323;
        FILE *energy;
        energy = fopen("energy", "a");








        double saveBuf[6];
        char saveName[32] = "mostRecentSession.dat";
        if(integrationTime/(StopTime/100)>10 && boola && integrationTime/(StopTime/100)<15){ //Stampo su file la traiettoria nello spazio delle fasi (x vs v) di 5 particelle
            double E_kinetic = 0, E_potential = 0;
            FILE *phaseSpace1, *energy1, *saveF;
            saveF = fopen(saveName, "wb");
            phaseSpace1 = fopen("plots/phaseSpace1", "w");
            for(i = 0; i < N_points; i++){
                fprintf(phaseSpace1, "%lf %lf\n", particle[i].x, particle[i].v);
                E_kinetic += 0.5 * particle[i].m * particle[i].v * particle[i].v;
                saveBuf[0]=particle[i].x;
                saveBuf[1]=particle[i].v;
                saveBuf[2]=particle[i].a;
                saveBuf[3]=particle[i].m;
                saveBuf[4]=integrationTime;
                saveBuf[5]= step;
                fwrite(saveBuf, sizeof(double), 6, saveF);
                for(j = 0; j < N_points; j++){
                    if(i!=j){
                    E_potential += particle[i].m*particle[j].m/(particle[j].x-particle[i].x);
                }
                }
            }
            fclose(phaseSpace1);
            fclose(saveF);
            boola = boola -1;
            E_potential = 0.5 * codeUnit_G * E_potential;
            fprintf(energy, "E_kin = %lf     E_pot = %lf    E_tot = %lf\n", E_kinetic, E_potential, E_kinetic + E_potential);
        }












        if(integrationTime/(StopTime/100)>20 && boolb && integrationTime/(StopTime/100)<22){
            double E_kinetic = 0, E_potential = 0;
            FILE *phaseSpace2, *saveF;
            saveF = fopen(saveName, "wb");
            phaseSpace2 = fopen("plots/phaseSpace2", "w");
                for(i=0; i<N_points; i++){
                    fprintf(phaseSpace2, "%lf %lf\n", particle[i].x, particle[i].v);
                    E_kinetic += 0.5 * particle[i].m * particle[i].v * particle[i].v;
                    saveBuf[0]=particle[i].x;
                    saveBuf[1]=particle[i].v;
                    saveBuf[2]=particle[i].a;
                    saveBuf[3]=particle[i].m;
                    saveBuf[4]=integrationTime;
                    saveBuf[5]= step;
                fwrite(saveBuf, sizeof(double), 6, saveF);
                for(j = 0; j < N_points; j++){
                    if(i!=j){
                    E_potential += particle[i].m*particle[j].m/(particle[j].x-particle[i].x);
                }
                }
            }
            fclose(phaseSpace2);
            fclose(saveF);
            boolb = boolb -1;
            E_potential = 0.5 * codeUnit_G * E_potential;
            fprintf(energy, "E_kin = %lf     E_pot = %lf      E_tot = %lf\n", E_kinetic, E_potential, E_kinetic + E_potential);
        }

        if(integrationTime/(StopTime/100)>40 && boolc && integrationTime/(StopTime/100)<42){
            double E_kinetic = 0, E_potential = 0;
            FILE *phaseSpace3, *saveF;
            saveF = fopen(saveName, "wb");
            phaseSpace3 = fopen("plots/phaseSpace3", "w");
                    for(i=0; i<N_points; i++){
                        fprintf(phaseSpace3, "%lf %lf\n", particle[i].x, particle[i].v);
                        E_kinetic += 0.5 * particle[i].m * particle[i].v * particle[i].v;
                        saveBuf[0]=particle[i].x;
                        saveBuf[1]=particle[i].v;
                        saveBuf[2]=particle[i].a;
                        saveBuf[3]=particle[i].m;
                        saveBuf[4]=integrationTime;
                        saveBuf[5]= step;
                        fwrite(saveBuf, sizeof(double), 6, saveF);
                        for(j = 0; j < N_points; j++){
                            if(i!=j){
                            E_potential += particle[i].m*particle[j].m/(particle[j].x-particle[i].x);
                        }
                }
            }
            fclose(phaseSpace3);
            fclose(saveF);
            boolc = boolc -1;
            E_potential = 0.5 * codeUnit_G * E_potential;
            fprintf(energy, "E_kin = %lf     E_pot = %lf      E_tot = %lf\n", E_kinetic, E_potential, E_kinetic + E_potential);
        }

        if(integrationTime/(StopTime/100)>60 && boold && integrationTime/(StopTime/100)<62){
            double E_kinetic = 0, E_potential = 0;
            FILE *phaseSpace4, *saveF;
            saveF = fopen(saveName, "wb");
            phaseSpace4 = fopen("plots/phaseSpace4", "w");
            for(i=0; i<N_points; i++){
            fprintf(phaseSpace4, "%lf %lf\n", particle[i].x, particle[i].v);
            E_kinetic += 0.5 * particle[i].m * particle[i].v * particle[i].v;
                        saveBuf[0]=particle[i].x;
                        saveBuf[1]=particle[i].v;
                        saveBuf[2]=particle[i].a;
                        saveBuf[3]=particle[i].m;
                        saveBuf[4]=integrationTime;
                        saveBuf[5]= step;
                        fwrite(saveBuf, sizeof(double), 6, saveF);
                for(j = 0; j < N_points; j++){
                    if(i!=j){
                    E_potential += particle[i].m*particle[j].m/(particle[j].x-particle[i].x);
                }
                }
            }
            fclose(phaseSpace4);
            fclose(saveF);
            boold = boold -1;
            E_potential = 0.5 * codeUnit_G * E_potential;
            fprintf(energy, "E_kin = %lf     E_pot = %lf      E_tot = %lf\n", E_kinetic, E_potential, E_kinetic + E_potential);
        }

        if(integrationTime/(StopTime/100)>80 && boolf && integrationTime/(StopTime/100)<81){
            double E_kinetic = 0, E_potential = 0;
            FILE *phaseSpace5, *saveF;
            saveF = fopen(saveName, "wb");
            phaseSpace5 = fopen("plots/phaseSpace5", "w");
            for(i=0; i<N_points; i++){
            fprintf(phaseSpace5, "%lf %lf\n", particle[i].x, particle[i].v);
            E_kinetic += 0.5 * particle[i].m * particle[i].v * particle[i].v;
                        saveBuf[0]=particle[i].x;
                        saveBuf[1]=particle[i].v;
                        saveBuf[2]=particle[i].a;
                        saveBuf[3]=particle[i].m;
                        saveBuf[4]=integrationTime;
                        saveBuf[5]= step;
                        fwrite(saveBuf, sizeof(double), 6, saveF);
                for(j = 0; j < N_points; j++){
                    if(i!=j){
                    E_potential += particle[i].m*particle[j].m/(particle[j].x-particle[i].x);
                }
                }
            }
            fclose(phaseSpace5);
            fclose(saveF);
            boolf = boolf -1;
            E_potential = 0.5 * codeUnit_G * E_potential;
            fprintf(energy, "E_kin = %lf     E_pot = %lf      E_tot = %lf\n", E_kinetic, E_potential, E_kinetic + E_potential);
        }
        if( integrationTime/(StopTime/100)>98 && boolg && integrationTime/(StopTime/100)<99){
            double E_kinetic = 0, E_potential = 0;
            FILE *phaseSpace6, *saveF;
            saveF = fopen(saveName, "wb");
            phaseSpace6 = fopen("plots/phaseSpace6", "w");
            for(i=0; i<N_points; i++){
                        fprintf(phaseSpace6, "%lf %lf\n", particle[i].x, particle[i].v);
                        E_kinetic += 0.5 * particle[i].m * particle[i].v * particle[i].v;
                        saveBuf[0]=particle[i].x;
                        saveBuf[1]=particle[i].v;
                        saveBuf[2]=particle[i].a;
                        saveBuf[3]=particle[i].m;
                        saveBuf[4]=integrationTime;
                        saveBuf[5]= step;
                        fwrite(saveBuf, sizeof(double), 6, saveF);
                        for(j = 0; j < N_points; j++){
                            if(i!=j)
                            {
                                E_potential += particle[i].m*particle[j].m/(particle[j].x-particle[i].x);
                            }
                        }
            }
            boolg = boolg -1;
            fclose(phaseSpace6);
            fclose(saveF);
            E_potential = 0.5 * codeUnit_G * E_potential;
            fprintf(energy, "E_kin = %lf     E_pot = %lf      E_tot = %lf\n", E_kinetic, E_potential, E_kinetic + E_potential);
        }
        if(leapCounter%10==0){
        printf("%.1lf/100 step: %lf\n", integrationTime/(StopTime/100), step);
        }
        leapCounter++;
        integrationTime += unit_Time  * step;
        fclose(energy);
        }
    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec);      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    printf("\n---------------------------------------------\n|Integration completed in %lf seconds|\n---------------------------------------------\n\n", elapsedTime);
        return 0;
    }

//Fare plot dell'energia in funzione del tempo cinetica e potenziale e totale.
//Fare energia potenziale vs tempo, calcolarla in 2 modi, uno utilizzando il potenziale datoci dalle fftw
// e calcolare come in fisica 1.
