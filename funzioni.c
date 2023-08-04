//rinomina funzioni_IC, metti massa, fai parameter file, numpunti in rejection deve generare 4096 punti non provarne 4096, metti rhocrit in define
#include "funzioni.h"
#include "allvars.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793
#define e 2.718281828459045
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum){
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    if (*idum < 0 || iff == 0)
        {
        iff=1;
        mj=labs(MSEED-labs(*idum));
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) 
            {
            ii=(21*i) % 55; 
            ma[ii]=mk;       
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
            }
        for (k=1;k<=4;k++)  
            for (i=1;i<=55;i++) 
                {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
                }
        inext=0;
        inextp=31;
        *idum=1;
        }
    if (++inext == 56) inext=1;   
    if (++inextp == 56) inextp=1; 
    mj=ma[inext]-ma[inextp];      
    if (mj < MZ) mj += MBIG;      
    ma[inext]=mj;                 
    return mj*FAC;                
    }

double delta(double A, double x, double L){
    double y;
    y = A*sin(2*PI/L*x-PI/2);
    return y;
    }

double rho(double rho_crit, double A, double x, double L){
    double z;
    z = delta(A, x, L)*rho_crit + rho_crit;
    return z;
    }

void makeHist(char *name, int Nbin){ //crea istogramma da file formattato
    FILE *Hp, *ptr;
    int *HISTO = malloc(sizeof(int) * Nbin);
    int index, i;
    i=0;
    while(i<Nbin){
        HISTO[i]=0;
        i++;
    }
    double MIN, MAX, buf;
    ptr = fopen(name, "r");
    fscanf(ptr, "%*f %lf %*f %*f %*d", &buf);
    MIN = 1.e33;
    MAX = -1.e33;
    //trova minimo e massimo
    while(fscanf(ptr, "%*f %lf %*f %*f %*d", &buf) == 1){
        if(buf<MIN){
            MIN=buf;
        }
        if(buf>MAX){
            MAX=buf;
        }
    }
    fclose(ptr);
    ptr = fopen(name, "r");
    while(fscanf(ptr, "%*f %lf %*f %*f %*d\n", &buf) == 1){
        buf = (buf - MIN)/MAX * Nbin;
        index = (int)buf;
        HISTO[index]=HISTO[index]+1;
    }
    i=0;
    Hp = fopen("Hist.txt", "w+");
    while(i<Nbin){
        fprintf(Hp, "%lf %d\n", i*MAX/Nbin + MIN, HISTO[i]);
        i=i+1;
    }
    fclose(Hp);
    fclose(ptr);
    }
void readParamDouble(char *name){
    char buf[256];
    double x;
    FILE *ptrToParam;
    ptrToParam = fopen(name, "r"); 
    //printf("%s", name);
    while(fscanf(ptrToParam, "%s %lf\n", buf, &x)!=EOF){
        if(!strcmp(buf, "unit_vel")){
            unit_vel = x;
        }
        if(!strcmp(buf, "unit_mass")){
            unit_mass = x;
        }
        if(!strcmp(buf, "unit_lenght")){
            unit_lenght = x;
        }
        if(!strcmp(buf, "BoxLenght")){
            BoxLenght = x;
        }
        if(!strcmp(buf, "A_deltaPar")){
            A_deltaPar = x;
        }
        if(!strcmp(buf, "H_0")){
            H_0 = x;
        }
        if(!strcmp(buf, "StopTime")){
            StopTime = x;
        }
        }
        fclose(ptrToParam);
    }
void readParamInt(char *name){
    char buf[256];
    int x;
    FILE *ptrToParam;
    ptrToParam = fopen(name, "r");  
    while(fscanf(ptrToParam, "%s %d\n", buf, &x)!=EOF){
        if(!strcmp(buf, "N_points")){
            N_points = x;
        }
        if(!strcmp(buf, "N_grid")){
            N_grid = x;
        }
        if(!strcmp(buf, "outText")){
            outText = x;
        }
        if(!strcmp(buf, "outBin")){
            outBin = x;
        }
    }
    fclose(ptrToParam);
    }

/*void readParamChar(char *name){
    char buf[256];
    char *x;
    FILE *ptrToParam;
    ptrToParam = fopen(name, "r");  
    while(fscanf(ptrToParam, "%s %s\n", buf, x)!=EOF){
        if(!strcmp(buf, "Dist_Name")){
            Dist_Name = x;
    }
    fclose(ptrToParam);
    }
    }
*/

void rejection(double ymin, double ymax, int numpunti, double intervallo, long *dumm, double rho_crit, double A, double massa, double i_speed, char *name, int TEXT, int BINAR){
	FILE *fileT, *fileB;
    double puntiBin[5];
	fileT = fopen(name, "w");
    char bina[6] = "B.dat";
    strcat(name, bina);
    fileB = fopen(name, "wb");
	double xtry, ytry, y; //
	int i = 0;
	while(i < numpunti){ //
		xtry = 1.*(intervallo * ran3(dumm));
		ytry = 1.*((ymax-ymin)*ran3(dumm));
		y = rho(rho_crit, A, xtry, intervallo); //
        if(ytry < y){
            if (TEXT){ //
                fprintf(fileT, "%lf %lf %lf %lf\n", xtry, i_speed, 0.0, massa);
            }
            if (BINAR){ //
                puntiBin[0] = xtry;
                puntiBin[1] = i_speed;
                puntiBin[2] = 0.0;
                puntiBin[3] = massa;
                fwrite(puntiBin, sizeof(double), 4, fileB); //f_try(x_0) < f(x_0) accepted, x, f(x_0), numero del punto
                /*stampa su file x FOPEN(ecc)*/
                i++;
            }
        }
	}
	fclose(fileT);
    fclose(fileB);
    }

 //Funzioni non utilizzate

/*double * creaPunti(int Npnt, double (*f)(double), long int *dumm){
    int i;
    double w;
    double *z=malloc(sizeof(double) * Npnt);
    for(i=0; i<Npnt; i++){
        z[i] = f(ran3(dumm));
        //printf("\n%lf", z[i]);
    }
    return z;
    }

double anglePHI(double phi){
    phi=2*PI*phi;
    return phi;
    }
double angleTHETA(double costheta){
    double theta;
    costheta=2*costheta-1;
    theta = acos(costheta);
    return theta;
    }
double radialC(double x){
    x = pow(x, 0.5);
    return x;
    }

double radialS(double x){
    x = pow(x, (1.*1/3));
    return x;
    }


double areaCerchio(double r){
    double Area = PI*r*r;
    return Area;
    }

double volSfera(double r){
        double volume = 4/3*PI*r*r*r;
        return volume;
    }

double density(double *punti, double (*A)(double), double dist1, double dist2, int Npnt){
    int i, j=0;
    double a, rho;
    a=(A(dist2)-A(dist1));
    for(i=0; i<Npnt; i++){
            if(*(punti + i)>=dist1 && *(punti + i)<=dist2){
                j++;
            }
        rho = j/a;
        }
    return rho;
    }

    */          
