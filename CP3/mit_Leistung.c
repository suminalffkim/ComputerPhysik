/*************************************
*kompiliere mit
*gcc -Wall -pedantic mit_Les.c -o leistung -lm -L./ -lh3
*aufrufen mit
*./les
**********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "h3.h"
  
/*****************************************
*DGL-Funktion
********************************/
double f(double t, double x, double v) {
    return v;
}

double g(double t, double x, double v,double omega_0, double gamma) {
   return (-omega_0 * omega_0 * x - 2 * gamma * v+ fext(t));
}

/****************************************************************
*Runge-Kutta-Verfahren 4.Ordnung um xdot(=v) zu bestimmen und
*momentane Leistung P=fext(t)*xdot(t) wird berechnet
*********************************************************************/
double* leistung(double v, double omega_0, double gamma, double t_0,double x, double h, unsigned int n){
    double k1v, k2v, k3v, k4v;
    double k1x, k2x, k3x, k4x;
    double t;
    double* P_res=malloc(n*sizeof(double));
    
    P_res[0]=fext(0)*v;
    
    for (int i = 1; i < n; i++) {
        t=t_0+i*h;
        
        k1v = h * g(t, x, v,omega_0,gamma);
        k1x = h * f(t, x, v);

        k2v = h * g(t + h/2, x + k1x/2, v + k1v/2,omega_0,gamma);
        k2x = h * f(t + h/2, x + k1x/2, v + k1v/2);

        k3v = h * g(t + h/2, x + k2x/2, v + k2v/2,omega_0,gamma);
        k3x = h * f(t + h/2, x + k2x/2, v + k2v/2);

        k4v = h * g(t + h, x + k3x, v + k3v,omega_0,gamma);
        k4x = h * f(t + h, x + k3x, v + k3v);

        x +=(k1x + 2*k2x + 2*k3x + k4x) / 6.;
        v +=(k1v + 2*k2v + 2*k3v + k4v) / 6.;

        P_res[i]=fext(t)*v;
    }
      return P_res;
}


/********************************************************************
*Die Funktion findet mittlere Leistung einzelne Periode 
**************************************************************************/
double* mitt_leistung(double* mom_leistung ,unsigned int n) {
    int counter = 1;
    int index = 0;
    int* max_ind=malloc(13*sizeof(int));
    double sum=0.;
    double* mitt_leistung=malloc(12*sizeof(double));

       //die Maximumstelle werden gespeichert.
int j=0;
    for (int i = 0; i < n - 1 && index < 13; i++) {
        switch (counter % 2) {
            case 1:
                
                
                while(mom_leistung[i+1]>mom_leistung[i]){ //Maximum erreicht
                    sum+=mom_leistung[i];
                    j++;
                }
                mitt_leistung[index]=sum/j;
                counter++;
                j=0; 
                sum=0;

                break;
                
            case 0:
                while(mom_leistung[i+1]<=mom_leistung[i]){// Minimum erreicht
                    sum+=mom_leistung[i];
                    j++;
                }
                counter++;
                break;
        }
    
    }
    
    return mitt_leistung;
}
    
/*******************************************
* Die Funktion berechnet die exakten Werte
********************************************/

double* schwingung(double omega_0, double gamma, double alpha, double omega, double x_0, double t_0, double h, unsigned int n) {
    double t, Q;
    double* schwing_res = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        t = t_0 + i * h;
        Q=omega_0/(2*gamma);
        schwing_res[i]=alpha*cos(omega*t-alpha)/(sqrt(pow(omega_0*omega_0-omega*omega,2)+(omega_0*omega_0*omega*omega/(Q*Q))));
    }
        return schwing_res;
    }



/*******************************
*main
*********************************/
int main(int argc, char **argv){
    double v_0=2.75;
    double omega_0=1.;
    double x_0=1.5;
    double gamma=0.1;
    double t_0=0.;
    int n=100000; //n Schritte werden berechnet
    double h=0.001; //h Werte
    
     double* mom_leistung=leistung(v_0,omega_0,gamma,t_0,x_0,h,n); //momentane Leistung werden gespeichert.

    
    printf("File open\n");//zur kontrolle
   FILE *f = fopen("momentane_leistung.txt", "w");
    
        
    //in Datei wird n, t, P_n Werte gespeichert.
    
    //FILE *f = fopen("mittlere_leistung.txt", "w");
   // double* mitt_res= mitt_leistung(mom_leistung,n);
      
    
    //in Datei wird n, t, P_n Werte gespeichert.
    fprintf(f,"n\tt\t\t\t\t\t\t\t\tP_n\n");
    for(int i=0;i<n;i++){
        double t=t_0+h*i;
        fprintf(f, "%d\t%25.19e\t%25.19e\n",i,t,mom_leistung[i]);
        }

        
        
   /* for(int i=0;i<12;i++){
        double t=t_0+h*i;
        fprintf(f, "%d\t%25.19e\t%25.19e\n",i,t,mitt_res[i]);
        }
*/
        
    fclose(f);
        
    printf("File close\n");
    
    //Differenze pro Periode 
    
    
    free(mom_leistung);
    //free(mitt_res);
        
        
    
    
}
