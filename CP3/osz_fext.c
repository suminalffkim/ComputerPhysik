/*************************************
*kompiliere mit
*gcc -Wall -pedantic osz_fext.c -o fext -lm -L./ -lh3
*aufrufen mit
*./fext
**********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "h3.h"
  
/*****************************************
*DGL-Funktion
********************************/
double f(double t, double x, double xdot, double omega_0, double gamma, double omega) {
    return (-omega_0 * omega_0 * x - 2 * gamma * xdot+ fext(t));
}

/*****************************************
*Runge-Kutta-Verfahren
*******************************************/
double* RK4(double v_0, double omega_0, double gamma, double omega, double t_0, double x_0, double h, unsigned int n) {
    double x = x_0; // Anfangswert für x(t)
    double xdot = v_0; // Anfangswert für xdot(t)
    double j_1, j_2, j_3, j_4, k_1, k_2, k_3, k_4, t;
    double* x_res = malloc(n * sizeof(double)); // x_n wird gespeichert.
    
    x_res[0]=x;

    for (int k = 1; k < n; k++) {
        t = t_0 + k * h;
        j_1 = h * f(t, x, xdot, omega_0, gamma, omega);
        k_1 = h * xdot;
        j_2 = h * f(t + h / 2, x + k_1 * 0.5, xdot + j_1 * 0.5, omega_0, gamma,omega);
        k_2 = h * (xdot + j_1 * 0.5);
        j_3 = h * f(t + h / 2, x + k_2 * 0.5, xdot + j_2 * 0.5, omega_0, gamma, omega);
        k_3 = h * (xdot + j_2 * 0.5);
        j_4 = h * f(t + h, x + k_3, xdot + j_3, omega_0, gamma, omega);
        k_4 = h * (xdot + j_3);

        xdot += (1. / 6.) * (j_1 + 2 * j_2 + 2 * j_3 + j_4); // neue xdot, xdot_(k+1)
        x += (1. / 6.) * (k_1 + 2 * k_2 + 2 * k_3 + k_4); // neue x, x_(k+1)
        x_res[k] = x; // neue x wird in array gespeichert.
    }

    return x_res;
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
    double omega=1.5;
    double t_0=0.;
    unsigned int n=100000;
    double h=0.001; //drei h Werte werden eingesetzt und Ergebnisse vergliechen.
    
 
        printf("File open\n");//zur kontrolle
    
        FILE *fp = fopen("Runge_Kutta_fext.txt", "w");


        double* result=RK4(v_0,omega_0,gamma,omega,t_0,x_0,h,n); //Ergebnisse die mit Runge Kutta Verfahren berechnet wurden
    
        //double* exact=schwingung(omega_0,gamma,alpha,omega,x_0,t_0,h[ih],n[ih]); //Ergebnisse, die mit exakte Loesung berechnet wurden.
    
        
        //in Datei wird n, t, x_n, exact_n, abs(x_n-exakt_n) Werte gespeichert.

        fprintf(fp,"n\tt\t\t\t\t\t\t\t\tx_n\n");
        
        for(int i=0;i<n;i++){
            double t=t_0+h*i;
            
            fprintf(fp, "%d\t%25.19e\t%25.19e\n",i,t,result[i]);
        }

        fclose(fp);
        printf("File close\n");        
    
    
        free(result);
    
        
    
    
}
