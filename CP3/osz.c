/*************************************
*kompiliere mit
*gcc -O2 -Wall -pedantic osz.c -o osz -lm
*aufrufen mit
*./osz
**********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
  
/*****************************************
*DGL-Funktion
********************************/
double f(double t_0, double x, double xdot,double omega_0, double gamma){    
    return (-omega_0*omega_0*x-2*gamma*xdot);
}

/*****************************************
*Runge-Kutta-Verfahren
*******************************************/
double* RK4(double v_0, double omega_0, double gamma,double t_0, double x_0,  double h,unsigned int n){
    double x=x_0; //Anfangswert fuer x(t)
    double xdot=v_0; //Anfangswert fuer xdot(t)
    double j_1,j_2,j_3,j_4,k_1,k_2,k_3,k_4,t;
    double* x_res=malloc(n*sizeof(double));//x_n wird gespeichert.

    
        for(int k=0; k<n; k++){
            t=t_0+k*h;
            j_1=h*f(t,x,xdot,omega_0,gamma);
            k_1=h*xdot;
            j_2=h*f(t+h/2,x+k_1*0.5, xdot+j_1*0.5,omega_0,gamma);
            k_2=h*(xdot+j_1*0.5);
            j_3=h*f(t+h/2,x+k_2*0.5, xdot+j_2*0.5,omega_0,gamma);
            k_3=h*(xdot+j_2*0.5);
            j_4=h*f(t+h,x+k_3,xdot+j_3,omega_0,gamma);
            k_4=h*(xdot+j_3);
            
            xdot+=(1./6.)*(j_1+2*j_2+2*j_3+j_4); //neue xdot, xdot_(k+1)
            x+=(1./6.)*(k_1+2*k_2+2*k_3+k_4); //neue x, x_(k+1)
            x_res[k]=x; //neue x wird in array gespeichert.
            //printf("%25.19e\n",j_4);
    }
    
    return x_res; 
}

/*******************************************
*Die Funktion wird exact Wert ermitteln
************************************************/
double* schwingfall(double omega_0, double gamma,double x_0, double t_0, double h, unsigned int n){
    double t;
    double *schwing_res=malloc(n*sizeof(double));
    
    for(int i=0;i<n;i++){
        t=t_0+i*h;
        schwing_res[i]=x_0*exp(-gamma*t)*cos(omega_0*t);
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
    unsigned int n=60;
    double h=1.;
    
    
    double* result=RK4(v_0,omega_0,gamma,t_0,x_0,h,n); //Ergebnisse die mit Runge Kutta Verfahren berechnet wurden
    double* exact=schwingfall(omega_0,gamma,x_0,t_0,h,n); //Ergebnisse, die mit exakte Loesung berechnet wurden.
    
    
    FILE *fp = fopen("Runge_Kutta.txt", "w"); // File offenen 
    
    //in Datei wird n, x_n, abs(x_n-exakt_n) gespeichert.

    fprintf(fp,"n\tx_n\t\t\t\t\t\t\t\t|x_n-exakt_n|\n");
    for(int i=0;i<n;i++){
        fprintf(fp, "%d\t%25.19e\t%25.19e\n",i+1,result[i],fabs(result[i]-exact[i]));
    }

    fclose(fp);
    
    free(result);
    
    
    
}
