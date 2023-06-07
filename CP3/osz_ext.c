/*************************************
*kompiliere mit
*gcc -Wall -pedantic osz_ext.c -o ext -lm
*aufrufen mit
*./ext
**********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
  
/*****************************************
*DGL-Funktion
********************************/
double f(double t, double x, double xdot, double omega_0, double gamma, double alpha, double omega) {
    return (-omega_0 * omega_0 * x - 2 * gamma * xdot + alpha * cos(omega * t));
}


/*****************************************
*Runge-Kutta-Verfahren 4.Ord.
*******************************************/
double* RK4(double v_0, double omega_0, double gamma,double omega, double t_0, double x_0, double alpha, double h,unsigned int n){
    
    double x=x_0; //Anfangswert fuer x(t)
    double xdot=v_0; //Anfangswert fuer xdot(t)
    double j_1,j_2,j_3,j_4,k_1,k_2,k_3,k_4,t;
    double* x_res=malloc(n*sizeof(double));//x_n wird gespeichert.
    
    x_res[0]=x; //0-te Wert fuer x ist x_0

    for (int k = 1; k < n; k++) {
            t=t_0+k*h; //zeit nach Schritte berechnet
        
            j_1=h * f(t,x,xdot,omega_0,gamma,alpha,omega);
            k_1=h * xdot;
            j_2=h * f(t+h/2, x+k_1*0.5, xdot+j_1*0.5, omega_0, gamma,alpha,omega);
            k_2=h * (xdot+j_1*0.5);
            j_3=h * f(t+h/2, x+k_2*0.5, xdot+j_2*0.5, omega_0, gamma,alpha,omega);
            k_3=h * (xdot+j_2*0.5);
            j_4=h * f(t+h, x+k_3, xdot+j_3, omega_0,gamma,alpha,omega);
            k_4=h * (xdot+j_3);
            
            xdot+=(1./6.) * (j_1+2*j_2+2*j_3+j_4); //neue xdot= xdot_(k+1)
            x+=(1./6.) * (k_1+2*k_2+2*k_3+k_4); //neue x= x_(k+1)
        
            x_res[k]=x; //neue x wird in array gespeichert.
            //printf("%25.19e\n",j_4);
    }
    
    return x_res; //gib das Ergebnis zurueck
}



/*******************************************
* Die Funktion berechnet die exakten Werte
********************************************/
double* schwingung_gl(double v_0, double omega_0, double gamma, double alpha, double omega, double x_0, double t_0, double h, unsigned int n) {
    double t;
    double* schwing_res = malloc(n * sizeof(double));

             //      double y=sqrt(pow((omega_0*omega_0)-(omega*omega),2)+pow((2*gamma*omega),2));
        double phi=atan((2*gamma*omega)/(omega_0*omega_0-omega*omega));
                double omega_hat=sqrt(omega_0*omega_0-gamma*gamma);                     
    for (int i = 0; i < n; i++) {
        
        t = t_0 + i * h;


        schwing_res[i]=exp(-gamma*t)*(x_0*cos(omega_0*t)+(v_0+gamma*x_0)/omega_0*sin(omega_0*t))+alpha/(omega_0*omega_0-omega*omega)*cos(omega*t+phi);    
    }

                    
    return schwing_res;
    }

/********************************************************************
*Die Funktion findet Differenze einzelne Periode  =>Segmentation Fault
**************************************************************************/
double* diff_pro_Periode(double* res, double* ext ,unsigned int n) {
    int counter = 1;
    int index = 0;
    int* max_ind=malloc(20*sizeof(int));
    

       //die Maximumstelle werden gespeichert.
    for (int i = 0; i < n - 1 && index < 20; i++) {
        switch (counter % 2) {
            case 1:
                if (res[i + 1] < res[i]) { // Maximum erreicht
                    max_ind[index] = i;
                    counter++;
                    index++;
                }
                break;
            case 0:
                if (res[i + 1] > res[i]) { // Minimum erreicht
                    counter++;
                }
                break;
        }
    
    }
    
    double* diff=malloc(20*sizeof(double));
    int k=0;
    while(k<20){
        double rk4=0;
        double exact=0;
        
        
        for(int j=max_ind[k];j<max_ind[k+1];j++){ //fuer jede Schwingvorgang wird die Differenz berechnet.
            rk4+=res[j];
            exact+=ext[j];           
        }
        diff[k]=fabs(rk4-exact);
        k++;
    }
    return diff;
}
/*******************************
*main
*********************************/
int main(int argc, char **argv){
    double v_0=2.75;
    double omega_0=1.;
    double x_0=1.5;
    double gamma=0.1;
    double alpha=0.5;
    double omega=1.5;
    double t_0=0.;
    int n[]={10000,20000,20000};
    double h[]={0.01,0.005,0.001}; //drei h Werte werden eingesetzt und Ergebnisse vergliechen.
    
    for (int ih=0;ih<3;ih++){
        
        printf("File open\n");//zur kontrolle
        
        char filename[32];
        sprintf(filename, "Runge_Kutta_ext_%.3lf.txt",h[ih]); //erstelle fuer jede h Wert ein Dokument
    
        FILE *fp = fopen(filename, "w"); //File offenen
        
        
    
        double* result=RK4(v_0,omega_0,gamma,omega,t_0,x_0,alpha,h[ih],n[ih]); //Ergebnisse die mit Runge Kutta Verfahren berechnet wurden
    
        double* exact=schwingung_gl(v_0, omega_0,gamma,alpha,omega,x_0,t_0,h[ih],n[ih]); //Ergebnisse, die mit exakte Loesung berechnet wurden.
    
        
        //in Datei wird n, t, x_n, exact_n, abs(x_n-exakt_n) Werte gespeichert.

        fprintf(fp,"n\tt\t\t\t\t\t\t\t\tx_n\t\t\t\t\t\t\texact\t\t\t\t\t\t\t|x_n-exakt_n|\n");
        
        for(int i=0;i<n[ih];i++){
            double t=t_0+h[ih]*i;
            
            fprintf(fp, "%d\t%25.19e\t%25.19e\t%25.19e\t%25.19e\n",i,t,result[i],exact[i],fabs(result[i]-exact[i]));
        }

        fclose(fp);
        printf("File close\n");
        //double* diff=diff_pro_Periode(result,exact,n);
        free(result);
        free(exact);
        //free(diff);
        }
    
    
}
