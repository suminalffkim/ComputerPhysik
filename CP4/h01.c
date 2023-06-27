/***************************************
*kompilieren mit 
*gcc h01.c -Wall -pedantic -o h1 -lm
*aufrufen mit
*./h1
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//s(x)=0
double sf(double x){
    return 0;
}

//g(x)=lambda
double gf(double x, double lambda){
    return lambda;
}

/*Numerov_Verfahren. Hier wird u-Werte berechnet*/
double* numerov(double x0,double x1,double u0 ,double u1,double h,double lambdawert){
    //Schritte
    int N=(x1-x0)/h;

    // Allocate memory for the array
    double* x = (double *)malloc((N + 1) * sizeof(double));
    // Allocate memory and initialize with zeros
    double* u = (double *)calloc(N + 1, sizeof(double));

    // Generiere Werte fuer x
    for (int i = 0; i <N+1; i++) {
        x[i] = x0 + i * h;
    }
    
    double factor = h * h / 12;
    u[0] = u0;  //Anfangswert
    u[1] = u1; 
    
    for (int i=1;i<N;i++){ //Vorwaertsiteration
        u[i+1] = ((2 - 10 * factor * gf(x[i], lambdawert)) * u[i] - (1 + factor * gf(x[i-1], lambdawert)) * u[i-1] + factor * (sf(x[i+1]) + 10 * sf(x[i]) + sf(x[i-1]))) / (1 + factor * gf(x[i+1], lambdawert));
    }
    
    return u;
}


/*theoretische Wert fuer die Eigenwerte*/
double* lamb_theoretisch(double x0, double x1, double u0, double u1, int num){
    double* lamb_guess=(double*) malloc((num+1)*sizeof(double));
    
    for(int i=0;i<num;i++){
        lamb_guess[i]=i*i*M_PI*M_PI/((x1-x0)*(x1-x0));
        
    }
    return lamb_guess;
}


/*findet maximum*/
double find_maximum(double* u,int N){
    double max_value=u[0];
    
    for(int i=1;i<N;i++){
        if(u[i]>max_value){
            max_value=u[i];
        }
    }
    return max_value;
}


/*Main Funktion*/
int main( int argc, char ** argv ){
//Gegebene Parametern
    double x0 = 0.0;
    double x1 = 60.0;
    double u0 = 0.0;
    double u1 = 0.0;
    double h = 0.01;
    int num = 10; //Anzahl eigenwerte
    
    double eigenwerte[num]; //eigenwerte werden in einer Array gespeichert.
    double eps=5e-4;
    
    //Fuer lambda
    double lambwert=0.001;
    double lambstep=0.001;

    int counter=0;
    
    /*Lambda Werte werden bestimmt*/
    while(counter<num){
        double* u = numerov(x0,x1,u0,1e-6,h,lambwert); 
        int len=(x1-x0)/h;
        double uwert=u[len];
    
        if (fabs(uwert-u1)<eps){ //wenn die anfangsbed erfuellt ist, wird lambda gespeichert
            eigenwerte[counter]=lambwert;
            printf("%.3lf\n",lambwert);
            counter+=1;
        }//Ende if schleife
        
        lambwert=lambwert+lambstep; //erhoeht lambdawert
    }
    
    //theoretische Loeseung fuer Lambda
    printf("Theoretische Eigenwerte sind: ");
    double* lamb_theo=lamb_theoretisch(x0,x1,u0,u1,num);
    for(int i=0;i<num;i++){
        printf("%.4lf\t",lamb_theo[i]);
    }
    printf("\n\n");
    
    //mit u0=1,u1=1
    u0=1.;
    for(int i=0;i<num;i++){
        double* u=numerov(x0,x1,u0,1e-6,h,lamb_theo[i]); //berechnet u-Werte fuer verschidene lambda
        double max_value=find_maximum(u,(int)(x1-x0)/h); //max. Wert bestimmen
        printf("Eigenwert lambda: %.4lf    Max Naeherungswert:%.4lf\n",lamb_theo[i],max_value); //Ausgabe
    }   
    
    

}
