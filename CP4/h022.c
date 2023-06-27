/*******************************************
*Komplieren mit
*gcc h022.c -pedantic -Wall -o h22 -lm
*aufrufen mit
*./h22
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
    
    /*s(x)_Funktion=0*/
double sf(double x){
    return 0;
}
    /*g(x)-Funktion*/
double gf(double x,double lambda){
    double res=2*lambda-2*60*pow(cos(M_PI*x),16);
    return res;
}
/*Numerov-Verfahren*/
double* numerov (double x0, double L, int N, double V, double y0, double lambdawert){
    double h=L/N; //Schritte
    double* x=malloc((N+1) *sizeof(double));
    double* y =malloc((N+1)*sizeof(double));
    
    //x Werte
    for(int i=0; i<N+1;i++){
        x[i]=x0+i*L/N;
    }
    
    y[0]=y0;//anfangswert 
    y[1]=1e-8; //Schaetzwert
    
    double factor = h * h / 12;

    
    for(int i=1;i<N;i++){ //Vorwaertsiteration
        y[i+1] = ((2 - 10 * factor * gf(x[i], lambdawert)) * y[i] - (1 + factor * gf(x[i-1], lambdawert)) * y[i-1] + factor * (sf(x[i+1]) + 10 * sf(x[i]) + sf(x[i-1]))) / (1 + factor * gf(x[i+1], lambdawert));
    }
    
    /*speicher freigeben*/
    free(x);
    
    return y;

}   


/*Main Funktion*/
int main( int argc, char ** argv ){
    //Gegebene Parametern
    double L[]={16.0,32.0,64.0};
    double Vmax=60.;
    double x0 = 0.0;
    double y0 =0.0;
    double yL=0.0;
    
    int N=1000;
    double eps=1e-5;

    /*Fuer jede L-Werte*/
    for(int i=0;i<3;i++){    
        double* y =malloc((N+1)*sizeof(double));

        //Fuer lambda
        double lambwert=0.01;
        double lambstep=0.01;

        //um eigenwerte zum speichern
        int lambda_size = (int)(Vmax / lambstep) + 1;
        double *eigenwert = (double*)malloc(lambda_size * sizeof(double)); 
    
        char filename[32];
        sprintf(filename,"h22_%2.0f_eigenwerte.txt",L[i]);
        
        /*Eigenwerte werden in einer Txt datei gespeichert*/
        FILE*fp=fopen(filename,"w");
        
        int counter=0;
        
        /*Lambda Werte werden bestimmt*/
        while(lambwert<= Vmax){
            y = numerov(x0, L[i],N,Vmax,y0,lambwert); //y wird berechnet

            if(fabs(y[N]-yL)<eps){ //wenn Anfangsbed erfuellt ist.
                eigenwert[counter]=lambwert; //is lambda Eigenwert
                printf("%.2lf\n",lambwert);
            
            
                fprintf(fp,"%.2lf\n",eigenwert[counter]);
                counter++;
            }//ende if schliefe
        
            lambwert+=lambstep;//lambda erhoehen
        }//ende while
        
        fclose(fp); //datei schliessen
        printf("File Close\n");
        
        free(eigenwert);//speicher freigeben
        
    }//Ende for schliefe
    
    return 0;
}

        
   