/*******************************************
*Komplieren mit
*gcc h021.c -Wall -pedantic -o h21 -lm
*aufrufen mit
*./h21
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

/*auf 1 nomieren*/
double* normierung(double x0, double L, int N, double V, double y0,double lambda){
    double h=L/N;
    double* y =malloc((N+1)*sizeof(double));
    y=numerov(x0,L,N,V,y0,lambda);

    double integral=0.;
    
    for(int i=0; i<N+1;i++){
        integral+=(y[i]*y[i])*h;
    }
    
    double norm=1/sqrt(integral);
    
    for(int i=0;i<N+1;i++){
        y[i]*=norm;
    }
    
    return y;
}

/*Main Funktion*/
int main( int argc, char ** argv ){
    //Gegebene Parametern
    double L=8.;
    double Vmax=60.;
    double x0 = 0.0;
    double y0 =0.0;
    double yL=0.0;
    
    int N=1000;
    double eps=1e-5;
    //um y-werte zu speichern
    double* y =malloc((N+1)*sizeof(double));


    
    //Fuer lambda
    double lambwert=0.01;
    double lambstep=0.01;

    //um eigenwerte zu speichern
    int lambda_size = (int)(Vmax / lambstep) + 1;
    double *eigenwert = (double*)malloc(lambda_size * sizeof(double)); 
    
    int counter=0;
    
        /*Eigenwerte  werden in einer Txt datei gespeichert*/
    FILE*fp=fopen("h2_eigenwerte.txt","w");
    
    /*Lambda Werte werden bestimmt*/
    while(lambwert<= Vmax){
        y = numerov(x0, L,N,Vmax,y0,lambwert); //y wird berechnet

        if(fabs(y[N]-yL)<eps){ //wenn Anfangsbed erfuellt ist.
            eigenwert[counter]=lambwert; //is die lambda Eigenwert
            printf("%.2lf\n",lambwert);
            
            fprintf(fp,"%.2lf\n",eigenwert[counter]); //in Datei speichern
            counter++;
        }
        
        lambwert+=lambstep;//lambda erhoehen
    }
    
    fclose(fp); //datei schliessen
    printf("File Close\n");
    
    /*lambda-Grenze der Luecken werden bestimmt*/
    double lamb_grenze[6];
    lamb_grenze[0]=eigenwert[0]; //Anfangswert
    lamb_grenze[5]=eigenwert[counter-1];//Endwert
    int index=1;
    
    while(index<5){ //lamb_grenze[index]
        for(int j =1;j<counter-2;j++){ //eigenwert[j]
            if((eigenwert[j+1]-eigenwert[j])>0.02){
                lamb_grenze[index]=eigenwert[j];
                index++;
                
                lamb_grenze[index]=eigenwert[j+1];
                index++;
            }//ende if
        }//ende for schliefe
    }//ende while schliefe
    
    
    /*nomierte Wellenfunktion fuer die Grenze berechnet*/
    double* ynorm =malloc((N+1)*sizeof(double));
    double* x=malloc((N+1) *sizeof(double));
    
    //x Werte
    for(int i=0; i<N+1;i++){
        x[i]=i*L/N;
    }
    

    printf("Grenzwerte fuer Bankluecken sind: ");
    for(int i=0;i<6;i++){        
        printf("%.2lf\t",lamb_grenze[i]);//Datei erstellen

        char filename[32];
        sprintf(filename,"h2_Band_%d.txt",i+1);
        FILE*f=fopen(filename,"w");
        
        ynorm=normierung(x0,L,N,Vmax,y0,lamb_grenze[i]);//nomierte Werte rechnen.
        
        for(int j=0;j<N+1;j++){
            fprintf(f,"%.4lf\t\t%1.15e\n",x[j],ynorm[j]);//in Datei speichern.
        }
        
        fclose(f);//Datei schliessen
        
    }
    printf("\n");

    /*speicher frei geben*/
    free(eigenwert);
    free(x);
    free(ynorm);
    
    return 0;
}

        
   