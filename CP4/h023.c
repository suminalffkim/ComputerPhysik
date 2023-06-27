/*******************************************
*Komplieren mit
*gcc h023.c -pedantic -Wall -o h23 -lm
*aufrufen mit
*./h23
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
    
    /*s(x)-Funktion =0*/
double sf(double x){
    return 0;
}

    /*g(x)-Funktion*/
double gf(double x,double lambda,double epsilon){
    double res=2*lambda-2*60*pow(cos(M_PI*x),16)-x*epsilon;
    return res;
}

 /*Numerov_Verfahren*/
double* numerov (double x0, double L, int N, double V, double y0, double lambdawert,double epsilon){
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
        y[i+1] = ((2 - 10 * factor * gf(x[i], lambdawert,epsilon)) * y[i] - (1 + factor * gf(x[i-1], lambdawert,epsilon)) * y[i-1] + factor * (sf(x[i+1]) + 10 * sf(x[i]) + sf(x[i-1]))) / (1 + factor * gf(x[i+1], lambdawert,epsilon));
    }
    
    /*speicher freigeben*/
    free(x);
    
    return y;

}  


/*auf 1 nomieren*/
double* normierung(double x0, double L, int N, double V, double y0,double lambda,double epsilon){
    double h=L/N;
    double* y =malloc((N+1)*sizeof(double));
    y=numerov(x0, L,N,V,y0,lambda,epsilon);

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


    //Fuer Feldstaerke epsilon
    double epsilon[3]={0.1,1.0,2.0};
    
    //fuer alle epsilon-werte
    for(int i=0;i<3;i++){

         //Fuer lambda
        double lambwert=0.01;
        double lambstep=0.01;
        
        //um eigenwerte zu speichern
        int lambda_size = (int)(Vmax / lambstep) + 1;
        double *eigenwert = (double*)malloc(lambda_size * sizeof(double)); 
        
        //um y-werte zu speichern
        double* y =malloc((N+1)*sizeof(double));


        
        /*Eigenwerte werden in einer Txt datei gespeichert*/        
        char filename[32];
        sprintf(filename,"h23_Bandluecke_%.1lf.txt",epsilon[i]);
        FILE*fp=fopen(filename,"w");
        printf("File Open\n");
                
        int counter=0; //fuer while schliefe  
        
        /*Lambda Werte werden bestimmt*/
        while(lambwert<= Vmax){
            y = numerov(x0,L,N,Vmax,y0,lambwert,epsilon[i]); //y wird berechnet

            if(fabs(y[N]-yL)<eps){ //wenn Anfangsbed erfuellt ist.
                eigenwert[counter]=lambwert;
                //printf("%.2lf\n",lambwert);
                fprintf(fp,"%.2lf\n",lambwert);
            
                //fprintf(fp,"%.2lf\n",eigenwert[counter]);
                counter++;
            }        
            lambwert+=lambstep;
        }//Ende erste While Schliefe
        
        free(y);//Speicher freigeben
        
        
        /*Grenzewerte fuer Band luecke werden bestimmt*/
        double lamb_grenze[4]; //fuer erste beide Bandluecke
        int index=0;
        printf("Grenzwerte fuer die Bandluecke sind: ");
        
        while(index<4){
            for(int j =1;j<counter-2;j++){
                if((eigenwert[j+1]-eigenwert[j])>0.02){
                    lamb_grenze[index]=eigenwert[j];
                    //fprintf(fp,"%.2f\n",lamb_grenze[index]);
                    index++;
                    
                    lamb_grenze[index]=eigenwert[j+1];
                    //fprintf(fp,"%.2f\n",lamb_grenze[index]);
                    index++;
                }
            }
        }//ende zweite while schleife
        
        //Grenzwerte ausgeben.
        for(int i=0;i<4;i++){
            printf("%.2lf\t",lamb_grenze[i]);
        }
        printf("\n");
        
        free(eigenwert); //speicher freigeben
            
            
        fclose(fp); //Datei Schliessen
        printf("File Close\n");
        
    }//Ende For Schliefe
    
    
    return 0;
}

        
   