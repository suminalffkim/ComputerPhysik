/****************************************
*kompileren mit 
*gcc -Wall -pedantic madelung.c -o madelung -lm

*aufruf mit
*./madelung
*******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>  //um zufallsvariable zu generieren.
    
#define PI (3.141592653589793)
#define epsilon (8.8541878128*pow(10,-12)) //Vakuum permittivitaet in Einheit F/m
#define elek (1.602176634*pow(10,-19)) //elektron Ladung
                 
    
/*************************************************
*Hier wird es bestimmt ob es Natrium oder Chlorid handelt.
*Sei auf der Position (0,0,0) ein Chlorid.
*****************************************************/
int NaCl(long long int pos[3]){ 
    int possum=pos[0]+pos[1]+pos[2]; // Summe alle Indices der Positionen
    if (possum%2==0){ //Chlorid
        return -1;
    }
    else{ //Natrium
        return 1;
    }
}

                 
/************************************************************
*Abstand zwischen zwei Baustein
*************************************************************/
double getAbstand(long long int elem1[3], long long int elem2[3]){
    return(sqrt(pow(elem1[0]-elem2[0],2)+pow(elem1[1]-elem2[1],2)+pow(elem1[2]-elem2[2],2))); //Gitterabstand "a" wurde weggelassen.
}
           
           
/***********************************************
*Funktion getEinzelEnergie ermittelt einzelne Madelung-Energie V_ij
*************************************************/
double getEinzelEnergie(long long int elem1[3], long long int elem2[3]){
    double r=getAbstand(elem1,elem2);
    
    if (r==0){ //fuer flache Kristall. wenn r=0 ist wird Energie inf, deswegen wird r=0 weggelassen
    return 0;
    }
    
    else{
    return (1/(4*PI*epsilon)*NaCl(elem1)*elek*elek*NaCl(elem2)/r);    
    }
}


/*********************************************************************
*Funktion getGesEnergie berechnert gesamte Energie einer Baustein V_i
***********************************************************************/
double getGesEnergie(int Anzahl){
    long long int N= pow(Anzahl,3); //Anzahl der gesamte Atome
    long long int elem[N][3]; //Plaetze der jeweile Atomen
    
    //die Position von alle einzelne Bausteine erstellen.
    int row=0;
        for (long long int x=0;x<Anzahl;x++){
            //printf("x=%ld\n",x);
            for(long long int y=0;y<Anzahl;y++){
                //printf("y=%ld\n",y);
                for(long long int z=0;z<Anzahl;z++){
                   // printf("z=%ld\n",z);
                    elem[row][0]=x;
                    elem[row][1]=y;
                    elem[row][2]=z;
                    row++;
                    
                }
            }
        } 
    
    //Ausgabe einzelne Atomposition
//    for (int a=0;a<N;a++){
//        for (int b=0;b<3;b++){
//            printf("%d",elem[a][b]);
//        }
//        printf("\n");
//    }
    
    
    long long int elem1[3],elem2[3];
    srand(time(NULL));   // Initialization, should only be called once.
    long long int zufall = rand()%N+1; //generiert zufaellige Nummer zwischen 1 und N
    long long int ion=zufall; //Waehle ein Atom zufaellig  
    double ges_energie=0.;
    
    //printf("%lf gesamte Energie vorher  %d\n",ges_energie,Anzahl);
    
    for (long long int j=1;j<N;j++){//j_1=j_2=j_3=0 wird weggelassen. 
        if(j!=ion){
            for (long long int k=0;k<3;k++){
                elem1[k]=elem[ion][k];
                elem2[k]=elem[j][k];
             }
            ges_energie+=getEinzelEnergie(elem1,elem2); //Gesamt energie des i-te Baustein

        }
    }
    //printf("%25.16e\n",ges_energie);
    return (ges_energie);
}
                 
     
/************************************************************************
*Flache Kristall
*j_3=i_3
****************************************************************************/
double getGesEnergieFlach(int Anzahl){
    long long int N= pow(Anzahl,3); //Anzahl der gesamte Atome
    long long int elem[N][3]; //Plaetze der jeweile Atomen
    
    //die Position von alle einzelne Bausteine erstellen.
    long long int row=0;
        for (long long int x=0;x<Anzahl;x++){
            //printf("x=%d\n",x);
            for(long long int y=0;y<Anzahl;y++){
                //printf("y=%d\n",y);
                for(long long int z=0;z<Anzahl;z++){
                   // printf("z=%d\n",z);
                    elem[row][0]=x;
                    elem[row][1]=y;
                    elem[row][2]=z;
                    row++;
                    
                }
            }
        } 
    long long int elem1[3],elem2[3];
    srand(time(NULL));   // Initialization, should only be called once.
    long long int zufall = rand()%N+1; //generiert zufaellige Nummer zwischen 1 und N
    long long int ion=zufall; //Waehle ein Atom zufaellig  
    double ges_energie=0.;
    
    //printf("%lf gesamte Energie vorher  %d\n",ges_energie,Anzahl);
    
    for (long long int j=1;j<N;j++){//j_1=j_2=j_3=0 wird weggelassen. 
        if(j!=ion){
            for (long long int k=0;k<3;k++){
                elem1[k]=elem[ion][k];
            }
            for (long long int l=0;l<2;l++){
                elem2[l]=elem[j][l];
             }
            elem2[2]=elem1[2];
            ges_energie+=getEinzelEnergie(elem1,elem2); //Gesamt energie des i-te Baustein

        }
    }
    //printf("%25.16e\n",ges_energie);
    return (ges_energie);
}
           

/***********************************************************************
 * MAIN PROGRAMM
 ***********************************************************************/    
int main(int argc, char** argv) {
    int Anzahl=70; //Die Anzahl der Atome pro Linie  !Ab 7x Segmentation fault! 
    double V_i[Anzahl-1][2]; // Array fuer Madelung Energie
    double V_ifl[Anzahl-1][2]; //fuer flache kristall

    
    //Energie berechnen.
    long long int row=0;
    for(long long int counter=2;counter<=Anzahl;counter++){
        V_i[row][0]=pow(counter,3); //in erste column wird die gesamte Anzahl der Atome gespeichert.
        V_i[row][1]=getGesEnergie(counter); // in zweite Column die Madelungsenergie
        row++;
    }
    //printf("Anzahl  Energie\n");
    /*for(long long int i=0;i<row;i++){
       // printf("%d  %25.16e\n",i+2,V_i[i][1]);
    }
    */
    //Die Ergebnisse werden in einer Text Datei gespeichert. 
    FILE *f = fopen("madelung.txt", "w");
    for(long long int i=0;i<row;i++) {
        for(int j=0;j<2;j++) {
            fprintf(f,"%25.13e",V_i[i][j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
    
    
    //Flache Kristall
    row=0;
    for(long long int counter=2;counter<=Anzahl;counter++){
        V_ifl[row][0]=pow(counter,3);
        V_ifl[row][1]=getGesEnergieFlach(counter);
        row++;
    }
    
   //Die Ergebnisse werden in einer Text Datei gespeichert. 
    FILE *fp = fopen("madelungFlach.txt", "w");
    for(long long int i=0;i<row;i++) {
        for(int j=0;j<2;j++) {
            fprintf(fp,"%25.13e",V_ifl[i][j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp); 
    
    
        
return 0;
    
}