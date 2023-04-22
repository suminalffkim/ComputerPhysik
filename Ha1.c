/****************************************
*kompileren mit 
*gcc -Wall -pedantic XXXX.c -o xxx -lm
*aufruf mit
*./ xxx
*******************************************/
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
    
#define PI (3.141592653589793)
#define epsilon (8.8541878128*pow(10,-12) //Vakuum permittivitaet in Einheit F/m

                 
/*************************************************
*Hier wird es bestimmt ob es Natrium oder Chlorid handelt.
*Sei auf der Position (0,0,0) ein Chlorid ist. Mehr dazu auf der PDF Dati zu sehen.
*****************************************************/
int NaCl(int pos[3]){ 
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
double Abstand(int elem1[3], int elem2[3]){
    return(sqrt(pow(elem1[0]-elem2[0],2)+pow(elem1[1]-elem2[1],2)+pow(elem1[2]-elem2[2],2))); //Gitterabstand a wurde erstens weggelassen.
}
           
           
/***********************************************
*Funktion Energie ermittelt einzelne Madelung-Energie
*************************************************/
double Energie(int elem1[3], int elem2[3]){
    double r=Abstand(elem1,elem2);
    return (1/(4*PI*epsilon)*NaCl(elem1)*NaCl(elem2)/r));    
}
                 
           
           
/***********************************************************************
 * MAIN PROGRAMM
 ***********************************************************************/    
int main(int argc, char **argv) {
    int Anzahl=2; //Die Anzahl der Atome pro Linie
    int N= pow(Anzahl,3); //Anzahl der gesamte Atome
   // double r[N]; //Abstand zwischen Atomen
    int elem[N][3]; //Plaetze der jeweile Atomen
   // double V[N]; // Madelung Energie
    
    
   //die Position von alle einzelne Bausteine erstellen.
    int row=0;
        for (int x=0;x<Anzahl;x++){
            //printf("x=%d\n",x);
            for(int y=0;y<Anzahl;y++){
                //printf("y=%d\n",y);
                for(int z=0;z<Anzahl;z++){
                   // printf("z=%d\n",z);
                    elem[row][0]=x;
                    elem[row][1]=y;
                    elem[row][2]=z;
                    row++;
                    
                }
            }
        }  
    
    //Ausgabe einzelne Atomposition
    for (int a=0;a<N;a++){
        for (int b=0;b<3;b++){
            printf("%d",elem[a][b]);
        }
        printf("\n");
    }
    
    row=0;
    int elem1[3],elem2[3];
    double energie;

    for(int i=0;i<N;i++){ //i-te Baustein
        for (int j=1;j<N;j++){//j-te Baustein
            for (int k=0;k<3;k++){
                elem1[k]=elem[i][k];
                elem2[k]=elem[j][k];
            }
            energie=Energie(elem1,elem2); //여기서 부터 다시하기. 맞는 답이 나오는지 확인할것.
            printf("%lf\n",energie);
        }

    }
        
return 0;
    
}