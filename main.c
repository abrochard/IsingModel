//defines the situation
//#define TRIANGULAR
#define SQUARE
//#define HEX
#define M 200
#define N 200
#define NTIMES 100000000
#define X 100
#define Y 100
#define STABLETIME 10000000

//debug mod
//#define SHOW
//#define FLIP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "HexLat.h"
#include "lattice.h"

/*
Ising Model on Triangular/Hex Lattice
Implementation of Metropolis algorithm
*/

int main()
{
    //time recording
    clock_t begin, end;
    double time_spent;
    begin = clock();
    //initialisation de rand
    srand(time(NULL));
    //create lattice with random spin and + b.c.
#ifdef TRIANGULAR
    printf("Triangular lattice, spin on vertices, random spin, +b.c.\n\n");
#endif
#ifdef SQUARE
    printf("Square lattice, spin on vertices, random spin, +b.c.\n\n");
#endif // SQUARE
#ifdef HEX
    printf("Hex lattice, spin on vertices, random spin, +b.c.\n\n");
#endif // HEX
    printf("Size:%dx%d\nIterations:%d\nStable time:%d\n",M,N,NTIMES,STABLETIME);
    printf("Position: (%d,%d)\n\n",X,Y);
    int** lattice=createIntLattice(M,N);
    lattice=setRandomSpin(lattice, M, N);
    lattice=posBoundary(lattice,M,N);
    int ** flipLattice=createIntLattice(M,N);
#ifdef SHOW
    showTriLattice(lattice,M,N);
    printHexNeighbors(lattice,2,2);
    printf("\n%d\n",getHexNeighborsSpin(lattice,2,2));
#endif
#ifdef SQUARE
    float t=2/log(sqrt(2)+1); //for square
#endif
#ifdef TRIANGULAR
    float t=4/log(3); //critical for triangular
#endif
#ifdef HEX
    float t=2/log(sqrt(3)+2); //for hex
#endif // HEX
    //create the probability table
#ifdef SQUARE
    float** probas=createFloattLattice(2,5);
    probas=createProbaTable(probas,4,t);
    printf("Proba table:\n");
    showLattice(probas,2,5);
#endif
#ifdef TRIANGULAR
    float** probas=createFloattLattice(2,7);
    probas=createProbaTable(probas,6,t);
    printf("Proba table:\n");
    showLattice(probas,2,7);
#endif // TRIANGULAR
#ifdef HEX
float** probas=createFloattLattice(2,4);
    probas=createProbaTable(probas,3,t);
    printf("Proba table:\n");
    showLattice(probas,2,4);
#endif // HEX
    printf("For t=%f\n",t);
    //launch the sweeping
    lattice=sweepMetropolis(lattice,probas,NTIMES,M,N,X,Y,STABLETIME,flipLattice);
    showLatticeI(lattice,M,N);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Computation time = %f s\n\n",time_spent);
    return 0;
}

float** createProbaTable(float** probas, int n, float t) //sets the probability for n neighbors and temperature t
/*
Nombres de voisins à +1 (de 0 à n)           | 0            | 1            | 2           | ... | n/2 | ... | n-1            | n             |
---------------------------------------------------------------------------------------------------------------------------------------------
Prob d'aller à +1 sachant que l'on est -1    | exp(1/T)^S1  | exp(1/T)^S1  | exp(1/T)^S1 | ... | 0.5 | ... |  1             | 1             |
---------------------------------------------------------------------------------------------------------------------------------------------
Prob de rester  +1                           | 0            | 0            | 0           | ... | 0.5 | ... | 1-exp(1/T)^S2  | 1-exp(1/T)^S2 |

S1=2*(sum of negative spins) (S1<0)
S2=-2*sum of positive spins (S2<0)
*/
{
    int nUp=0;
    for(nUp=0;nUp<n+1;nUp++){
        if ((float)nUp<(float)(n)/2)
        {
//            printf("n=%d nUp=%d t=%f %f %f\n",n, nUp, t,2*(n-nUp*2)/t,exp(2*(n-nUp*2)/t));
            probas[0][nUp]=exp(-2*(n-nUp*2)/t);
            probas[1][nUp]=0;
        }
        if ((float)nUp==(float)(n)/2)
        {
            probas[0][nUp]=0.5;
            probas[1][nUp]=0.5;
        }
        if ((float)nUp>(float)(n)/2)
        {
            probas[0][nUp]=1;
            probas[1][nUp]=1-exp(-2*(2*nUp-n)/t);
        }
    }
    return probas;
}

int** sweepMetropolis(int** lat, float** probUp, double nTimes, int m, int n, int x, int y, double stableTime, int** flipLattice)
{
    int time=0, i=0, j=0;
    int counter=0;
    int temp=0;
//    int edge=0;
    float prob=0;
    for(time=0;time<nTimes;time++){
            i=rand()%(m-2)+1;
            j=rand()%(n-2)+1;
            prob=rand()/(double)RAND_MAX;
//            printf("prob=%f spin=%d\n",prob,getSquareNeighborsSpin(lat,i,j));
#ifdef SQUARE
            if (prob<probUp[(lat[i][j])][getSquareNeighborsSpin(lat,i,j)])
            {
                temp=lat[i][j];
                lat[i][j]=1;
                if(temp != lat[i][j] && time>stableTime){
                    flipLattice[i][j]++;
                }
            }
            else
            {
                temp=lat[i][j];
                lat[i][j]=0;
                if(temp != lat[i][j] && time>stableTime){
                    flipLattice[i][j]++;
                }
            }
#endif
#ifdef TRIANGULAR
            if (prob<probUp[(lat[i][j])][getTriNeighborsSpin(lat,i,j)])
            {
                temp=lat[i][j];
                lat[i][j]=1;
                if(temp != lat[i][j] && time>stableTime){
                    flipLattice[i][j]++;
                }
            }
            else
            {
                temp=lat[i][j];
                lat[i][j]=0;
                if(temp != lat[i][j] && time>stableTime){
                    flipLattice[i][j]++;
                }
            }
#endif // TRIANGULAR
#ifdef HEX
            if (prob<probUp[(lat[i][j])][getHexNeighborsSpin(lat,i,j)])
            {
                temp=lat[i][j];
                lat[i][j]=1;
                if(temp != lat[i][j] && time>stableTime){
                    flipLattice[i][j]++;
                }
            }
            else
            {
                temp=lat[i][j];
                lat[i][j]=0;
                if(temp != lat[i][j] && time>stableTime){
                    flipLattice[i][j]++;
                }
            }
#endif // HEX
        if (lat[x][y]==lat[x][y+1] && time>stableTime)
        {
            counter++;
        }
//        if (lat[x][y]==1 && time>stableTime)
//        {
//            counter++;
//        }
//        if (lat[7][7]==1 && time>stableTime)
//        {
//            edge++;
//        }
#ifdef SHOW
        showTriLattice(lat,m,n);
#endif
    }
    printf("\nRatio = %f\n",(float)counter/(nTimes-stableTime));
//    printf("Center: %d Border: %d\nRatio=%f\n", counter, edge, (float)counter/edge);
#ifdef FLIP
#ifdef SQUARE
    showLatticeI(flipLattice,M,N);
#endif
#ifdef TRIANGULAR
    showLatticeI(flipLattice,M,N);
#endif // TRIANGULAR
#ifdef HEX
    showLatticeI(flipLattice,M,N);
#endif // HEX
#endif //FLIP
    return lat;
}

