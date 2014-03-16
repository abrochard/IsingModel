#include <stdlib.h>
#include <stdio.h>
#include "lattice.h"
/*
lattice.c
----------
contains methods for dealing with lattices
*/

int** createIntLattice(int m, int n) //creates a lattice of size nxm for integers
{
    int** lattice;
    int i=0;
    lattice = (int**)malloc(m * sizeof(int *));
	if(lattice == NULL)
    {
		fprintf(stderr, "out of memory\n");
	}
	for(i = 0; i < m; i++){
		lattice[i] = malloc(n * sizeof(int));
		if(lattice[i] == NULL)
		{
			fprintf(stderr, "out of memory\n");
		}
    }
    lattice=zeroIntLattice(lattice,m,n);
    return lattice;
}

int** zeroIntLattice(int** lattice, int m, int n) //sets all values of an int lattice to zero
{
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            lattice[i][j]=0;
        }
    }
    return lattice;
}

float** createFloattLattice(int m, int n) //creates a lattice of size nxm for integers
{
    float** lattice;
    int i=0;
    lattice = (float**)malloc(m * sizeof(float *));
	if(lattice == NULL)
    {
		fprintf(stderr, "out of memory\n");
	}
	for(i = 0; i < m; i++){
		lattice[i] = malloc(n * sizeof(float));
		if(lattice[i] == NULL)
		{
			fprintf(stderr, "out of memory\n");
		}
    }
    lattice=zeroFloatLattice(lattice,m,n);
    return lattice;
}

float** zeroFloatLattice(float** lattice, int m, int n) //sets all flaot values to zero
{
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            lattice[i][j]=0;
        }
    }
    return lattice;
}

int** posBoundary(int** lattice, int m, int n) //imposes + boundary condition
{
    int i;
    for(i=0;i<m;i++){
        lattice[i][0]=1;
        lattice[i][n-1]=1;
    }
    for(i=0;i<n;i++){
        lattice[0][i]=1;
        lattice[m-1][i]=1;
    }
    return lattice;
}

int** setRandomSpin(int** lattice, int m, int n) //sets a random spin on every vertex
{
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            lattice[i][j]=rand()%2;
        }
    }
    return lattice;
}

int getTriNeighborsSpin(int** lat, int m, int n) //returns the 6 neighbors of a point
{
    int neighbors[6]={0,0,0,0,0,0};
    neighbors[0]=lat[m][n-1];
    neighbors[1]=lat[m+1][n-1];
    neighbors[2]=lat[m+1][n];
//    neighbors[2]=lat[m+1][n+1];
    neighbors[3]=lat[m][n+1];
//    neighbors[4]=lat[m-1][n+1];
    neighbors[4]=lat[m-1][n];
    neighbors[5]=lat[m-1][n-1];
    int i=0, sum=0;
    for(i=0; i<6; i++){
        sum +=neighbors[i];
    }
    return sum;
}

int getSquareNeighborsSpin(int** lat, int m, int n) //returns the 4 neighbors of a point
{
    int neighbors[4]={0,0,0,0};
    neighbors[0]=lat[m][n-1];
    neighbors[1]=lat[m+1][n];
    neighbors[2]=lat[m][n+1];
    neighbors[3]=lat[m-1][n];
    int i=0, sum=0;
    for(i=0; i<4; i++){
        sum +=neighbors[i];
    }
    return sum;
}

int getHexNeighborsSpin(int **lat, int m, int n) //returns the 3 neighbors of a point
{
    int neighbors[3]={0};
    if(m%2==0)
    {
        if(n%2==0)
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m+1][n];
            neighbors[2]=lat[m][n+1];
        }
        else
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m][n-1];
            neighbors[2]=lat[m+1][n];
        }
    }
    else
    {
        if(n%2==0)
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m][n-1];
            neighbors[2]=lat[m+1][n];
        }
        else
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m+1][n];
            neighbors[2]=lat[m][n+1];
        }
    }
    int i=0, sum=0;
    for(i=0;i<3;i++){
        sum+=neighbors[i];
    }
    return sum;
}

void showLattice(float **lat, int m, int n) //to display the lattice
{
    int i=0, j=0;
    for(i=0;i<m;i++) {
        for(j=0;j<n;j++){
            printf("%f ",lat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void showLatticeI(int **lat, int m, int n) //to display the lattice
{
    int i=0, j=0;
    for(i=0;i<m;i++) {
        for(j=0;j<n;j++){
            printf("%d ",lat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void showTriLattice(int** lat, int m, int n) //to display a triangular lattice
{
    int i=0,j=0;
    for(i=0;i<m;i++){
        if(i%2==0)
            printf(" ");
        for(j=0;j<n;j++){
            printf("%d ",lat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printTriNeighbors(int** lat, int m, int n) //check function for triangular lattice
{
    int neighbors[6]={0,0,0,0,0,0};
    neighbors[0]=lat[m][n-1];
    neighbors[1]=lat[m+1][n];
    neighbors[2]=lat[m+1][n+1];
    neighbors[3]=lat[m][n+1];
    neighbors[4]=lat[m-1][n+1];
    neighbors[5]=lat[m-1][n];
    int i=0;
    for(i=0; i<6; i++){
        printf("%d ", neighbors[i]);
    }
}

void printSquareNeighbors(int** lat, int m, int n) //check function for square lattice
{
    int neighbors[4]={0,0,0,0};
    neighbors[0]=lat[m][n-1];
    neighbors[1]=lat[m+1][n];
    neighbors[2]=lat[m][n+1];
    neighbors[3]=lat[m-1][n];
    int i=0;
    for(i=0; i<4; i++){
        printf("%d ", neighbors[i]);
    }
}

void printHexNeighbors(int **lat, int m, int n) //check function for hex lattice
{
    int neighbors[3]={0,0,0};
    if(m%2==0)
    {
        if(n%2==0)
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m+1][n];
            neighbors[2]=lat[m][n+1];
        }
        else
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m][n-1];
            neighbors[2]=lat[m+1][n];
        }
    }
    else
    {
        if(n%2==0)
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m][n-1];
            neighbors[2]=lat[m+1][n];
        }
        else
        {
            neighbors[0]=lat[m-1][n];
            neighbors[1]=lat[m+1][n];
            neighbors[2]=lat[m][n+1];
        }
    }
    int i=0;
    for(i=0; i<3; i++){
        printf("%d ", neighbors[i]);
    }
}

void eraseLattice(float** lat, int m, int n) //erase the lattice from the memory
{
    int i;
    for(i=0;i<m;i++){
        free(lat[i]);
    }
    free(lat);
}
