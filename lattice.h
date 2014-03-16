/*
lattice.h
*/

int** createIntLattice(int m, int n);

int** zeroIntLattice(int** lattice, int m, int n);

float** createFloattLattice(int m, int n);

float** zeroFloatLattice(float** lattice, int m, int n);

int** posBoundary(int** lattice, int m, int n);

int** setRandomSpin(int** lattice, int m, int n);

int getTriNeighborsSpin(int** lat, int m, int n);

int getSquareNeighborsSpin(int** lat, int m, int n);

int getHexNeighborsSpin(int **lat, int m, int n);

void showLattice(float **lat, int m, int n);

void showLatticeI(int **lat, int m, int n);

void showTriLattice(int** lat, int m, int n);

void printTriNeighbors(int** lat, int m, int n);

void printSquareNeighbors(int** lat, int m, int n);

void printHexNeighbors(int **lat, int m, int n);

void eraseLattice(float** lat, int m, int n);
