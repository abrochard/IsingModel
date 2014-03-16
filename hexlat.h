/*
HexLat.h
--------
*/

float** createProbaTable(float** probas, int n, float t);

int** sweepMetropolis(int** lat, float** probUp, double nTimes, int m, int n, int x, int y, double stableTime, int** flipLattice);

