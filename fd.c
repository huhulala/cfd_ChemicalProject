#include "fd.h"
#include <math.h>

double d2udx2(int i, int j, double **U, double dx) {
	return (U[i - 1][j] - 2 * U[i][j] + U[i + 1][j]) / (dx * dx);
}

double d2udy2(int i, int j, double **U, double dy) {
	return (U[i][j - 1] - 2 * U[i][j] + U[i][j + 1]) / (dy * dy);
}

double du2dx(int i, int j, double **U, double dx, double alpha) {
	return (1 / dx * (((U[i][j] + U[i + 1][j]) / 2) * ((U[i][j] + U[i + 1][j])
			/ 2) - ((U[i - 1][j] + U[i][j]) / 2)
			* ((U[i - 1][j] + U[i][j]) / 2)) + alpha / dx * (fabs(
			U[i][j] + U[i + 1][j]) / 2 * (U[i][j] - U[i + 1][j]) / 2 - fabs(
			U[i - 1][j] + U[i][j]) / 2 * (U[i - 1][j] - U[i][j]) / 2));
}

double duvdy(int i, int j, double **U, double **V, double dy, double alpha) {
	return (1 / dy * ((V[i][j] + V[i + 1][j]) * (U[i][j] + U[i][j + 1]) / 4
			- (V[i][j - 1] + V[i + 1][j - 1]) * (U[i][j - 1] + U[i][j]) / 4)
			+ alpha / dy * (fabs(V[i][j] + V[i + 1][j]) / 2 * (U[i][j] - U[i][j
					+ 1]) / 2 - fabs(V[i][j - 1] + V[i + 1][j - 1]) / 2
					* (U[i][j - 1] - U[i][j]) / 2));
}

double d2vdx2(int i, int j, double **V, double dx) {
	return (V[i - 1][j] - 2 * V[i][j] + V[i + 1][j]) / (dx * dx);
}

double d2vdy2(int i, int j, double **V, double dy) {
	return (V[i][j - 1] - 2 * V[i][j] + V[i][j + 1]) / (dy * dy);
}

double  duvdx( int i,  int j,  double **U, double **V, double dx, double alpha ){
	return (1/dx*((U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j])/4 -
			(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j])/4) +
			alpha/dx * (fabs(U[i][j] + U[i][j+1])*(V[i][j] - V[i+1][j])/4 -
					fabs(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] - V[i][j])/4));
}


double  dv2dy( int i,  int j,  double **V, double dy, double alpha ){
	return (1/dy*(((V[i][j] + V[i][j+1])/2)*((V[i][j] + V[i][j+1])/2) -
			((V[i][j-1] + V[i][j])/2)*((V[i][j-1] + V[i][j])/2)) +
			alpha/dy * (fabs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1])/4 -
					fabs(V[i][j-1] + V[i][j])*(V[i][j-1] - V[i][j])/4));
}
