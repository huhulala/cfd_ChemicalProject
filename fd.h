#ifndef FD_H_
#define FD_H_

/**
 * d²u/dx²
 */
double  d2udx2( int i,  int j,  double **U, double dx );

/**
 * d²u/dy²
 */
double  d2udy2( int i,  int j,  double **U, double dy );

/**
 * du²/dx
 */
double  du2dx( int i,  int j,  double **U, double dx, double alpha );

/**
 * duv/dy
 */
double  duvdy( int i,  int j,  double **U, double **V, double dy, double alpha );

/**
 * d²v/dx²
 */
double  d2vdx2( int i,  int j,  double **V, double dx );

/**
 * d²v/dy²
 */
double  d2vdy2( int i,  int j,  double **V, double dy );

/**
 * dv²/dy
 */
double  dv2dy( int i,  int j,  double **V, double dy, double alpha );

/**
 *  duv/dx
 */
double  duvdx( int i,  int j,  double **U, double **V, double dx,double alpha );


#endif /* FD_H_ */
