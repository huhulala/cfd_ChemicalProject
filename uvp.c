#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include "fd.h"
#include <math.h>
#include "NSDefinitions.h"

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G, int **Flag) {
	/* see formulas 9 and 10 in combination with formulas 4 and 5 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			/********** calculate F **********/
			/* calculate f and g only between two fluid cells */
			if(((Flag[i][j]&B_C)==B_C)&& i<imax ){
				F[i][j] = U[i][j] + dt * (
				/* 1/Re * (d²u/dx² + d²u/dy²) */
				1 / Re * (d2udx2(i, j, U, dx) + d2udy2(i, j, U, dy))
				/* - du²/dx */
				- du2dx(i, j, U, dx, alpha)
				/* - duv/dy */
				- duvdy(i, j, U, V, dy, alpha) + GX);
			}

			/********** calculate G **********/
			if(((Flag[i][j]&B_C)==B_C)&& i<imax ){
				G[i][j] = V[i][j] + dt * (
				/* 1/Re * (d²v/dx² + d²v/dy²) */
				1 / Re * (d2vdx2(i, j, V, dx) + d2vdy2(i, j, V, dy))
				/* - duv/dx */
				- duvdx(i, j, U, V, dx, alpha)
				/* - dv²/dy */
				- dv2dy(i, j, V, dy, alpha) + GY);
			}

			if((Flag[i][j]&31)==B_N){
				G[i][j]=V[i][j];
			}
			else if((Flag[i][j]&31)==B_S){
				G[i][j-1]=V[i][j-1];
			}
			else if((Flag[i][j]&31)==B_W){
				F[i-1][j]=U[i-1][j];
			}
			else if((Flag[i][j]&31)==B_O){
				F[i][j]=U[i][j];
			}
			else if((Flag[i][j]&31)==B_NO){
				F[i][j]=U[i][j];
				G[i][j]=V[i][j];
			}
			else if((Flag[i][j]&31)==B_NW){
				F[i-1][j]=U[i-1][j];
				G[i][j]=V[i][j];
			}
			else if((Flag[i][j]&31)==B_SO){
				F[i][j]=U[i][j];
				G[i][j-1]=V[i][j-1];
			}
			else if((Flag[i][j]&31)==B_SW){
				F[i-1][j]=U[i-1][j];
				G[i][j-1]=V[i][j-1];
			}

		}

	/* calculate boundary values -  see formula 17 */
	for (j = 1; j <= jmax; j++) {
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for (i = 1; i <= imax; i++) {
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V) {
	/* See formula 13 */
	double umax = fabs(U[1][1]);
	double vmax = fabs(V[1][1]);
	double dtcon, dxcon, dycon;
	double min;
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			if (fabs(U[i][j]) > umax)
				umax = fabs(U[i][j]);
			if (fabs(V[i][j]) > vmax)
				vmax = fabs(V[i][j]);
		}

	/* conditions */
	dtcon = Re / (2 * (1 / (dx * dx) + 1 / (dy * dy)));
	dxcon = dx / fabs(umax);
	dycon = dy / fabs(vmax);

	/* determine smalles condition */
	min = dtcon;
	if (min > dxcon)
		min = dxcon;
	if (min > dycon)
		min = dycon;

	/* calculate dt */
	*dt = tau * min;
}

void calculate_dt1(
		double Re,
		double tau,
		double *dt,
		double dx,
		double dy,
		int imax,
		int jmax,
		double **U,
		double **V,
		int **Flag
) {
	/*calculates maximum absolute velocities in x and y direction*/
	double umax=0, vmax=0;
	double a,b,c;
	int i, j;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j<=jmax; j++) {
			if((Flag[i][j]&B_C)==B_C){
				if(abs(U[i][j])>umax)
					umax = abs(U[i][j]);

				if(abs(V[i][j])>vmax)
					vmax = abs(V[i][j]);

			}
		}
	}

	/*Determines the minimum of dt according to stability criteria and multiply it by safety factor tau if tau is positive, otherwise uses the default value of dt*/
	if (tau>0){
		a = Re/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
		b = dx/umax;
		c = dy/vmax;
		if(a < b && a < c){
			*dt = tau * a;
		}
		else if(b < a && b < c){
			*dt = tau * b;
		}
		else{
			*dt = tau * c;
		}
	}
}




void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P, int **Flag) {

	int i;
	int j;
	for(i = 1; i <= imax; i++)
	for(j = 1; j <= jmax; j++)
	{
			/*
			 * Check that the cell is a fluid cell.
			 */
			if((Flag[i][j]&B_C)==B_C){
				/*Calculate the new velocity U according to the formula above*/
				if(i<imax){
					U[i][j] = F[i][j]-(dt/dx)*(P[i+1][j]-P[i][j]);
				}
				/*Calculate the new velocity V according to the formula above*/
				if(j<jmax){
					V[i][j] = G[i][j]-(dt/dy)*(P[i][j+1]-P[i][j]);
				}
			}
	}
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS, int **Flag) {
	/* right hand side of formula 11 */
	int i;
	int j;

	for(i = 1; i <= imax; i++)
	for(j = 1; j <= jmax; j++)
	{
			RS[i][j] = 1 / dt*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy);
	}
}

void calculate_fg1(
		double Re,
		double GX,
		double GY,
		double alpha,
		double dt,
		double dx,
		double dy,
		int imax,
		int jmax,
		double **U,
		double **V,
		double **F,
		double **G,
		int **Flag
)

{
	int i ;
	int j ;
	double du2dx ;
	double duvdy ;
	double d2udx2 ;
	double d2udy2 ;

	double dv2dy ;
	double duvdx ;
	double d2vdx2 ;
	double d2vdy2 ;

	/*Determines the value of F according to the formula above with the help of temporary variables*/
	for ( i = 1 ; i <= imax ; i++ )
	{
		for( j = 1 ; j <= jmax ; j++ )
		{
			/*
			 * We need to check that the cell is actually a fluid cell.
			 */
			if(((Flag[i][j]&B_C)==B_C)&& i<imax ){
				d2udx2 = ( U[i+1][j]  - 2*U[i][j] + U[i-1][j] ) / ( dx * dx) ;

				d2udy2 = ( U[i][j+1]  - 2*U[i][j] + U[i][j-1]) / (dy * dy )  ;

				du2dx = (1/dx) * ( ( (U[i][j] + U[i+1][j])/2 )*( (U[i][j] + U[i+1][j])/2 ) - ( (U[i-1][j] + U[i][j])/2 )*( (U[i-1][j] + U[i][j])/2 ) ) +
						alpha/dx * ( abs( U[i][j] + U[i+1][j] ) / 2  * ( U[i][j] - U[i+1][j] ) / 2 - abs( U[i-1][j] + U[i][j] ) / 2  * ( U[i-1][j] - U[i][j] ) / 2   ) ;

				duvdy = (1/dy) * ( ( V[i][j] + V[i+1][j] ) /2  *  ( U[i][j] + U[i][j+1] )/2 - (V[i][j-1] + V[i+1][j-1])/2 * (U[i][j-1] + U[i][j])/2  ) +
						alpha/dy * (abs( V[i][j] + V[i+1][j] ) /2  *  ( U[i][j] - U[i][j+1] )/2 - abs(V[i][j-1] + V[i+1][j-1])/2 * (U[i][j-1] - U[i][j])/2 ) ;

				F[i][j] = U[i][j]  + dt * ( 1/Re * ( (d2udx2 ) + (d2udy2) ) - (du2dx)  - duvdy + GX ) ;

			}
			/*Determines the value of G according to the formula above with the help of temporary variables*/
			if(((Flag[i][j]&B_C)==B_C) && j<jmax ){
				d2vdx2 = ( V[i+1][j]  - 2*V[i][j] + V[i-1][j] ) / ( dx * dx) ;

				d2vdy2 = ( V[i][j+1]  - 2*V[i][j] + V[i][j-1]) / (dy * dy )  ;

				duvdx = (1/dx) * ( ( V[i][j] + V[i+1][j] ) /2  *  ( U[i][j] + U[i][j+1] )/2 - (U[i-1][j] + U[i-1][j+1])/2 * (V[i-1][j] + V[i][j])/2  ) +
						alpha/dx * (( V[i][j] - V[i+1][j] ) /2  *  abs( U[i][j] + U[i][j+1] )/2 - abs(U[i-1][j] + U[i-1][j+1])/2 * (V[i-1][j] - V[i][j])/2 ) ;

				dv2dy = (1/dy) * ( ( (V[i][j] + V[i][j+1])/2 )*( (V[i][j] + V[i][j+1])/2 ) - ( (V[i][j-1] + V[i][j])/2 )* (V[i][j-1] + V[i][j])/2 )  +
						alpha/dy * ( abs( V[i][j] + V[i][j+1] ) / 2  * ( V[i][j] - V[i][j+1] ) / 2 - abs( V[i][j-1] + V[i][j] ) / 2  * (  V[i][j-1] - V[i][j]  ) / 2   ) ;

				G[i][j] = V[i][j]  + dt * ( 1/Re * ( (d2vdx2 ) + (d2vdy2) ) - (duvdx)  - dv2dy + GY ) ;
			}
			/*
			 * In case its a boundary cell, then we check it by comparing the flags and calculate
			 * only the useful values of F and G.
			 */
			if((Flag[i][j]&31)==B_N){
				G[i][j]=V[i][j];
			}
			else if((Flag[i][j]&31)==B_S){
				G[i][j-1]=V[i][j-1];
			}
			else if((Flag[i][j]&31)==B_W){
				F[i-1][j]=U[i-1][j];
			}
			else if((Flag[i][j]&31)==B_O){
				F[i][j]=U[i][j];
			}
			else if((Flag[i][j]&31)==B_NO){
				F[i][j]=U[i][j];
				G[i][j]=V[i][j];
			}
			else if((Flag[i][j]&31)==B_NW){
				F[i-1][j]=U[i-1][j];
				G[i][j]=V[i][j];
			}
			else if((Flag[i][j]&31)==B_SO){
				F[i][j]=U[i][j];
				G[i][j-1]=V[i][j-1];
			}
			else if((Flag[i][j]&31)==B_SW){
				F[i-1][j]=U[i-1][j];
				G[i][j-1]=V[i][j-1];
			}
		}
	}


	/*Set boundary values along the columns*/
	for (j = 1; j <= jmax; j++){
		/*F values on right and left boundaries*/
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}

	/*Set boundary values along the rows*/
	for (i = 1; i <= imax; i++){
		/*G values on top and bottom boundaries*/
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
}
