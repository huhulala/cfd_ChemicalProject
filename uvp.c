#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include "fd.h"
#include <math.h>
#include <limits.h>
#include "NSDefinitions.h"

void calculate_fg(double Re, double GX, double GY, double alpha, double beta,
		double dt, double dx, double dy, int imax, int jmax, double **U,
		double **V, double **F, double **G, double **T, int **Flag) {
	/* see formulas 9 and 10 in combination with formulas 4 and 5 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			/********** calculate F **********/
			/* calculate f and g only between two fluid cells */
			if (((Flag[i][j] & B_C) == B_C) && i < imax) {
				F[i][j] = U[i][j] + dt * (
				/* 1/Re * (d²u/dx² + d²u/dy²) */
				1 / Re * (d2udx2(i, j, U, dx) + d2udy2(i, j, U, dy))
				/* - du²/dx */
				- du2dx(i, j, U, dx, alpha)
				/* - duv/dy */
				- duvdy(i, j, U, V, dy, alpha) + GX)
				/* temperature */
				- beta * dt / 2 * (T[i][j] + T[i + 1][j]);
			}

			/********** calculate G **********/
			if (((Flag[i][j] & B_C) == B_C) && i < imax) {
				G[i][j] = V[i][j] + dt * (
				/* 1/Re * (d²v/dx² + d²v/dy²) */
				1 / Re * (d2vdx2(i, j, V, dx) + d2vdy2(i, j, V, dy))
				/* - duv/dx */
				- duvdx(i, j, U, V, dx, alpha)
				/* - dv²/dy */
				- dv2dy(i, j, V, dy, alpha) + GY)
				/* temperature */
				- beta * dt / 2 * (T[i][j] + T[i][j + 1]);
			}

			if ((Flag[i][j] & 31) == B_N) {
				G[i][j] = V[i][j];
			} else if ((Flag[i][j] & 31) == B_S) {
				G[i][j - 1] = V[i][j - 1];
			} else if ((Flag[i][j] & 31) == B_W) {
				F[i - 1][j] = U[i - 1][j];
			} else if ((Flag[i][j] & 31) == B_O) {
				F[i][j] = U[i][j];
			} else if ((Flag[i][j] & 31) == B_NO) {
				F[i][j] = U[i][j];
				G[i][j] = V[i][j];
			} else if ((Flag[i][j] & 31) == B_NW) {
				F[i - 1][j] = U[i - 1][j];
				G[i][j] = V[i][j];
			} else if ((Flag[i][j] & 31) == B_SO) {
				F[i][j] = U[i][j];
				G[i][j - 1] = V[i][j - 1];
			} else if ((Flag[i][j] & 31) == B_SW) {
				F[i - 1][j] = U[i - 1][j];
				G[i][j - 1] = V[i][j - 1];
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

/*
void calculate_dt(double Re, double Pr, double tau, double *dt, double dx,
		double dy, int imax, int jmax, int s_max, double **U, double **V, double ***C, double lambda) {
	// See formula 13
	double umax = fabs(U[1][1]);
	double vmax = fabs(V[1][1]);
	double dtcon, dxcon, dycon, dttcon, dccon;
	double min;
	int i, j, s;

	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			if (fabs(U[i][j]) > umax)
				umax = fabs(U[i][j]);
			if (fabs(V[i][j]) > vmax)
				vmax = fabs(V[i][j]);
		}
	// Compute maximal concentration difference
	dccon = fabs(C[0][1][1] - C[0][1][2]);
	for(s = 0; s < s_max; s++)
		for(j = 1; j<= jmax; j++)
			for(i = 1; i <= imax; i++)
			{
				if(fabs(C[s][i][j] - C[s][i][j+1]) > dccon)
					dccon = fabs(C[s][i][j] - C[s][i][j+1]);
				if(fabs(C[s][i][j] - C[s][i+1][j]) > dccon)
					dccon = fabs(C[s][i][j] - C[s][i+1][j]);
			}
	dccon = dx*(1-dccon);

	// conditions
	dtcon = Pr * Re / (2 * (1 / (dx * dx) + 1 / (dy * dy)));
    dttcon = Re/(2*(1/(dx*dx) + 1/(dy*dy)));

    if(umax != 0 && vmax != 0)
    {
		dxcon = dx / fabs(umax);
		dycon = dy / fabs(vmax);
    }
    else
    {
    	// values that should never be used
    	dxcon = LONG_MAX;
    	dycon = LONG_MAX;
    }

	// determine smalles condition
	min = dtcon;
	if (min > dxcon)
		min = dxcon;
	if (min > dycon)
		min = dycon;
	if (min > dttcon)
		min = dttcon;
	if (min > dccon)
		min = dccon;

	// calculate dt
	*dt = tau * min;
}
*/

void calculate_dt(double Re, double Pr, double tau, double *dt, double dx,
		double dy, int imax, int jmax, int s_max, double **U, double **V, double ***C, double lambda) {
	/* See formula 13 */
	double umax = fabs(U[1][1]);
	double vmax = fabs(V[1][1]);
	double dtcon, dxcon, dycon, dttcon, dtdcond;
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
	dtdcond = lambda/(2*(1/(dx*dx) + 1/(dy*dy)));
	dtcon = Pr * Re / (2 * (1 / (dx * dx) + 1 / (dy * dy)));
	dttcon = Re/(2*(1/(dx*dx) + 1/(dy*dy)));

	dxcon = dx / fabs(umax);
	dycon = dy / fabs(vmax);

	/* determine smalles condition */
	min = dtcon;
    if (min > dxcon && dxcon > 0)
		min = dxcon;
	if (min > dycon && dycon > 0)
		min = dycon;
	if (min > dttcon && dttcon > 0)
		min = dttcon;
	if (min > dtdcond && dtdcond > 0)
		min = dtdcond;

	/* calculate dt */
	*dt = tau * min;
}


void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P, int **Flag) {

	int i;
	int j;
	for (i = 1; i <= imax; i++)
		for (j = 1; j <= jmax; j++) {
			/*
			 * Check that the cell is a fluid cell.
			 */
			if ((Flag[i][j] & B_C) == B_C) {
				/*Calculate the new velocity U according to the formula above*/
				if (i < imax) {
					U[i][j] = F[i][j] - (dt / dx) * (P[i + 1][j] - P[i][j]);
				}
				/*Calculate the new velocity V according to the formula above*/
				if (j < jmax) {
					V[i][j] = G[i][j] - (dt / dy) * (P[i][j + 1] - P[i][j]);
				}
			}
		}
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS, int **Flag) {
	/* right hand side of formula 11 */
	int i;
	int j;

	for (i = 1; i <= imax; i++)
		for (j = 1; j <= jmax; j++) {
			RS[i][j] = 1 / dt * ((F[i][j] - F[i - 1][j]) / dx + (G[i][j]
					- G[i][j - 1]) / dy);
		}
}

void calculate_fg1(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G, int **Flag, double **T, double beta)

{
	int i;
	int j;
	double du2dx;
	double duvdy;
	double d2udx2;
	double d2udy2;

	double dv2dy;
	double duvdx;
	double d2vdx2;
	double d2vdy2;

	/*Determines the value of F according to the formula above with the help of temporary variables*/
	for (i = 1; i <= imax; i++) {
		for (j = 1; j <= jmax; j++) {
			/*
			 * We need to check that the cell is actually a fluid cell.
			 */
			if (((Flag[i][j] & B_C) == B_C) && i < imax) {

				d2udx2 = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (dx * dx);

				d2udy2 = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dy * dy);

				du2dx = (1 / dx) * (((U[i][j] + U[i + 1][j]) / 2) * ((U[i][j]
						+ U[i + 1][j]) / 2) - ((U[i - 1][j] + U[i][j]) / 2)
						* ((U[i - 1][j] + U[i][j]) / 2)) + alpha / dx * (abs(
						U[i][j] + U[i + 1][j]) / 2 * (U[i][j] - U[i + 1][j])
						/ 2 - abs(U[i - 1][j] + U[i][j]) / 2 * (U[i - 1][j]
						- U[i][j]) / 2);

				duvdy = (1 / dy) * ((V[i][j] + V[i + 1][j]) / 2 * (U[i][j]
						+ U[i][j + 1]) / 2 - (V[i][j - 1] + V[i + 1][j - 1])
						/ 2 * (U[i][j - 1] + U[i][j]) / 2) + alpha / dy * (abs(
						V[i][j] + V[i + 1][j]) / 2 * (U[i][j] - U[i][j + 1])
						/ 2 - abs(V[i][j - 1] + V[i + 1][j - 1]) / 2 * (U[i][j
						- 1] - U[i][j]) / 2);


				F[i][j] = U[i][j] + dt
						* ((1 / Re) * (d2udx2 + d2udy2) - du2dx - duvdy + GX)
						        -dt*beta*GX*(T[i][j]+T[i+1][j])/2;


			}
		    else
		    	 F[i][j] = U[i][j];

			/*Determines the value of G according to the formula above with the help of temporary variables*/
			if (((Flag[i][j] & B_C) == B_C) && j < jmax) {
				d2vdx2 = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (dx * dx);

				d2vdy2 = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dy * dy);

				duvdx = (1 / dx) * ((V[i][j] + V[i + 1][j]) / 2 * (U[i][j]
						+ U[i][j + 1]) / 2 - (U[i - 1][j] + U[i - 1][j + 1])
						/ 2 * (V[i - 1][j] + V[i][j]) / 2) + alpha / dx
						* ((V[i][j] - V[i + 1][j]) / 2 * abs(
								U[i][j] + U[i][j + 1]) / 2 - abs(
								U[i - 1][j] + U[i - 1][j + 1]) / 2
								* (V[i - 1][j] - V[i][j]) / 2);

				dv2dy = (1 / dy) * (((V[i][j] + V[i][j + 1]) / 2) * ((V[i][j]
						+ V[i][j + 1]) / 2) - ((V[i][j - 1] + V[i][j]) / 2)
						* (V[i][j - 1] + V[i][j]) / 2) + alpha / dy * (abs(
						V[i][j] + V[i][j + 1]) / 2 * (V[i][j] - V[i][j + 1])
						/ 2 - abs(V[i][j - 1] + V[i][j]) / 2 * (V[i][j - 1]
						- V[i][j]) / 2);

				double lala = dt*beta*GY*(T[i][j]+T[i][j+1])/2;


				G[i][j] = V[i][j] + dt
						* ((1 / Re) * (d2vdx2 + d2vdy2) - duvdx - dv2dy + GY)
					       -dt*beta*GY*(T[i][j]+T[i][j+1])/2;



			}
			else
				G[i][j] = V[i][j];
			/*
			 * In case its a boundary cell, then we check it by comparing the flags and calculate
			 * only the useful values of F and G.
			 */
			if ((Flag[i][j] & 31) == B_N) {
				G[i][j] = V[i][j];
			} else if ((Flag[i][j] & 31) == B_S) {
				G[i][j - 1] = V[i][j - 1];
			} else if ((Flag[i][j] & 31) == B_W) {
				F[i - 1][j] = U[i - 1][j];
			} else if ((Flag[i][j] & 31) == B_O) {
				F[i][j] = U[i][j];
			} else if ((Flag[i][j] & 31) == B_NO) {
				F[i][j] = U[i][j];
				G[i][j] = V[i][j];
			} else if ((Flag[i][j] & 31) == B_NW) {
				F[i - 1][j] = U[i - 1][j];
				G[i][j] = V[i][j];
			} else if ((Flag[i][j] & 31) == B_SO) {
				F[i][j] = U[i][j];
				G[i][j - 1] = V[i][j - 1];
			} else if ((Flag[i][j] & 31) == B_SW) {
				F[i - 1][j] = U[i - 1][j];
				G[i][j - 1] = V[i][j - 1];
			}
		}
	}

	/*Set boundary values along the columns*/
	for (j = 1; j <= jmax; j++) {
		/*F values on right and left boundaries*/
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}

	/*Set boundary values along the rows*/
	for (i = 1; i <= imax; i++) {
		/*G values on top and bottom boundaries*/
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
}

void calculate_Temp(double **U, double **V, double **TEMP, int **Flag,
		int imax, int jmax, double dt, double dx, double dy, double gamma, double Re, double Pr, double**H) {

	double dutdx;
	double d2tdx2;
	double dvtdy;
	double d2tdy2;
	double firstOperand;
	double secondOperand;
	double indelx2;
	double indely2;

	double LAPLT, DUTDX, DVTDY;
	int i, j;

	LAPLT = 0.0;
	DUTDX = 0.0;
	DVTDY = 0.0;

	indelx2 = 1 / (dx * dx);
	indely2 = 1 / (dy * dy);

	//double q; /* heat source */
	/* TODO: Heat Source bei reactions setzen */
	//q = 0;
	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			/* See formula 9.20
			 T[i][j] = T[i][j] + dt / (Re * Pr) * ((T[i+1][j] - 2*T[i][j] + T[i-1][j])/(dx*dx) +
			 (T[i][j+1] - 2*T[i][j] + T[i][j-1])/(dy*dy)) + q
			 -1/dx*(U[i][j]*(T[i][j]+T[i+1][j])/2 - U[i-1][j]*(T[i-1][j] + T[i][j])/2)
			 + gamma/dx*(fabs(U[i][j])*(T[i][j] - T[i+1][j])/2 - fabs(U[i-1][j])*(T[i-1][j] - T[i][j])/2)
			 - 1/dy*(V[i][j]*(T[i][j]+T[i][j+1])/2 - V[i][j-1]*(T[i][j-1] + T[i][j])/2)
			 + gamma/dy*(fabs(V[i][j])*(T[i][j] - T[i][j+1])/2 - fabs(V[i][j-1])*(T[i][j-1] - T[i][j])/2);
			 */

			if (Flag[i][j] >= C_F) {


				double lala1 = TEMP[i + 1][j];
				double lala2 = TEMP[i - 1][j];
				double lala3 = TEMP[i][j+1];
				double lala4 = TEMP[i][j-1];
				double lala5 = TEMP[i][j];


				double huhu = (TEMP[i + 1][j] - 2.0 * TEMP[i][j] + TEMP[i - 1][j]);
				double huhu1 = (TEMP[i][j + 1] - 2.0 * TEMP[i][j] + TEMP[i][j - 1]);


				LAPLT = (TEMP[i + 1][j] - 2.0 * TEMP[i][j] + TEMP[i - 1][j])
						* indelx2 + (TEMP[i][j + 1] - 2.0 * TEMP[i][j]
						+ TEMP[i][j - 1]) * indely2;

				DUTDX = ((U[i][j] * 0.5 * (TEMP[i][j] + TEMP[i + 1][j]) - U[i
						- 1][j] * 0.5 * (TEMP[i - 1][j] + TEMP[i][j])) + gamma
						* (fabs(U[i][j]) * 0.5 * (TEMP[i][j] - TEMP[i + 1][j])
								- fabs(U[i - 1][j]) * 0.5 * (TEMP[i - 1][j]
										- TEMP[i][j]))) / dx;
				DVTDY = ((V[i][j] * 0.5 * (TEMP[i][j] + TEMP[i][j + 1])
						- V[i][j - 1] * 0.5 * (TEMP[i][j - 1] + TEMP[i][j]))
						+ gamma * (fabs(V[i][j]) * 0.5 * (TEMP[i][j]
								- TEMP[i][j + 1]) - fabs(V[i][j - 1]) * 0.5
								* (TEMP[i][j - 1] - TEMP[i][j]))) / dy;

				double lala12 = TEMP[i][j] + dt
									* (LAPLT / (Re * Pr) - DUTDX - DVTDY + H[i][j]);
							TEMP[i][j] = lala12;
				TEMP[i][j] = lala12;



				/*
				 firstOperand = (1 / dx) * (U[i][j] * (T[i][j] + T[i + 1][j])
				 / 2 - U[i - 1][j] * (T[i - 1][j] + T[i][j]) / 2);
				 secondOperand = (alpha / dx) * (fabs(U[i][j]) * (T[i][j] - T[i
				 + 1][j]) / 2 - fabs(U[i - 1][j]) * (T[i - 1][j]
				 - T[i][j]) / 2);

				 dutdx = firstOperand + secondOperand;


				 firstOperand = (1 / dy) * (V[i][j] * (T[i][j] + T[i][j + 1])
				 / 2 - V[i][j - 1] * (T[i][j - 1] + T[i][j]) / 2);
				 secondOperand = (alpha / dy) * (fabs(V[i][j]) * (T[i][j]
				 - T[i][j + 1]) / 2 - fabs(V[i][j - 1]) * (T[i][j - 1]
				 - T[i][j]) / 2);

				 dvtdy = firstOperand + secondOperand;



				 d2tdx2 = (T[i + 1][j] - 2 * T[i][j] + T[i - 1][j]) / (dx * dx);



				 d2tdy2 = (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) / (dy * dy);



				 T[i][j] = T[i][j] + dt * ((1 / (Re * Pr)) * (d2tdx2 + d2tdy2)
				 - dutdx - dvtdy);
				 */
			}

		}
}

void calculate_Concentrations(double **U, double **V, double ***C, double ***Q,
		int **Flag, int imax, int jmax, int s_max, double dt, double dx, double dy,
		double lambda, double gamma2)
{
	/* equal to energy transfer equation */
	int i,j,s;
	//printf("C 11/8: %f, C 11/9: %f\n", C[0][11][8], C[0][11][9]);
	//printf("lambda: %f, dx: %f, dt: %f\n", lambda, dx, dt);
	for(s=0; s<s_max; s++)
		for(j=1; j<=jmax; j++)
			for(i=1; i<=imax; i++)
				if(Flag[i][j] >= C_F)
				{


					double ducdx = 1/dx*(U[i][j]*(C[s][i][j]+C[s][i+1][j])/2 - U[i-1][j]*(C[s][i-1][j] + C[s][i][j])/2)
										+ gamma2/dx*(fabs(U[i][j])*(C[s][i][j] -C[s][i+1][j])/2 - fabs(U[i-1][j])*(C[s][i-1][j] - C[s][i][j])/2);

					double dvcdy = 1/dy*(V[i][j]*(C[s][i][j]+C[s][i][j+1])/2 - V[i][j-1]*(C[s][i][j-1] + C[s][i][j])/2)
										+ gamma2/dy*(fabs(V[i][j])*(C[s][i][j] - C[s][i][j+1])/2 - fabs(V[i][j-1])*(C[s][i][j-1] - C[s][i][j])/2);

					C[s][i][j] = C[s][i][j] + dt *
						(lambda *
							(
								(C[s][i+1][j] - 2*C[s][i][j] + C[s][i-1][j])/(dx*dx) +
								(C[s][i][j+1] - 2*C[s][i][j] + C[s][i][j-1])/(dy*dy)
							)
							- ducdx
							- dvcdy
							+ Q[s][i][j]
						);

					/*if(C[s][i][j] < 0)
						C[s][i][j] = 0;
					if(C[s][i][j] > 1)
						C[s][i][j] = 1;*/
				}
}

void chemical_reaction_irreversible(double ***C, double ***Q, double **H, int imax, int jmax, int s_max,
		int a, int b, int c, int d, double dH)
{
	int i,j,s;
	double reacted_concentration;

	for(j=0; j<=jmax; j++)
		for(i=0; i<=imax; i++)
		{
			/* loop through each cell to see if there are substances to react
			 * C[0] and C[1] store reactant, C[2] and C[3] product concentrations */
			if(C[0][i][j] != 0 && C[1][i][j] != 0)
			{
				/* Volume is constant, compute the portion that reacted */
				reacted_concentration = C[0][i][j]/a;
				if(C[1][i][j]/b < reacted_concentration)
					reacted_concentration = C[1][i][j]/b;
				/* Change concentrations */
				Q[0][i][j] = -a * reacted_concentration;
				Q[1][i][j] = -b * reacted_concentration;
				Q[2][i][j] = c * reacted_concentration;
				Q[3][i][j] = d * reacted_concentration;
				/* Compute generated/consumed energy */
				H[i][j] = dH * reacted_concentration;
			}
			else
			{
				/* No reaction */
				for(s=0; s<s_max; s++)
					Q[s][i][j] = 0;
				H[i][j] = 0;
			}
		}
}
