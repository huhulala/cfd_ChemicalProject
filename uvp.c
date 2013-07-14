#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include "fd.h"
#include <math.h>
#include <limits.h>
#include "NSDefinitions.h"

void calculate_dt(double Re, double Pr, double tau, double *dt, double dx,
		double dy, int imax, int jmax, int s_max, double **U, double **V, double ***C, double lambda) {
	// See formula 13
	double umax = fabs(U[1][1]);
	double vmax = fabs(V[1][1]);
	double dtcon, dxcon, dycon, dttcon, dccon;
	double min;
	int i, j, s;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		if (fabs(U[i][j]) > umax)
			umax = fabs(U[i][j]);
		if (fabs(V[i][j]) > vmax)
			vmax = fabs(V[i][j]);
	}
	// Compute maximal concentration difference



	// conditions
	dccon = lambda/(2*(1/(dx*dx) + 1/(dy*dy)));
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
	if (min > dxcon && dxcon!= 0)
		min = dxcon;
	if (min > dycon  && dycon!= 0)
		min = dycon;
	if (min > dttcon  && dxcon!= 0)
		min = dttcon;
	if (min > dccon  && dccon!= 0)
		min = dccon;

	// calculate dt
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

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
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

void calculate_Temp(double **U, double **V, double **T, int **Flag,
		int imax, int jmax, double dt, double dx, double dy, double gamma, double Re, double Pr, double**H) {

	double indelx2;
	double indely2;

	double laplt, dutdx, dvtdy;
	int i, j;

	laplt = 0.0;
	dutdx = 0.0;
	dvtdy = 0.0;

	indelx2 = 1 / (dx * dx);
	indely2 = 1 / (dy * dy);

	//double q; /* heat source */
	/* TODO: Heat Source bei reactions setzen */
	//q = 0;
	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			if (Flag[i][j] >= C_F) {


				laplt = (T[i + 1][j] - 2.0 * T[i][j] + T[i - 1][j])
						* indelx2 + (T[i][j + 1] - 2.0 * T[i][j]
						+ T[i][j - 1]) * indely2;

				dutdx = ((U[i][j] * 0.5 * (T[i][j] + T[i + 1][j]) - U[i
						- 1][j] * 0.5 * (T[i - 1][j] + T[i][j])) + gamma
						* (fabs(U[i][j]) * 0.5 * (T[i][j] - T[i + 1][j])
								- fabs(U[i - 1][j]) * 0.5 * (T[i - 1][j]
										- T[i][j]))) / dx;
				dvtdy = ((V[i][j] * 0.5 * (T[i][j] + T[i][j + 1])
						- V[i][j - 1] * 0.5 * (T[i][j - 1] + T[i][j]))
						+ gamma * (fabs(V[i][j]) * 0.5 * (T[i][j]
								- T[i][j + 1]) - fabs(V[i][j - 1]) * 0.5
								* (T[i][j - 1] - T[i][j]))) / dy;

				T[i][j] = T[i][j] + dt
						* (laplt / (Re * Pr) - dutdx - dvtdy + H[i][j]);

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
			 * C[0] and C[1] store reactants, C[2] and C[3] products concentrations */
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


void chemical_reaction_reversible(double ***C, double ***Q, double **H, int imax, int jmax, int s_max,
		int a, int b, int c, int d, double k1, double k2, double dH, double dt)
{
	int i,j,k;
	double dConcentration;
	double h;/* Euler step size */
	int euler_it = 5;/* Euler iterations */

	h = dt/euler_it;
	/* Solve dA/dt = dB/dt = -k1 * (c_A)^a * (c_B)^b
	 * and   dC/dt = dD/dt =  k1 * (c_A)^a * (c_B)^b */
	for(j=1; j<=jmax; j++)
		for(i=1; i<=imax; i++)
			if(C[0][i][j] != 0 && C[1][i][j] != 0)
			{
				Q[0][i][j] = C[0][i][j];
				Q[1][i][j] = C[1][i][j];
				Q[2][i][j] = C[2][i][j];
				Q[3][i][j] = C[3][i][j];

				/* Euler iterations */
				for(k=0; k<euler_it; k++)
				{
					/* left-to-right reaction aA + bB -> cC + dD
					 * (irreversible reactions use only this one) */
					dConcentration = h * k1 * pow(Q[0][i][j], a) * pow(Q[1][i][j], b);
					Q[0][i][j] -= dConcentration;
					Q[1][i][j] -= dConcentration;
					Q[2][i][j] += dConcentration;
					Q[3][i][j] += dConcentration;



					/* right-to-left reaction aA + bB <- cC + dD
					 * (only for reversible reactions, irreversible reactions have k2 = 0 */
					dConcentration = h * k2 * pow(Q[2][i][j], c) * pow(Q[3][i][j], d);
					Q[0][i][j] += dConcentration;
					Q[1][i][j] += dConcentration;
					Q[2][i][j] -= dConcentration;
					Q[3][i][j] -= dConcentration;
				}

				Q[0][i][j] -= C[0][i][j];
				Q[1][i][j] -= C[1][i][j];
				Q[2][i][j] -= C[2][i][j];
				Q[3][i][j] -= C[3][i][j];

				/* Energy */
				H[i][j] = dH *10* (-Q[0][i][j]);
				//printf("H[i][j]  %f \n",H[i][j]);
			}
			else /* no reaction */
			{
				for(k=0; k<s_max; k++)
					Q[k][i][j] = 0;
				H[i][j] = 0;
			}
}
