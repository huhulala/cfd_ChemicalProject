#include "boundary_val.h"
#include "NSDefinitions.h"
#include <string.h>

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax, int jmax, double dx, double dy, int wl, int wr,
		int wt, int wb, double **U, double **V, double **F, double **G,
		double **P, double **T, int** Flag, double ***C,int s_max,
		  double cl,
		  double cr,
		  double cb,
		  double ct) {

	int i = 0;
	int j = 0;
	int k = 0;

	for (i = 1; i <= imax; i++) {
		/* lower bounder */
		switch (wb) {
		case NO_SLIP:
			U[i][0] = -U[i][1];
			V[i][0] = 0.0;
			break;
		case FREE_SLIP:
			U[i][0] = U[i][1];
			V[i][0] = 0.0;
			break;
		case OUTFLOW:
			U[i][0] = U[i][1];
			V[i][0] = V[i][1];
			break;
		}

		/* upper bounder */
		switch (wt) {
		case NO_SLIP:
			U[i][jmax + 1] = -U[i][jmax];
			V[i][jmax] = 0.0;
			break;
		case FREE_SLIP:
			U[i][jmax + 1] = U[i][jmax];
			V[i][jmax] = 0.0;
			break;
		case OUTFLOW:
			U[i][jmax + 1] = U[i][jmax];
			V[i][jmax] = V[i][jmax - 1];
			break;
		}

		T[i][jmax+1] = T[i][jmax];
	}
	for (j = 1; j <= jmax; j++) {
		/* left border */
		switch (wl) {
		case NO_SLIP:
			U[0][j] = 0.0;
			V[0][j] = -V[1][j];
			break;
		case FREE_SLIP:
			U[0][j] = 0.0;
			V[0][j] = V[1][j];
			break;
		case OUTFLOW:
			U[0][j] = U[1][j];
			V[0][j] = V[1][j];
			break;
		}

		/* right border */
		switch (wr) {
		case NO_SLIP:
			U[imax][j] = 0.0;
			V[imax + 1][j] = -V[imax][j];
			break;
		case FREE_SLIP:
			U[imax][j] = 0.0;
			V[imax + 1][j] = V[imax][j];
			break;
		case OUTFLOW:
			U[imax][j] = U[imax - 1][j];
			V[imax + 1][j] = V[imax][j];
			break;
		}

		T[imax+1][j] = T[imax][j];
	}

	/* set the values of the inner obstacles */
	/* treat obstacle boundaries, see 1.4 to 1.6 */
	for (j = 1; j <= jmax; j++)
		for (i = 1; i <= imax; i++) {
			switch (Flag[i][j]) {
			case B_N:
				V[i][j] = 0;
				U[i - 1][j] = -U[i - 1][j + 1];
				U[i][j] = -U[i][j + 1];
				break;
			case B_S:
				V[i][j - 1] = 0;
				U[i - 1][j] = -U[i - 1][j - 1];
				U[i][j] = -U[i][j - 1];
				break;
			case B_W:
				U[i - 1][j] = 0;
				V[i][j - 1] = -V[i - 1][j - 1];
				V[i][j] = -V[i - 1][j];
				break;
			case B_O:
				U[i][j] = 0;
				V[i][j - 1] = -V[i + 1][j - 1];
				V[i][j] = -V[i + 1][j];
				break;
			case B_NO:
				U[i][j] = 0;
				U[i - 1][j] = -U[i - 1][j + 1];
				V[i][j] = 0;
				V[i][j - 1] = -V[i + 1][j - 1];
				break;
			case B_NW:
				U[i - 1][j] = 0;
				U[i][j] = -U[i][j + 1];
				V[i][j] = 0;
				V[i][j - 1] = -V[i - 1][j - 1];
				break;
			case B_SO:
				U[i][j] = 0;
				U[i - 1][j] = -U[i - 1][j - 1];
				V[i][j - 1] = 0;
				V[i][j] = -V[i + 1][j];
				break;
			case B_SW:
				U[i - 1][j] = 0;
				U[i][j] = -U[i][j - 1];
				V[i][j - 1] = 0;
				V[i][j] = -V[i - 1][j];
				break;
			}
		}


	for (i = 1; i <= imax; ++i) {
		for (j = 1; j <= jmax; ++j) {
			if (Flag[i][j] == B_N) {
				T[i][j] = T[i][j + 1];
			}

			else if (Flag[i][j] == B_W) {
				T[i][j] = T[i - 1][j];
			}

			else if (Flag[i][j] == B_S) {
				T[i][j] = T[i][j - 1];
			}

			else if (Flag[i][j] == B_O) {
				T[i][j] = T[i + 1][j];
			}

			else if (Flag[i][j] == B_NO) {
				T[i][j] = (T[i + 1][j] + T[i][j + 1]) / 2;
			}

			else if (Flag[i][j] == B_NW) {
				T[i][j] = (T[i - 1][j] + T[i][j + 1]) / 2;
			}

			else if (Flag[i][j] == B_SO) {
				T[i][j] = (T[i + 1][j] + T[i][j - 1]) / 2;
			}

			else if (Flag[i][j] == B_SW) {
				T[i][j] = (T[i - 1][j] + T[i][j - 1]) / 2;
			}
		}
	}

	 /** concentration goes here;
	  * either boundary cons or inflow **/
    for (k = 0; k < s_max; ++k)
    {
        /** left and right wall **/
        for (j = 0; j <= jmax; ++j)
        {
         	if(cl == 0.0)
         	  C[k][0][j] = -C[k][1][j];
         	else if(cl > 0.0)
         	  C[k][0][j] = 2*cl - C[k][1][j];

         	if(cr == 0.0)
         	  C[k][imax+1][j] = C[k][imax][j];
         	else if(cr > 0.0)
         	  C[k][imax+1][j] = 2*cr-C[k][imax][j];
        }
        /** top and bottom wall **/
        for(i = 0; i <= imax; ++i)
        {
        	if(cb == 0.0)
        	  C[k][i][0] = C[k][i][1];
        	else if(cb > 0.0)
              C[k][i][0] = 2*cb - C[k][i][1];

        	if(ct == 0.0)
              C[k][i][jmax+1] = C[k][i][jmax];
        	else if(ct > 0.0)
              C[k][i][0] = 2*ct - C[k][i][jmax];
        }
    }
}

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues1(int imax, int jmax, double **U, double **V, const int wl,
		const int wr, const int wt, const int wb, int **Flag) {

	int i, j;
	/*Initialize corners*/
	U[0][0] = 0.0;
	U[0][jmax + 1] = 0.0;
	U[imax + 1][0] = 0.0;
	U[imax + 1][jmax + 1] = 0.0;
	V[0][0] = 0.0;
	V[0][jmax + 1] = 0.0;
	V[imax + 1][0] = 0.0;
	V[imax + 1][jmax + 1] = 0.0;

	/* Set values for all the outside boundary depending on the value that
	 * the data file is set to. NO_SLIP = 1, FREE_SLIP=2 and OUTFLOW=3
	 * We start with the left boundary. */
	switch (wl) {
	case NO_SLIP:
		for (j = 1; j < jmax + 1; j++) {
			/*U velocity on left boundary */
			U[0][j] = 0;
			/*V velocity left boundary */
			V[0][j] = -1 * V[1][j];
		}
		break;
	case FREE_SLIP:
		for (j = 1; j < jmax + 1; j++) {
			/*U velocity on left boundary */
			U[0][j] = 0;
			/*V velocity left boundary */
			V[0][j] = V[1][j];
		}
		break;
	case OUTFLOW:
		for (j = 1; j < jmax + 1; j++) {
			/*U velocity on left boundary */
			U[0][j] = U[1][j];
			/*V velocity left boundary */
			V[0][j] = V[1][j];
		}
		break;
	default:
		break;
	}

	/*Set values for the right boundary*/
	switch (wr) {
	case NO_SLIP:
		for (j = 1; j < jmax + 1; j++) {
			U[imax][j] = 0;
			V[imax + 1][j] = -1 * V[imax][j];
		}
		break;
	case FREE_SLIP:
		for (j = 1; j < jmax + 1; j++) {
			U[imax][j] = 0;
			V[imax + 1][j] = V[imax][j];
		}
		break;
	case OUTFLOW:
		for (j = 1; j < jmax + 1; j++) {
			U[imax][j] = U[imax - 1][j];
			V[imax + 1][j] = V[imax][j];
		}
		break;
	default:
		break;
	}

	/*Set values for the top boundary*/
	switch (wt) {
	case NO_SLIP:
		for (i = 1; i < imax + 1; i++) {
			V[i][jmax] = 0;
			U[i][jmax + 1] = -1 * U[i][jmax];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i < imax + 1; i++) {
			V[i][jmax] = 0;
			U[i][jmax + 1] = U[i][jmax];
		}
		break;
	case OUTFLOW:
		for (i = 1; i < imax + 1; i++) {
			V[i][jmax] = V[i][jmax - 1];
			U[i][jmax + 1] = U[i][jmax];
		}
		break;
	default:
		break;
	}

	/*Set values for the bottom boundary*/
	switch (wb) {
	case NO_SLIP:
		for (i = 1; i < imax + 1; i++) {
			V[i][0] = 0;
			U[i][0] = -1 * U[i][1];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i < imax + 1; i++) {
			V[i][0] = 0;
			U[i][0] = U[i][1];
		}
		break;
	case OUTFLOW:
		for (i = 1; i < imax + 1; i++) {
			V[i][0] = V[i][1];
			U[i][0] = U[i][1];
		}
		break;
	default:
		break;
	}

	/**
	 * Loop to check for boundary cells in the inner domain, and assign
	 * correct values of U and V.
	 */
	for (i = 1; i < imax + 1; i++) {
		for (j = 1; j < jmax + 1; j++) {
			if ((Flag[i][j] & 31) == B_N) {
				V[i][j] = 0;
				U[i - 1][j] = -1 * U[i - 1][j + 1];
				U[i][j] = -1 * U[i][j + 1];
			} else if ((Flag[i][j] & 31) == B_S) {
				V[i][j - 1] = 0;
				U[i - 1][j] = -1 * U[i - 1][j - 1];
				U[i][j] = -1 * U[i][j - 1];
			} else if ((Flag[i][j] & 31) == B_W) {
				U[i - 1][j] = 0;
				V[i][j - 1] = -1 * V[i - 1][j - 1];
				V[i][j] = -1 * V[i - 1][j];
			} else if ((Flag[i][j] & 31) == B_O) {
				U[i][j] = 0;
				V[i][j - 1] = -1 * V[i + 1][j - 1];
				V[i][j] = -1 * V[i + 1][j];
			} else if ((Flag[i][j] & 31) == B_NO) {
				U[i][j] = 0;
				U[i - 1][j] = -1 * U[i - 1][j + 1];
				V[i][j] = 0;
				V[i][j - 1] = -1 * V[i + 1][j - 1];
			} else if ((Flag[i][j] & 31) == B_NW) {
				U[i - 1][j] = 0;
				U[i][j] = -1 * U[i][j + 1];
				V[i][j] = 0;
				V[i][j - 1] = -1 * V[i - 1][j - 1];
			} else if ((Flag[i][j] & 31) == B_SO) {
				U[i][j] = 0;
				U[i - 1][j] = -1 * U[i - 1][j - 1];
				V[i][j - 1] = 0;
				V[i][j] = -1 * V[i + 1][j];
			} else if ((Flag[i][j] & 31) == B_SW) {
				U[i - 1][j] = 0;
				U[i][j] = -1 * U[i][j - 1];
				V[i][j] = -1 * V[i - 1][j];
				V[i][j - 1] = 0;
			}
		}
	}

}

void spec_boundary_val(char *problem, int imax, int jmax, int s_max, double dx, double dy,
		double Re, double deltaP, double **U, double **V, double **P, double **T,
		double ***C) {
	int s;
	int j = 0;
	if (strcmp(problem, "karman") == 0 || strcmp(problem, "baffle") == 0) {
		for (j = 1; j <= jmax; ++j) {
			/* set karman inflow */
			U[0][j] = 1.0;
			V[0][j] = 0.0;
		}
	}
	else if (strcmp(problem, "karman_diffusion") == 0) {
		for (j = 1; j <= jmax; ++j) {
			/* set karman inflow */
			U[0][j] = 1.0;
			V[0][j] = 0.0;
		}
	}
	else if (strcmp(problem, "step") == 0)
	{
		/* set step inflow - only upper half */
		for (j = (jmax + 2) / 2 + 1; j < jmax + 1; ++j) {
			U[0][j] = 1.0;
			V[0][j] = 0.0;
		}
		for (j = 1; j < (jmax + 2) / 2; ++j) {
			U[0][j] = 0.0;
			V[0][j] = 0.0;
		}
	}
	//  Rayleigh-Benard flow: top T = -0.5 bottom T = 0.5
	//  left and right walls adiabatic,
	//  lower wall heated, upper wall cooled.
	//
	 else if ( strcmp(problem, "rayleigh") ==0)
	 {
	     for(j=0; j<=jmax+1; j++)
	     {
	    	 T[0][j] = T[1][j];
	    	 T[imax+1][j] = T[imax][j];
	      }

	     for(j=0;j<=imax+1;j++)
	     {
	    	 T[j][0] = 2*(0.5)-T[j][1];
	    	 T[j][jmax+1] = 2*(-0.5)-T[j][jmax];
	     }
	     return;
	  }
	 else if (strcmp(problem, "fluidTrap")==0)
	 {
	     for(j=0; j<=jmax+1; j++)
	     {
	    	 T[0][j] = 2*(2.5)-T[1][j];
	    	 T[imax+1][j] = 2*(-2.5)-T[imax][j];
	      }

	     for(j=0;j<=imax+1;j++)
	     {
	    	 T[j][0] = T[j][1];
	    	 T[j][jmax+1] = T[j][jmax] ;
	     }
	     return;
	  }
	 else if (strcmp(problem, "rayleigh_plane")==0)
	 {
	     for(j=0; j<=jmax+1; j++)
	     {
	    	 T[0][j] = T[1][j];
	    	 T[imax+1][j] = T[imax][j];
	      }

	     for(j=0;j<=imax+1;j++)
	     {
	    	 T[j][0] = 2*(293.5)-T[j][1];
	    	 T[j][jmax+1] = 2*(292.5)-T[j][jmax];
	     }
	     return;
	  }
}

