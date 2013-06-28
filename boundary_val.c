#include "boundary_val.h"
#include "NSDefinitions.h"
#include <string.h>

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax,int jmax,  double dx,
		  double dy, int wl, int wr, int wt, int wb, double **U, double **V,
	 double **F, double **G,double **P, double **T,int** Flag,double tl,double tr,double tb,double tt)
{
    int i = 0;
    int j = 0;

	for(i=1; i<=imax; i++)
	{
		/* lower bounder */
		switch(wb)
		{
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
		switch(wt)
		{
			case NO_SLIP:
				U[i][jmax+1] = -U[i][jmax];
				V[i][jmax] = 0.0;
				break;
			case FREE_SLIP:
				U[i][jmax+1] = U[i][jmax];
				V[i][jmax] = 0.0;
				break;
			case OUTFLOW:
				U[i][jmax+1] = U[i][jmax];
				V[i][jmax] = V[i][jmax-1];
				break;
		}
	}
	for(j=1; j<=jmax; j++)
	{
		/* left border */
		switch(wl)
		{
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
		switch(wr)
		{
			case NO_SLIP:
				U[imax][j] = 0.0;
				V[imax+1][j] = -V[imax][j];
				break;
			case FREE_SLIP:
				U[imax][j] = 0.0;
				V[imax+1][j] = V[imax][j];
				break;
			case OUTFLOW:
				U[imax][j] = U[imax-1][j];
				V[imax+1][j] = V[imax][j];
				break;
		}
	}

	/** left and right wall **/
	for ( j = 1; j <= jmax; ++j ) {
            /* left wall */
            if (tl > 0)
            {
                T[0][j] = 2*tl - T[1][j];
            }
            else
            {
            	double dTdn = 0;
                T[0][j] = T[1][j] + (dTdn * (j - 0.5) * dy * dx);
            }
            /* right wall */
            if (tr > 0)
            {
                T[imax+1][j] = 2*tr - T[imax][j];
            }
            else
            {
            	double dTdn = 0;
                T[imax+1][j] = T[imax][j] + (dTdn * (j - 0.5) * dy * dx);
            }
	}

	/** bottom and top wall **/
	for (i = 1; i <= imax; ++i) {
            /* bottom wall */
            if (tb > 0)
            {
            	T[i][0] = 2*tb - T[i][1];
            }
            else
            {
            	double dTdn = 0;
                T[i][0] = T[i][1] + (dTdn*(i - 0.5) *dx * dy);
            }

            /* top wall */
            if (tt > 0)
            {
            	T[i][jmax+1] = 2*tt - T[i][jmax];
            }
            else
            {
            	double dTdn = 0;
                T[i][jmax+1] = T[i][jmax] + (dTdn * (i-0.5)*dx*dy);
            }
	}

	/**** TEMPERATURE END ****/s


    /* set the values of the inner obstacles */
	/* treat obstacle boundaries, see 1.4 to 1.6 */
	for(j = 1; j <= jmax; j++)
	for(i = 1; i <= imax; i++)
	{
		switch(Flag[i][j])
		{
		case B_N:
			V[i][j] = 0;
			U[i-1][j] = -U[i-1][j+1];
			U[i][j] = -U[i][j+1];
			break;
		case B_S:
			V[i][j-1] = 0;
			U[i-1][j] = -U[i-1][j-1];
			U[i][j] = -U[i][j-1];
			break;
		case B_W:
			U[i-1][j] = 0;
			V[i][j-1] = -V[i-1][j-1];
			V[i][j] = -V[i-1][j];
			break;
		case B_O:
			U[i][j] = 0;
			V[i][j-1] = -V[i+1][j-1];
			V[i][j] = -V[i+1][j];
			break;
		case B_NO:
			U[i][j] = 0;
			U[i-1][j] = -U[i-1][j+1];
			V[i][j] = 0;
			V[i][j-1] = -V[i+1][j-1];
			break;
		case B_NW:
			U[i-1][j] = 0;
			U[i][j] = -U[i][j+1];
			V[i][j] = 0;
			V[i][j-1] = -V[i-1][j-1];
			break;
		case B_SO:
			U[i][j] = 0;
			U[i-1][j] = -U[i-1][j-1];
			V[i][j-1] = 0;
			V[i][j] = -V[i+1][j];
			break;
		case B_SW:
			U[i-1][j] = 0;
			U[i][j] = -U[i][j-1];
			V[i][j-1] = 0;
			V[i][j] = -V[i-1][j];
			break;
		}
	}
}

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues1(
		int imax,
		int jmax,
		double **U,
		double **V,
		const int wl,
		const int wr,
		const int wt,
		const int wb,
		int **Flag
) {

	int i,j;
	/*Initialize corners*/
	U[0][0]=0.0;
	U[0][jmax+1]=0.0;
	U[imax+1][0]=0.0;
	U[imax+1][jmax+1]=0.0;
	V[0][0]=0.0;
	V[0][jmax+1]=0.0;
	V[imax+1][0]=0.0;
	V[imax+1][jmax+1]=0.0;

	/* Set values for all the outside boundary depending on the value that
	 * the data file is set to. NO_SLIP = 1, FREE_SLIP=2 and OUTFLOW=3
	 * We start with the left boundary. */
	switch(wl){
	case NO_SLIP:
		for (j = 1; j < jmax + 1; j++){
			/*U velocity on left boundary */
			U[0][j] = 0;
			/*V velocity left boundary */
			V[0][j]=-1*V[1][j];
		}
		break;
	case FREE_SLIP:
		for (j = 1; j < jmax + 1; j++){
			/*U velocity on left boundary */
			U[0][j] = 0;
			/*V velocity left boundary */
			V[0][j] = V[1][j];
		}
		break;
	case OUTFLOW:
		for (j = 1; j < jmax + 1; j++){
			/*U velocity on left boundary */
			U[0][j] = U[1][j];
			/*V velocity left boundary */
			V[0][j]= V[1][j];
		}
		break;
	default:
		break;
	}

	/*Set values for the right boundary*/
	switch(wr){
	case NO_SLIP:
		for (j = 1; j < jmax + 1; j++){
			U[imax][j] = 0;
			V[imax+1][j]=-1*V[imax][j];
		}
		break;
	case FREE_SLIP:
		for (j = 1; j < jmax + 1; j++){
			U[imax][j] = 0;
			V[imax+1][j] = V[imax][j];
		}
		break;
	case OUTFLOW:
		for (j = 1; j < jmax + 1; j++){
			U[imax][j] = U[imax-1][j];
			V[imax+1][j]= V[imax][j];
		}
		break;
	default:
		break;
	}

	/*Set values for the top boundary*/
	switch(wt){
	case NO_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = 0;
			U[i][jmax+1]= -1*U[i][jmax];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = 0;
			U[i][jmax+1] = U[i][jmax];
		}
		break;
	case OUTFLOW:
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = V[i][jmax-1];
			U[i][jmax+1] = U[i][jmax];
		}
		break;
	default:
		break;
	}

	/*Set values for the bottom boundary*/
	switch(wb){
	case NO_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][0] = 0;
			U[i][0]=-1*U[i][1];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][0] = 0;
			U[i][0] = U[i][1];
		}
		break;
	case OUTFLOW:
		for (i = 1; i < imax + 1; i++){
			V[i][0] = V[i][1];
			U[i][0]= U[i][1];
		}
		break;
	default:
		break;
	}

	/**
	 * Loop to check for boundary cells in the inner domain, and assign
	 * correct values of U and V.
	 */
	for(i = 1; i < imax+1; i++){
		for(j = 1; j < jmax+1; j++){
			if((Flag[i][j]&31)==B_N){
				V[i][j]=0;
				U[i-1][j]=-1*U[i-1][j+1];
				U[i][j]=-1*U[i][j+1];
			}
			else if((Flag[i][j]&31)==B_S){
				V[i][j-1]=0;
				U[i-1][j]=-1*U[i-1][j-1];
				U[i][j]=-1*U[i][j-1];
			}
			else if((Flag[i][j]&31)==B_W){
				U[i-1][j]=0;
				V[i][j-1]=-1*V[i-1][j-1];
				V[i][j]=-1*V[i-1][j];
			}
			else if((Flag[i][j]&31)==B_O){
				U[i][j]=0;
				V[i][j-1]=-1*V[i+1][j-1];
				V[i][j]=-1*V[i+1][j];
			}
			else if((Flag[i][j]&31)==B_NO){
				U[i][j]=0;
				U[i-1][j]=-1*U[i-1][j+1];
				V[i][j]=0;
				V[i][j-1]=-1*V[i+1][j-1];
			}
			else if((Flag[i][j]&31)==B_NW){
				U[i-1][j]=0;
				U[i][j]=-1*U[i][j+1];
				V[i][j]=0;
				V[i][j-1]=-1*V[i-1][j-1];
			}
			else if((Flag[i][j]&31)==B_SO){
				U[i][j]=0;
				U[i-1][j]=-1*U[i-1][j-1];
				V[i][j-1]=0;
				V[i][j]=-1*V[i+1][j];
			}
			else if((Flag[i][j]&31)==B_SW){
				U[i-1][j]=0;
				U[i][j]=-1*U[i][j-1];
				V[i][j]=-1*V[i-1][j];
				V[i][j-1]=0;
			}
		}
	}

}


void spec_boundary_val(char *problem,int imax, int jmax, double dx, double dy,
    double Re, double deltaP, double **U, double **V, double **P )
{
    int j = 0;
    if (strcmp(problem, "karman") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
        	/* set karman inflow */
            U[0][j] = 1.0;
            V[0][j] = 0.0;
        }
    }
	else if(strcmp(problem, "step")==0)
	{
		/* set step inflow - only upper half */
		for(j=(jmax+2)/2+1; j<jmax+1; ++j)
		{
			U[0][j] = 1.0;
			V[0][j] = 0.0;
		}
		for(j=1; j<(jmax+2)/2; ++j)
		{
			U[0][j] = 0.0;
			V[0][j] = 0.0;
		}

	}
}

