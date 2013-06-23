#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int wl,
  int wr,
  int wt,
  int wb,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  int **Flag
);

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
);
/**

/*Set special boundary*/

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    double dx,
    double dy,
    double Re,
    double deltaP,
    double **U,
    double **V,
    double **P
);

#endif
