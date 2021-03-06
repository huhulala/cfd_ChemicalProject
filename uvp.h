#ifndef __UVP_H__
#define __UVP_H__


/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
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
  int **Flag,
  double **T,
  double beta
  );


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS,
  int **Flag
);


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double Pr,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  int s_max,
  double **U,
  double **V,
  double ***C,
  double lambda
);


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  int **Flag
);

/* Computes the new temperature values */
void calculate_Temp(double **U, double **V, double **T, int **Flag, int imax, int jmax,
		double dt, double dx, double dy, double gamma,
		double Re, double Pr, double**H);

/* Calculates the concentrations of each substance */
void calculate_Concentrations(double **U, double **V, double ***C, double ***Q,
		int **Flag, int imax, int jmax, int s_max, double dt, double dx, double dy,
		double lambda, double gamma2);

/* Computes the changement of temperature and concentration for each cell caused by
 * an irreversible chemical reaction
 * @param ***C	Array of concentrations
 * @param ***Q	Change of concentrations
 * @param **H	By reaction consumed/generated heat
 * @param imax	Number of cells in x-direction
 * @param jmax	Number of cells in y-direction
 * @param s_max	Number of substances
 * @param a,b,c,d	Coefficients of the reaction aA + bB -> cC + dD
 * @param dH	Enthalpy of reaction, heat that gets free/consumed by one reaction */
void chemical_reaction_irreversible(double ***C, double ***Q, double **H, int imax, int jmax, int s_max,
		int a, int b, int c, int d, double dH);


/* Computes the changement of temperature and concentration for each cell caused by
 * an reversible chemical reaction
 * @param ***C	Array of concentrations, 0 and 1 store reactant concentrations, 2 and 3 product concentrations
 * @param ***Q	Change of concentrations
 * @param **H	By reaction consumed/generated heat
 * @param imax	Number of cells in x-direction
 * @param jmax	Number of cells in y-direction
 * @param s_max	Number of substances (always 4)
 * @param a,b,c,d	Coefficients of the reaction aA + bB -> cC + dD
 * @param k1	Speed coefficient for the left-to-right reaction
 * @param k2	Speed coefficient for the right-to-left reaction
 * @param dH	Enthalpy of reaction, heat that gets free/consumed by one reaction
 * @param dt	Time step size */
void chemical_reaction_reversible(double ***C, double ***Q, double **H, int imax, int jmax, int s_max,
		int a, int b, int c, int d, double k1, double k2, double dH, double dt);

#endif
