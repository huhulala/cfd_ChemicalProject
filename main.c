#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args)
{
	/**************** Variable definition goes here  ****************/
	double Re, UI, VI, TI , PI, GX, GY;
	double t_end, xlength, ylength;
	double dt, dx, dy;
	double alpha, gamma, omg, tau;
	double eps, dt_value, res, t, deltaP;
	int itermax, it, n;
	double t_print;
	/**************** dimensionless quantities for temperature  ****************/
	double beta; /* cooefficient for termal expansion beta */
	double Pr;   /* Prandtl number Pr */

    double tl;
    double tr;
    double tb;
    double tt;

	int imax = 0;
	int jmax = 0;
	int wl = 0;
	int wr = 0;
	int wt = 0;
	int wb = 0;

	/* Output filename */
	char* output_filename;
	char output_filename_array[64];
	char* inputDir;
	char inputDirCharArray[64];
	char* problem;

	/* arrays */
	double **U = NULL;
	double **V = NULL;
	double **P = NULL;
	double **RS = NULL;
	double **F = NULL;
	double **G = NULL;

	/*  array for temperature  */
    double**  T = NULL;

	int **Problem = NULL;
	int **Flag = NULL;
	/*  flag field for the sources  */
    int **ChemicalSources = NULL;

	char inputString[64];
	char problemImageName[64];
	char szFileName[64];

	int verbose = 1; /* verbose flag */
	int debug = 1; /* verbose flag */
	/**************** Variable definition ends here  ****************/
	/* check arguments */
	if (argn <= 1)
	{
		printf("ERROR: you need to specify a problem \n");
		return 1;
	}
	if (!(strcmp(args[1], "karman") == 0
		|| strcmp(args[1], "plane") == 0
		|| strcmp(args[1], "step") == 0
		|| strcmp(args[1], "cavity") == 0
		|| strcmp(args[1], "rayleigh") == 0
		|| strcmp(args[1], "rayleigh_plane")== 0
		|| strcmp(args[1], "diffusion") == 0))
	{
		printf("ERROR: pass cavity, rayleigh, rayleigh_plane, fluidTrap, karman, plane or step\n");
		return 1;
	}

	/*************** parameter loading and problem input goes here *************************/

	t = 0.1;
	n = 0;
	problem = args[1];
	/* assemble parameter file string */
	output_filename= "./output_";
	inputDir= "./input/";

	strcpy(inputDirCharArray, inputDir);

	strcpy(output_filename_array, output_filename);
	strcat(output_filename_array, args[1]);
	strcat(output_filename_array, "/");
	strcat(output_filename_array, args[1]);
	strcpy(inputString, args[1]);
	strcpy(szFileName, inputString);
	strcat(szFileName, ".dat");
	strcat(inputDirCharArray, szFileName);

	/* load parameters from "problem".dat file */
	/* grid size (dx,dy,imax,ymax) should now be read from the image */
    read_parameters(inputDirCharArray,&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,&dt,&alpha,
    		        &omg,&tau,&itermax,&eps,&wl,&wr,&wt, &wb, &dt_value, &deltaP, &TI, &beta, &gamma,
    		        &Pr,&tl,&tr, &tb,&tt);

    /* assemble problem file string */
	strcpy(problemImageName, inputString);
	strcat(problemImageName, ".pgm");
	strcpy(inputDirCharArray, inputDir);
	strcat(inputDirCharArray, problemImageName);

	/*  load problem description from "problem".pgm file */
	/*  read imax and jmax (and calculate dx und dy) from problem description now */
	Problem = read_pgm(inputDirCharArray, &imax, &jmax);
	/* calculute dx/dy */
	dx = xlength / (double) (imax);
	dy = ylength / (double) (jmax);

	/* allocate memory for the arrays */
	U = matrix(0, imax + 1, 0, jmax + 1);
	V = matrix(0, imax + 1, 0, jmax + 1);
	P = matrix(0, imax + 1, 0, jmax + 1);
	F = matrix(0, imax + 1, 0, jmax + 1);
	G = matrix(0, imax + 1, 0, jmax + 1);
	RS = matrix(0, imax + 1, 0, jmax + 1);

	T   = matrix(0, imax + 1, 0, jmax + 1);

	/* allocate memory for the fields */
	Flag = imatrix(0, imax + 1, 0, jmax + 1);
	ChemicalSources = imatrix(0, imax + 1, 0, jmax+1);

	/*************** the algorithm starts here *************************/

	//init_flag1(Problem,problem, imax, jmax, 0, 0, 0, Flag);

	/* init flag field*/
	if (init_flag(Problem, imax, jmax, Flag,ChemicalSources) != 1)
	{
		printf("Invalid obstacle. Program quits ...\n");
		free_imatrix(Problem, 0, imax + 1, 0, jmax + 1);
		free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);
		return -1;
	}

	/* init uvp - here the signature was extended to the problem
	 * string to init u&v in the step case to 0 */
	init_uvp(TI, UI, VI, PI, imax, jmax, problem, Flag, U, V, P, T);

	printf("\n");
	printf("Geometric Domain:\n");
	print_matrix(Flag,0, imax + 1, 0, jmax + 1);
	printf("\n");

	printf("\n");
	printf("Chemical Domain:\n");
	print_matrix(ChemicalSources,0, imax + 1, 0, jmax + 1);
	printf("\n");



//	printf("init U\n");
//	print_matrixD(U,0, imax + 1, 0, jmax + 1);

	t_print = 0;
	while (t < t_end)
	{
		/*calculate the timestep */
        calculate_dt(Re, Pr, tau, &dt, dx, dy, imax,jmax, U,V);

        /*calculate the boundary values   */
		boundaryvalues( imax, jmax,dx,dy, wl, wr, wt, wb, U, V, F, G, P, T, Flag,
				tl,tr, tb,tt);

    	/* set special boundary values*/
	    spec_boundary_val( problem, imax, jmax, dx, dy, Re, deltaP, U, V, P,T);

	    /* calculate new temperature values */
	    calculate_Temp(U, V, T, Flag, imax, jmax, dt, dx, dy, alpha, Re, Pr);


		//printf("calculate_Temp T\n");
		//print_matrixD(T,0, imax + 1, 0, jmax + 1);

		/*calculate F&G* - here the signature was extended to the FLAG
	     * matrix to calculate values only for fluid cells */

	    // some error in fg helper
	    //calculate_fg( Re, GX, GY, alpha, dt, dx,dy, imax, jmax, U, V, F, G, Flag);
		calculate_fg1(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V ,F , G, Flag ,T,beta);

		//printf("G\n");
		//print_matrixD(G,0, imax + 1, 0, jmax + 1);

		/*calculate righthand site - here the signature was extended to the FLAG
	     * matrix to calculate values only for fluid cells */
        calculate_rs(dt,dx, dy, imax, jmax, F, G, RS, Flag);

		/* set initial residual*/
		it = 0;
		res = eps + 1;
		while (it < itermax && res > eps)
		{
			/* here the signature was extended to the FLAG
	         * matrix to calculate values only for fluid cells */
            sor(omg, dx, dy, imax, jmax, P, RS, &res, Flag, problem, deltaP);
			it++;
		}
		/* calculate uv - here the signature was extended to the FLAG
	     * matrix to calculate values only for fluid cells */
	    calculate_uv(dt,dx,dy,imax,jmax, U, V, F, G, P, Flag);


		if (t > t_print)
	 	{
	 		write_vtkFile(output_filename_array, n, imax, jmax, dx, dy, U, V, P, T);
	 		t_print += dt_value;
	 	}

	 	t = t + dt;
		n++;
		if (verbose)
		{
			if (t > t_print)
				printf("write outputfile - step-counter: %i, time: %f, sor-interations: %i  \n",n, t, it);
		}
	}
	//write_vtkFile(output_filename_array, n, imax, jmax, dx, dy, U, V, P);

	/* free arrays */
	free_matrix(U, 0, imax + 1, 0, jmax + 1);
	free_matrix(V, 0, imax + 1, 0, jmax + 1);
	free_matrix(P, 0, imax + 1, 0, jmax + 1);
	free_matrix(F, 0, imax + 1, 0, jmax + 1);
	free_matrix(G, 0, imax + 1, 0, jmax + 1);
	free_matrix(RS, 0, imax + 1, 0, jmax + 1);
	free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);
	return 1;
}
