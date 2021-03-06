#include "helper.h"
#include "init.h"
#include "NSDefinitions.h"

int read_parameters(
					const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                                        /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
                    int    *wl,
                    int    *wr,
                    int    *wt,
                    int    *wb,
        		    double *dt_value,            /* time for output */
                    double *deltaP,
                    double *TI,
                    double *beta,
                    double *gamma,
                    double *Pr,
                    int *a,
                    int *b,
                    int *c,
                    int *d,
                    double *dH,
                    double *k1,
                    double *k2,
                    double *lambda,
                    int *static_substances,
                    double    *cl,
                    double    *cr,
                    double    *ct,
                    double    *cb
)
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );
   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );
   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );
   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *dt_value );
   READ_DOUBLE( szFileName, *deltaP );
   READ_DOUBLE( szFileName, *TI );
   READ_DOUBLE( szFileName, *beta );
   READ_DOUBLE( szFileName, *gamma );
   READ_DOUBLE( szFileName, *lambda );
   READ_DOUBLE( szFileName, *Pr );

   READ_DOUBLE( szFileName, *cl );
   READ_DOUBLE( szFileName, *cr );
   READ_DOUBLE( szFileName, *ct );
   READ_DOUBLE( szFileName, *cb );

   READ_INT( szFileName, *wl );
   READ_INT( szFileName, *wr );
   READ_INT( szFileName, *wt );
   READ_INT( szFileName, *wb );
   READ_INT   ( szFileName, *itermax );

   /*read source count*/
   READ_INT( szFileName, *static_substances );

   /* read chemical reaction coefficients and enthalpy of reaction */
   READ_INT( szFileName, *a );
   READ_INT( szFileName, *b );
   READ_INT( szFileName, *c );
   READ_INT( szFileName, *d );
   READ_DOUBLE( szFileName, *dH );
   READ_DOUBLE( szFileName, *k1 );
   READ_DOUBLE( szFileName, *k2 );
   return 1;
}

void read_numberOfSources(const char *fileName, int *numberOfSources)
{
	FILE *fp;
    char line[1024];

    if ((fp=fopen(fileName,"r"))==0)
    {
       char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", fileName );
	   ERROR( szBuff );
    }
    printf("Read sources from %s:\n", fileName);
    fgets(line, 1024, fp);
    sscanf(line, "%d", numberOfSources);
    printf("Number of sources: %d\n", *numberOfSources);
    fclose(fp);
}

int read_sources(const char *fileName, double **sources)
{
	int i, numberOfSources;
	double foo, bar;
	/* try to read file and find number of sources */
	FILE *fp;
    char line[1024];

    if ((fp=fopen(fileName,"r"))==0)
    {
       char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", fileName );
	   ERROR( szBuff );
    }
    printf("Read sources from %s:\n", fileName);
    fgets(line, 1024, fp);
    /* dummy line, we need the ones from the second */
    sscanf(line, "%d", &numberOfSources);
	for(i=1; i<=numberOfSources; i++)
	{
		printf("i: %d\n", i);
		if(!feof(fp))
		{
			fgets(line, 1024, fp);
			sscanf(line, "%lf %lf\n", &foo, &bar);
			sources[i-1][0] = foo;
			sources[i-1][1] = bar;
			printf("Source %d: ", i);
			if(!sources[i-1][0])
				printf("start concentration: %f\n", sources[i-1][1]);
			else
				printf("source output %f/sec\n", sources[i-1][1]);
		}
		else
		{
			ERROR("Number of Sources is not equal to the number of lines in the .sources-file!");
		}
	}
	fclose(fp);
	return 0;
}

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_uvp(double TI, double UI, double VI, double PI, int imax, int jmax,
		char* problem, int **Flag, double **U, double **V, double **P, double **T,
		  double*** C, double ***Q, int s_max)
{
	int i;
	int j;
	int k;

	for(j=1; j<=jmax; j++)
	for(i=1; i<=imax; i++)
	{
		if((Flag[i][j]&B_C)==B_C)
		{
			U[i][j] = UI;
			V[i][j] = VI;
			P[i][j] = PI;
			T[i][j] = TI;
		}
		else
		{
			U[i][j] = 0.0 ;
			V[i][j] = 0.0 ;
			P[i][j] = 0.0 ;
		}
	}
	/* set step initial values -
	 * half of the domain must set to zero  */
    if (strcmp(problem, "step") == 0)
    {
        for (i = 0; i <= imax + 1; ++i)
        for (j = 0; j < (jmax+2)/2; ++j)
        	U[i][j] = 0;
    }

    for (k = 0; k < s_max; k++)
    {
        init_matrix(C[k], 1, imax, 1, jmax, 0);
        init_matrix(Q[k], 1, imax, 1, jmax, 0);
    }
}

void init_staticConcentrations(double*** C,int **Sources,
		double **sourceTypeArray, int s_max, int imax,int jmax)
{
    int i = 0;
    int j = 0;
    int k = 0;

    for (k = 0; k < s_max; ++k)
    for (i = 1; i <= imax; ++i)
    for (j = 1; j <= jmax; ++j)
    {
    	if( Sources[i][j] == k+1)
    	{
    		if(sourceTypeArray[k][0] == 1.0)
    		    C[k][i][j] = sourceTypeArray[k][1];
        }
    }
}

void adjust_Concentration(double*** C,int **Sources,double **sourceTypeArray, int s_max, int imax,int jmax){

    int i = 0;
    int j = 0;
    int k = 0;

    for (k = 0; k < s_max; ++k)
    for (i = 1; i <= imax; ++i)
    for (j = 1; j <= jmax; ++j)
    {
    	if( Sources[i][j] == k+1)
    	{
    		double lala = sourceTypeArray[k][0];
    		if(sourceTypeArray[k][0] == 0.0)
    		    C[k][i][j] += sourceTypeArray[k][1];
        }
    }
}

/**  Init the flag field **/
int init_flag(int **Problem,int imax,int jmax, int **Flag, int **Sources)
{
    int i, j;

    /* first read chemical source domain */
    for (i = 1; i < imax+1; i++) {
        for (j = 1; j < jmax+1; j++) {
            if(Problem[i][j] != 0 && Problem[i][j] != 255 && Problem[i][j] != 128)
            {
                Sources[i][j] = Problem[i][j];
            }
            else
            {
                Sources[i][j] = 0;
            }
	    }
    }

    /* normalize problem - due to read_pgm returns not always 1 and zero,
     * the problem needs to be normalized */
    for (i = 1; i < imax+1; i++)
    for (j = 1; j < jmax+1; j++)
    {
    	if(Problem[i][j] != 0 )
    		Problem[i][j] = 1;
    }

    for(i = 1; i < imax+1; i++)
    for(j = 1; j < jmax+1; j++)
    {
        if(Problem[i][j] == 1)
        {
        	/* set all fluid cells to C_F */
        	Flag[i][j] = C_F;
        }
        else
        {
        	/* calculate flag here: 8*east + 4*west + 2*south + 1*north */
         	Flag[i][j] = 8 * Problem[i+1][j] +  4 * Problem[i-1][j]  + 2 * Problem[i][j-1] + 1 * Problem[i][j+1];

        	/* test if flag is valid */
            if(Flag[i][j] == 3 || Flag[i][j] == 7 || Flag[i][j] == 11 || Flag[i][j] == 12 || Flag[i][j] == 13 || Flag[i][j] == 14 || Flag[i][j] == 15)
            	return -1;
        }
    }

    /*boundary cells*/
    for(i = 1; i < imax+1; i++)
    {
        Flag[i][0] = 1 * Problem[i][1];
        Flag[i][jmax+1] = 2 * Problem[i][jmax];
    }

    for(j = 1; j < jmax+1; j++)
    {
        Flag[0][j] = 8 * Problem[1][j];
        Flag[imax+1][j] = 4 * Problem[imax][j];
    }
    return 1;
}
