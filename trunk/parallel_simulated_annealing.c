/*
   This file contains implementation of Parallel Simulated Annealing optimization algorithm with 
   parallel compuation implemented using Massage Passing Interface (MPI).
   
   Please go to my Master's Thesis to understand the algorithm and its implementation: Shridhar Shah, "Musculoskeletal simulation of upper extremity motion: Effect of selective muscle weakness and application to rehabilitation" University of Delaware, Summer 2009.
 
   or visit: "https://sites.google.com/site/shridharshah/projects" for the algorithm flow chart.
 
   This implementation is adapted from Bill Goffe's simulated annealing algorithm implementation
   in fortran for serial computation. Please go to "http://netlib2.cs.utk.edu/opt/simann.f" to learn
   about the original algorithm.

   The present file will not compile as is. This is a portion of actual implementation and a function "COSTFUNCTION" is not implemented here. 
   This file is for reference who wants to implement Parallel Simulated Annealing algorithm with MPI. I am sharing this file with belief that
   someone can benefit from this implementation and adept for their own implementation.
   
   Copyright (C) 2008 - 2009, Shridhar Shah,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Any Feedback is welcome:
   shridhar [at] udel [dot] edu


*******************************************************************************/
#include "universal.h"
#include "model.h" 
#include "math.h"
#include "mpi.h"

/*************** DEFINES (for this file only) *********************************/

/* Tolerances and Integration parameters. */
#define BAUMGARTE_STAB 20.0
#define N 145
#define NEPS 4
/*************** STATIC GLOBAL VARIABLES (for this file only) *****************/
SystemInfo si;

/**************** GLOBAL VARIABLES (used in only a few files) *****************/
dpModelStruct* sdm;
MotionData* kinetics_data = NULL;
dpSimulationParameters params;
char current_dir[CHARBUFFER];

float expt[60][5]; 
double u[97], c, cd, cm; 
int ncalls, i97, j97; e
int gimbal_err = 0;
/*************** EXTERNED VARIABLES (declared in another file) ****************/
extern char buffer[];
extern dpBoolean verbose;

/*************** PROTOTYPES for STATIC FUNCTIONS (for this file only) *********/
static void init_model(void);
/*************** Addition*********/
void sim_drive(int *n, double xx[], int *opt, double *fcn_value);
double track_err(double track[60][6]);

/************************************/
// simulated annealing function definitions

void write_best(double xx[],int n,double best_perf,int ncpu,int nfunc_1, int nfunc_total);
void prtvec(double vector[],int ncols,char name[]);
void prt10();
void prt9(int max,int n,double tt,double xopt[],double vm[],double fopt,int nup,int ndown,int nrej,int lnobds,int nnew);
void prt8(int n,double vm[],double xopt[],double x[]);
void prt7(int max);
void prt6(int max);
void prt5(void);
void prt4(int max,int n,double xp[],double x[],double fp,double f);
void prt3(int max,int n,double xp[],double x[],double f);
void prt2(int max,int n,double x[],double f);
void prt1(void);
double ranmar(void);
double fcn(int n,double x[]);
double exprep(double rdum);
void rmarin(int ij,int kl);
double span(int n, double x[],int max,double rt,double eps,int ns,int nt,int neps,int maxevl,double lb[],double ub[],double c[],int iprint,int iprint2, int iseed1,int iseed2,double tt,double vm[],double xopt[],double fopt,int nacc,int nfcnevtot,int nobds,int ier,double fstar[],double xp[],int nacp[],int work[]);
void finalize(int thiscpu,int n,double xopt[],double vm[],double fopt,int nfcnevtot,int nacc,int nobds,double tt,int ier);
void early_communicate(int ncpu,int thiscpu,double ub[],double lb[],double c[],int nacp[],int total_nacp[],int m,int n,double x[],double xopt[],double fp,double vm[],int nt,int ns,double f,double fopt);
void span_communicate(int ncpu,int thiscpu,double ub[],double lb[],double c[],int nacp[],int total_nacp[],int m,int n,double x[],double xopt[],double fopt,double vm[],int nt,int ns,double f);
//FILE *sa;

int main(int argc, const char * argv[])
{
	

int n, neps;

int i, j, k, l, m;
double lb[N], ub[N], x[N];
double xopt[N], c[N], vm[N];
double fstar[NEPS], xp[N];
double tt, eps, rt, fopt;
int nacp[N], work[N], ns, nt;
int nfcnevtot, ier, iseed1,iseed2;
int maxevl, iprint, iprint2, nacc, nobds;
double fcn_value;
int max;

n = 145;
neps = 4;

FILE *pFile;
FILE *iniFile;

//Initialize MPI..
 
MPI_Init((int*)&argc,(char***) &argv);

/*Run it in simann - main program LOAD EXPERIMENTAL DATA  */

pFile = fopen("expdata.txt","r");
//rewind(pFile);
for (i=0; i<60; i++)
{
fscanf(pFile,"%f %f %f %f %f",&expt[i][0],&expt[i][1],&expt[i][2],&expt[i][3],&expt[i][4]);
//printf("%f %f %f %f %f\n",expt[i][0],expt[i][1],expt[i][2],expt[i][3],expt[i][4]);
}
fclose(pFile);

// Set input parameters
//Recommended values:  nt = 100, ns = even multiples of ncpu
// Recommended tolerance: 1.0E-03


max = FALSE; // 0 for , 1 for 
eps = 50;
rt = 0.75;
iseed1 = 1;
iseed2 = 2;
ns = 24;
nt = 4;
maxevl = 100000;
iprint = 3;
iprint2 = 1;
ncalls = 0;

for (i =  0; i < n; i++ ) 
{
    c[i] = 2.0;
}

//Excitation Magnitudes

for (i =  0; i < n; i++ ) 
{
    lb[i] = 0.0;
	ub[i] = 1.0;
	work[i] = 0;
}

// Upper and Lower bounds for all muscles.
for(i=90+3;i<93+24;i+=2)
{
	lb[i] = 0.0;
	ub[i] = 0.20;
	lb[i+1] = 0.70;
	ub[i+1] = 1.0;
	
}
 for(i=93+24;i<n;i+=2)
{
	lb[i] = 0.0;
	ub[i] = 0.30;
	lb[i+1] = 0.62;
	ub[i+1] = 1.0;
}
 
     
//Initial Guess

iniFile = fopen("initialx.txt","r");
//rewind(iniFile);
for (i=0; i< n; i++)
{
fscanf(iniFile,"%lf ",&x[i]);
//printf("%lf\n",x[i]);
}
fclose(iniFile); 


//Set input values of the input/output parameters.
tt = 10.0;

for(i=0; i<n; i++)
{
	vm[i] = 1.0;
}


nacc = 0;
nobds = 0;
ier = 99;
nfcnevtot = 0;
fopt = 0;



//sim_init();
span(n, x, max, rt, eps, ns, nt, neps, maxevl, lb, ub, c, iprint,iprint2, iseed1, iseed2, tt, vm, xopt, fopt, nacc, nfcnevtot, nobds, ier, fstar, xp, nacp, work);

//fcn(n, x);
return 0;
} // main function ends



double fcn(int n,double x[])
//This subroutine computes the cost function for pedaling simulation
//opt=0 to find cost function, 1 to generate motion file
//value=1e10 for incomplete trial, 2e10 for aborted trial (init error)
{

int opt, i;
double value;

// double x[n], value;
//int n;

opt = 0;
ncalls =  ncalls + 1;
value = 20000000000.0;
COSTFUNCTION(&n, x, &opt, &value); //COSTFUNCTION is not implemented. Once can implement their own cost function.
if(gimbal_err == 1 || isnan(value) != 0)
value = 10000000000;
else
value = value/25.0;
printf("value = %e\n", value);
printf("ncalls = %d\n",ncalls);
return value;
}




double span(int n, double x[],int max,double rt,double eps,int ns,int nt,int neps,int maxevl,double lb[],double ub[],double c[],int iprint,int iprint2, int iseed1,int iseed2,double tt,double vm[],double xopt[],double fopt,int nacc,int nfcnevtot,int nobds,int ier,double fstar[],double xp[],int nacp[],int work[])
{
	//Type all external variables.
	//double x[], lb[], ub[], c[], vm[], fstar[], xopt[], xp[],tt, eps, rt, fopt;
	//int nacp[], n, ns, nt, neps, nacc, maxevl, iprint, nobds, ier, nfcnev, iseed1, iseed2;
	//int max; //logical
//Type all internal variables.
	double f, fp, p, pp, ratio;
	int nup, ndown, nrej, nnew, lnobds, h, i, j, m,l, k, nfcnev;
	int quit; // logical

	// Include all internal CPU variables.
	int ncpu, thiscpu, totalcalcs, ncalcpcpu, ncalcleft;
	int ncalc;
	int rank, ierror, tag;

	ncpu = 0;
	rank = 0;
	FILE * sa;

	MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	thiscpu = rank + 1;

//check for ranmar function:
//	printf("1: %e CPU: %d \n",ranmar(iseed1+rank, iseed2 + 2*rank + 1), rank);


	if(thiscpu == 1)
	{

sa = fopen("sa.txt", "w");
fprintf(sa,"Start of File\n");
fprintf(sa,"Simulated Annealing\n");
fprintf(sa, "Number of Parameters : %d\nMaximization : %d\nInitial Temp : %e\n", n, max, tt);
fprintf(sa, "RT : %e\nEPS : %e\nNS : %d\nNT : %d\nNEPS : %d\nMAXEVL : %d\nIprint : %d\nIseed1 : %d\nIseed2 : %d\n", rt, eps, ns, nt, neps, maxevl, iprint, iseed1, iseed2);
fclose(sa);

prtvec(x, n, "Starting Values");
prtvec(vm, n, "Initial Step Length");
prtvec(lb, n, "Lower Bound");
prtvec(ub, n, "Upper Bound");
prtvec(c, n, "C Vector");


sa = fopen("sa.txt", "a");
fprintf(sa,"\n");
fprintf(sa, "End of driver routine output before call to SA\n");
fclose(sa);

	}

//Determine how many function evaluations each cpu will make
	totalcalcs = ns;
	ncalcpcpu = totalcalcs/ncpu;
	ncalcleft = totalcalcs - (ncalcpcpu * ncpu);
	if(thiscpu <= ncalcleft)
		ncalc = ncalcpcpu + 1;
	else
		ncalc = ncalcpcpu;


    ///Initialize the random number generator ranmar.
	rmarin((iseed1+rank), (iseed2+ 2*rank + 1));

//Set initial values.
	nacc = 0;
	nobds = 0;
	nfcnev = 0;
	nfcnevtot = 0;
	ier = 99;


for(i=0;i<n;i++)
{
	xopt[i] = x[i];
	nacp[i] = 0;
}
for(i=0;i<neps;i++)
{
	fstar[i] = 1.0e+20;
}

//If the initial temperature is not positive, notify the user and return to the calling routine.
if (tt <= 0.0)
{
	sa = fopen("sa.txt", "a");
	fprintf(sa, "The initial temperature is not positive, reset the variable tt\n");
	fclose(sa);
	ier = 3;
	finalize(thiscpu, n, xopt, vm, fopt, nfcnevtot, nacc, nobds, tt, ier);
	return ier;
}

//If the initial value is out of bounds, notify the user and return to the calling routine.
for (i=0;i<n;i++)
{
	if ((x[i] > ub[i]) || (x[i] < lb[i]))
	{
		sa = fopen("sa.txt", "a");
		fprintf(sa, "I: %d\n", i);
		fclose(sa);
		prt1();
		ier = 2;
		finalize(thiscpu, n, xopt, vm, fopt, nfcnevtot, nacc, nobds, tt, ier);
		return ier;
	}
}

// Evaluate the function with input x and return value as f.

gimbal_err = 0;
f = fcn(n, x);

/* If the function is to be minimized, switch the sign of the function.
Note that all intermediate and final output switches the sign back to eliminate any 
possible confusion for the user. */

if (max == FALSE) 
f = -f;

nfcnev = nfcnev+1;
nfcnevtot = nfcnevtot + 1;
fopt = f;
fstar[1] = f;
if ((iprint >= 1) && (thiscpu ==1))
{
	sa = fopen("sa.txt", "a");
	fprintf(sa, "CPU: %d \n", thiscpu);
	fclose(sa);
	prt2(max, n, x, f);
}

//Write best data to file

if((iprint2 >= 1) && (thiscpu == 1))
	{
		write_best(x, n, fopt, ncpu, 1, ncpu);
	}

/* Start the main loop. Note that it terminates if (i) the algorithm
succesfully optimizes the function or (ii) there are too many function 
evaluations (more than maxevl). */
mainloop:
nup =0;
nrej = 0;
nnew =0;
ndown = 0;
lnobds = 0;

for (m = 0;m < nt; m++)
{
	for(j=0;j<ncalc; j++)
	{

		for (h=0;h<n; h++)
		{
			//Generate xp, the trial value of x. Note use of vm to choose xp.
			for(i=0; i<n; i++)
			{
				for(k=0;k<n;k++){
					xp[i] = x[i];
				}

				if(i == h)
					xp[i] = x[i] + (ranmar()*2.0 - 1.0)*vm[i];
				else
					xp[i] = x[i];
			// If XP is out of bounds, select a point in bounds for the trial.
				if ((xp[i] < lb[i]) || (xp[i] > ub[i]))
				{
					if(iprint >=4)
					{
						sa = fopen("sa.txt", "a");
						fprintf(sa, "CPU: %d\n", thiscpu);
						fclose(sa);
						prt3(max, n, xp, x, f);
					}		
					xp[i] = lb[i] + (ub[i] - lb[i])*ranmar();
					lnobds = lnobds + 1;
					nobds = nobds + 1;
					
				}
			}

			//Evaluate the function with the trial point xp and return as fp.
			gimbal_err = 0;
			fp = fcn(n,xp);
			if (max == FALSE) 
				fp = -fp;
			
			nfcnev = nfcnev + 1;
			
			if (iprint >= 3 && fp >= f) 
			{
				sa = fopen("sa.txt", "a");
				fprintf(sa, "CPU: %d\n", thiscpu);
				fclose(sa);
				prt4(max, n, xp, x, fp, f);
			}

			//If too many function evaluations occur, terminate the algorithm.

			if (nfcnev >= maxevl) 
			{
				prt5();
				if (max == FALSE)
					fopt = -fopt;
				
				ier = 1;
				finalize(thiscpu, n, xopt, vm, fopt, nfcnevtot, nacc, nobds, tt, ier);
				return ier;
			}

			//Accept the new point if the function value increases.
			if(fp >= f)
			{
				if (iprint >= 3)
				{
					sa = fopen("sa.txt", "a");
					fprintf(sa, "CPU: %d\n", thiscpu);
					fprintf(sa, "Point Accepted\n");
					fclose(sa);
				}
                
				for(i=0;i<n;i++)
				{
					x[i] = xp[i];
				}
				f = fp;
				nacc = nacc + 1;
				nacp[h] = nacp[h] + 1;
				nup = nup + 1;

				//If greater than any other point, record as new optium.

				if (fp > fopt)
				{
					if ( iprint >= 3)
					{
						sa = fopen("sa.txt", "a");
						fprintf(sa, "CPU: %d\n", thiscpu);	
						fprintf(sa, "New Optimum\n");
						fclose(sa);
					}
				
					for(i=0;i<n;i++)
					{
						xopt[i] = xp[i];
					}
					fopt = fp;
					nnew = nnew + 1;
				}

				//If the point is lower, use the Metropolis criteria to decide on 
				//acceptance or rejection.

				else
				{
					p = exprep((fp-f)/tt);
					pp = ranmar();
					if (pp < p)
					{
						if (iprint >= 3)
						{
							sa = fopen("sa.txt", "a");
							fprintf(sa, "CPU: %d\n", thiscpu);	
							fclose(sa);
							prt6(max);
						}
						for (i=0;i<n;i++)
						{
							x[i] = xp[i];
						}

						f = fp;
						nacc = nacc + 1;
						nacp[h] = nacp[h] + 1;
						ndown = ndown + 1;
					}
					else
					{
						nrej = nrej + 1;
						if (iprint >= 3)
						{
							sa = fopen("sa.txt", "a");
							fprintf(sa, "CPU: %d\n", thiscpu);	
							fclose(sa);
							prt7(max);
						}
					}
				}
			}

			
				
						
				MPI_Barrier(MPI_COMM_WORLD);
			
				early_communicate(ncpu, thiscpu, ub, lb, c, nacp,work, m, n, x, xopt, fp, vm, nt, ns, f,fopt);
			
				//printf("CPU: %d\n",thiscpu);
				//for(i=0;i<n;i++)
				//printf("%e  ",x[i]);
				//printf("\n");
				
			}
	}

	// Communicate among processors to refine the neighborhood.
	MPI_Barrier(MPI_COMM_WORLD);
	span_communicate(ncpu, thiscpu, ub, lb, c, nacp, work, m, n, x, xopt, fopt, vm, nt, ns, f);
		
	// Compute total number of function evaluations
	
	nfcnevtot = nfcnevtot + ns*n;

//Write best data to file

	if ((iprint == 1 )&& (thiscpu == 1))
	{
		// Check - > system("cat best_perf.out >> best_perf.out.all");
		write_best(xopt, n, fopt,ncpu, nfcnev, nfcnevtot);

	}

	if (iprint >= 2)
	{
		if(thiscpu == 1)
		prt8(n,vm,xopt, x);
	}

}

	if (iprint >= 1)
	{
		if(thiscpu == 1)
			prt9(max, n,tt, xopt, vm, fopt, nup, ndown, nrej, lnobds, nnew);
	}
	// Check termination criteria.

	quit = FALSE;
	fstar[1] = f;

	if ((fopt - fstar[1]) <= eps)
		quit = TRUE;

	for ( i=0;i<neps;i++)
	{
		if (ABS(f - fstar[i]) > eps ) 
		{
			quit = FALSE;
			break;
		}
	}

	// Terminate SA if appropriate.

	if (quit == TRUE)
	{
		for (i=0;i<n;i++)
		{
			x[i] = xopt[i];
		}
		ier = 0;
		if (max == FALSE)
			fopt = -fopt;
		if(iprint >= 1) 
		{
			if(thiscpu == 1)
				prt10();
		}
		finalize(thiscpu, n, xopt, vm, fopt, nfcnevtot, nacc, nobds, tt, ier);
		return ier;
	}

	// If termination criteria is not met, prepare for another loop.

	tt = rt*tt;
	for (i=neps-1; i>=1; i--)
	{
		fstar[i] = fstar[i-1];
	}
	f = fopt;
	for (i=0;i<n;i++)
	{
		x[i] = xopt[i];
	}

	// Loop Again
	goto mainloop;



}

void early_communicate(int ncpu,int thiscpu,double ub[],double lb[],double c[],int nacp[],int total_nacp[],int m,int n,double x[],double xopt[],double fp,double vm[],int nt,int ns,double f, double fopt)
{
//local variables
	int size, i, j, k, best_cpu, new_best;
	int tag;
	MPI_Status status;
	double ratio, best_f;
//	printf("determine fp and cpu that had best FOPT\n");

	best_cpu = 0;
	new_best = 0;
	tag = 1;
	
	//printf("Early Communication started, Send data to cpu 1\n");
if(thiscpu > 1)
	{
		MPI_Ssend(&f, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	else
	{
		//Receive data from all other cpu's.
		best_cpu = thiscpu;
		best_f = f;
		for (j=2;j<=ncpu;j++)
		{
			MPI_Recv(&f, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD, &status);	
			if(status.MPI_ERROR != 0)
			{		
				printf("Error receiving FOPT from cpu %d\n", j);
				exit(0);
			}
			//Compare value of FOPT from CPU j to best_fopt
			if(f > best_f)
			{
				best_cpu = j;
				best_f = f;
				new_best = 1;
				printf("best cpu : %d\n",best_cpu);
				printf("best f : %e\n",best_f);
			}
		}
		f = best_f;
	}
MPI_Barrier(MPI_COMM_WORLD);
//printf("Communicate which CPU was best\n");
	tag = 2;
	if(thiscpu == 1)
	{
		for(j=2;j<=ncpu;j++)
		{
			MPI_Ssend(&best_cpu, 1, MPI_INTEGER, j-1, tag, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(&best_cpu, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving best_cpu from master\n");
			exit(0);
		}
	}

printf("Best cpu: %d, CPU: %d\n",best_cpu, thiscpu);

//Update all cpus
	tag = 4;
	//Best cpu sends new fopt, f, vm, x and xopt to other cpu's

	if(thiscpu == best_cpu)
	{
		for(j=1;j<=ncpu;j++)
		{
			if(j != best_cpu)
			{
				MPI_Ssend(&fp, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(&f, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(vm, n, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(x, n, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(xopt, n, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(&fopt, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
			}
		}
	}
	else if(thiscpu != best_cpu)
	{
		MPI_Recv(&fp, 1, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving fopt from best cpu\n");
			exit(0);
		}

		MPI_Recv(&f, 1, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving f from best cpu\n");
			exit(0);
		}

		MPI_Recv(vm, n, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving vm from best cpu\n");
			exit(0);
		}

		MPI_Recv(x, n, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving x from best cpu\n");
			exit(0);
		}

		MPI_Recv(xopt, n, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving xopt from best cpu\n");
			exit(0);
		}
		MPI_Recv(&fopt, 1, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving fopt from best cpu\n");
			exit(0);
		}

	}

printf("f: %e, CPU: %d\n",f, thiscpu);

}


void span_communicate(int ncpu,int thiscpu,double ub[],double lb[],double c[],int nacp[],int total_nacp[],int m,int n,double x[],double xopt[],double fopt,double vm[],int nt,int ns,double f)
{
//local variables
	int size, i, j, k, best_cpu, new_best;
	int tag;
	MPI_Status status;
	double ratio, best_fopt;

//	printf("determine fopt and cpu that had best FOPT\n");

	best_cpu = 0;
	new_best = 0;
	tag = 1;
	
//	printf("Send data to cpu 1\n");
	if(thiscpu > 1)
	{
		MPI_Ssend(&fopt, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	else
	{
		//Receive data from all other cpu's.
		best_cpu = thiscpu;
		best_fopt = fopt;
		for (j=2;j<=ncpu;j++)
		{
			MPI_Recv(&fopt, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD, &status);	
			if(status.MPI_ERROR != 0)
			{		
				printf("Error receiving FOPT from cpu %d\n", j);
				exit(0);
			}
			//Compare value of FOPT from CPU j to best_fopt
			if(fopt > best_fopt)
			{
				best_cpu = j;
				best_fopt = fopt;
				new_best = 1;
			}
		}
		fopt = best_fopt;
	}

	//printf("Communicate which CPU was best\n");
	tag = 2;
	if(thiscpu == 1)
	{
		for(j=2;j<=ncpu;j++)
		{
			MPI_Ssend(&best_cpu, 1, MPI_INTEGER, j-1, tag, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(&best_cpu, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving best_cpu from master\n");
			exit(0);
		}
	}

	//Adjust Step Length
	tag = 3;
	// Send nacp (number of accepted points) to best cpu
	if (thiscpu != best_cpu)
		MPI_Ssend(nacp, n, MPI_INTEGER, best_cpu-1, tag, MPI_COMM_WORLD);
	
	//Receive nacp from cpu's other than the best cpu
	else
	{
		for(j=0;j<n;j++)
		{
			total_nacp[j] = nacp[j];
		}

		for(j = 1;j<=ncpu; j++)
		{
			if(j != best_cpu)
			{
				MPI_Recv(nacp, n, MPI_INTEGER, j-1, tag, MPI_COMM_WORLD, &status);
				if(status.MPI_ERROR != 0)
				{
					printf("Error receiving nacp from cpu %d\n");
					exit(0);
				}
            
				for(k=0;k<n;k++)
				{
					total_nacp[k] = total_nacp[k] + nacp[k];
				}
			}
		}
	//Update NACP
		for(j=0;j<n;j++)
		{
			nacp[j] = total_nacp[j];
		}

	//Adjust vm so that approximately half of all evaluations are accepted.
	
		for(i = 0; i<n; i++)
		{
			ratio = (float)nacp[i]/(float)ns;
			if (ratio > 0.6)
				vm[i] = vm[i]*(1 + c[i]*(ratio - 0.6)/0.4);
			else if (ratio < 0.4)
				vm[i] = vm[i]/(1 + c[i]*(0.4 - ratio)/0.4);
		
			if (vm[i] > (ub[i] - lb[i]))
				vm[i] = ub[i] - lb[i];
		}
	}
MPI_Barrier(MPI_COMM_WORLD);
	//Update all cpus
	tag = 4;
	//Best cpu sends new fopt, f, vm, x and xopt to other cpu's

	if(thiscpu == best_cpu)
	{
		for(j=1;j<=ncpu;j++)
		{
			if(j != best_cpu)
			{
				MPI_Ssend(&fopt, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(&f, 1, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(vm, n, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(x, n, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
				MPI_Ssend(xopt, n, MPI_DOUBLE, j-1, tag, MPI_COMM_WORLD);
			}
		}
	}
	else if(thiscpu != best_cpu)
	{
		MPI_Recv(&fopt, 1, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving fopt from best cpu\n");
			exit(0);
		}

		MPI_Recv(&f, 1, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving f from best cpu\n");
			exit(0);
		}

		MPI_Recv(vm, n, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving vm from best cpu\n");
			exit(0);
		}

		MPI_Recv(x, n, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving x from best cpu\n");
			exit(0);
		}

		MPI_Recv(xopt, n, MPI_DOUBLE, best_cpu-1, tag, MPI_COMM_WORLD,&status);
		if(status.MPI_ERROR != 0)
		{
			printf("Error receiving xopt from best cpu\n");
			exit(0);
		}
	}

	//Reset NACP
	for(i=0;i<n;i++)
	{
		nacp[i] = 0;
	}

	
	
	//return;
}


double exprep(double rdum)
{
	/* This function replaces exp to avoid under - and overflows and is
	designed for IBM 370 type machines. It may be necessary to modify 
	it for other machines. Note that the maximum and minimum values of 
	exprep are such that they has no effect on the algorithm. */

	double exprep;

	if (rdum > 174.0)
		exprep = 3.69e+75;
	else if (rdum < -180.0)
		exprep = 0.0;
	else
		exprep = exp(rdum);
	return exprep;
}

void rmarin(int ij,int kl)
{
	/* This subroutine and the next function generate random numbers. See
	the comments for SA for more information. The only changes from the 
	original code is that (1) the test to make sure that rmarin runs first
	was taken out since SA assures that this is done (this test didn't
	compile under IBM's VS Fortran) and  (2) typing ivec as integer was 
	taken out since ivec isn't used. With these exceptions, all following
	lines are original.

	This is the initialization routine for the random number generator
	ranmar() 

	Note: The seed variables can have balues between 0<= ij <= 31328
	0 <= kl <= 30081 */

	//double u[97], c, cd, cm;
	//int i97, j97;
	int i, j, k, l, ii, jj;
	double s, tt;
	int m;
	//double uni, ranmar;
	if ( (ij < 0) || (ij > 31328) || (kl < 0) || (kl > 30081) )
	{
		printf("The first random number seed must have a value betweeen 0 and 31328\n");
		printf("The second random number seed must have a value betweeen 0 and 30081\n");
		exit(0);
	}

	i = ((ij/177)%177) + 2;
	j = (ij%177) + 2;
	k = ((kl/169)%178) + 1;
	l = (kl%169);

	for(ii = 0; ii<97; ii++)
	{
		s = 0.0;
		tt = 0.5;
		for (jj = 0;jj< 24; jj++)
		{
			m = ((((i*j)%179)*k)%179);
			i = j;
			j = k;
			k = m;
			l = ((53*l+1)%169);
			if (((l*m)%64) >= 32)
				s = s + tt;
			tt = 0.5*tt;
		}
		u[ii] = s;
	}
	c = 362436.0 / 16777216.0;
    cd = 7654321.0 / 16777216.0;
    cm = 16777213.0 /16777216.0;
    i97 = 97;
    j97 = 33;
//    return;
}



double ranmar()
{
	double uni;
	double ranmar;
	
	uni = u[i97-1] - u[j97-1];
	if (uni <= 0.0)
		uni = uni + 1.0;
	u[i97-1] = uni;
	i97 = i97 - 1;
	if (i97 == 0)
		i97 = 97;
	j97 = j97 - 1;
	if(j97 == 0)
		j97 = 97;
	c = c - cd;
	if (c < 0.0)
		c = c + cm;
	uni = uni - c;
	if (uni < 0.0 ) 
		uni = uni + 1.0;
	ranmar = uni;
	return(ranmar);
}

void prt1()
{
/* This subroutine prints intermediate output, as does prt2 through
  PRT10. Note that if SA is minimizing the function, the sign of the
  function value and the directions (up/down) are reversed in all
  output to correspond with the actual function optimization. This
  correction is because SA was written to maximize functions and
  it minimizes by maximizing the negative a function. */
	FILE *sa;
	sa = fopen("sa.txt", "a");
	
    fprintf(sa,"THE STARTING value (x) IS OUTSIDE THE BOUNDS (lb AND ub).\n");
	fprintf(sa,"EXECUTION TERMINATED WITHOUT ANY OPTIMIZATION. RESPECIFY x, ub OR lb SO THAT lb[i] .LT. x[i] .LT. ub[i], i = 1, n.\n");
	fclose(sa);


     
}

void prt2(int max,int n,double x[],double f)
{
	FILE *sa;
	sa = fopen("sa.txt", "a");
	
      prtvec(x,n,"INITIAL x");
      if (max == TRUE)
		  fprintf(sa, "Initial f: %1.18e\n",f);
	  else
         fprintf(sa, "Initial f: %1.18e\n",f);
	fclose(sa);
	
}


void prt3(int max,int n,double xp[],double x[],double f)
{
FILE *sa;
sa = fopen("sa.txt", "a");

      prtvec(x,n,"CURRENT x");
      if (max == TRUE)
         fprintf(sa, "Current f: %1.18e\n",f);
	  else
         fprintf(sa, "Current f: %1.18e\n",-f);

      prtvec(xp,n,"TRIAL x");
	
      fprintf(sa,"POINT REJECTED SINCE OUT OF BOUNDS\n");
fclose(sa);
     
}

void prt4(int max,int n,double xp[],double x[],double fp,double f)
{
	FILE *sa;
sa = fopen("sa.txt", "a");

     // prtvec(x,n,"CURRENT x");
      if (max == TRUE)
	  {
		  fprintf(sa, "Current f: %1.18e\n",f);
		  prtvec(xp,n,"TRIAL x");
		  fprintf(sa,"RESULTING f: %1.14e\n",fp);
	  }
	  else
	  {
         fprintf(sa, "Current f: %1.18e\n",-f);
		  prtvec(xp,n,"TRIAL x");
		  fprintf(sa,"RESULTING f: %1.14e\n",-fp);
	  }
fclose(sa);
      
}
 
void prt5()
{
FILE *sa;
sa = fopen("sa.txt", "a");

      fprintf(sa,"TOO MANY FUNCTION EVALUATIONS; CONSIDER INCREASING maxevl OR eps, OR DECREASING\n");
      fprintf(sa,"nt OR rt. THESE RESULTS ARE LIKELY TO BE POOR.\n");
fclose(sa);

}

void prt6(int max)
{
FILE *sa;
sa = fopen("sa.txt", "a");


      if (max == TRUE)
         fprintf(sa,"THOUGH LOWER, POINT ACCEPTED\n");
	  else
         fprintf(sa,"THOUGH HIGHER, POINT ACCEPTED\n");
fclose(sa);
}


void prt7(int max)
{
   FILE *sa;
sa = fopen("sa.txt", "a");

      if (max == TRUE)
         fprintf(sa,"LOWER POINT REJECTED\n");
	  else
         fprintf(sa,"HIGHER POINT REJECTED\n");
	  fclose(sa);

}

void prt8(int n,double vm[],double xopt[],double x[])
{
FILE *sa;
sa = fopen("sa.txt", "a");
  
      fprintf(sa,"INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT\n");
      prtvec(vm,n,"NEW STEP LENGTH (vm)");
      prtvec(xopt,n,"CURRENT OPTIMAL x");
      prtvec(x,n,"CURRENT x");
      fprintf(sa,"\n");
  fclose(sa);
	
}

void prt9(int max,int n,double tt,double xopt[],double vm[],double fopt,int nup,int ndown,int nrej,int lnobds,int nnew)
{

      int totmov;
	  FILE *sa;
	sa = fopen("sa.txt", "a");
 

      totmov = nup + ndown + nrej;
      
      fprintf(sa,"INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION\n");
      fprintf(sa,"CURRENT TEMPERATURE:        %1.14e\n",tt);

      if (max == TRUE)
	  {
         fprintf(sa,"max FUNCTION value SO FAR:  %1.14e\n",fopt);
         fprintf(sa,"TOTAL MOVES:                %d\n",totmov);
         fprintf(sa,"UPHILL:                  %d\n",nup);
         fprintf(sa,"ACCEPTED DOWNHILL:       %d\n",ndown);
         fprintf(sa,"REJECTED DOWNHILL:       %d\n",nrej);
         fprintf(sa,"OUT OF BOUNDS TRIALS:       %d\n",lnobds);
         fprintf(sa,"NEW MAXIMA THIS TEMPERATURE:%d\n",nnew);
	  }
	  else
	  {
	     fprintf(sa,"MIN FUNCTION value SO FAR:  %1.14e\n",fopt);
         fprintf(sa,"TOTAL MOVES:                %d\n",totmov);
         fprintf(sa,"UPHILL:                  %d\n",nup);
         fprintf(sa,"ACCEPTED DOWNHILL:       %d\n",ndown);
         fprintf(sa,"REJECTED DOWNHILL:       %d\n",nrej);
         fprintf(sa,"OUT OF BOUNDS TRIALS:       %d\n",lnobds);
         fprintf(sa,"NEW MINIMA THIS TEMPERATURE:%d\n",nnew);
	  }

	
      prtvec(xopt,n,"CURRENT OPTIMAL x");
	  prtvec(vm,n,"STEP LENGTH (vm)");
      fprintf(sa,"\n");

	    fclose(sa);

}

void prt10()
{
	FILE *sa;
	sa = fopen("sa.txt", "a");
    fprintf(sa,"SA ACHIEVED TERMINATION CRITERIA. ier = 0.\n");
  fclose(sa);
}


void finalize(int thiscpu,int n,double xopt[],double vm[],double fopt,int nfcnevtot,int nacc,int nobds,double tt,int ier)
{

	

	if(thiscpu == 1)
	{
		FILE *sa;
		sa = fopen("sa.txt", "a");
		fprintf(sa, " Results after SA \n");
		prtvec(xopt, n, "Solution");
		prtvec(vm, n, "Final Step");
		fprintf(sa, "%e %d %d %d %e %d\n", fopt, nfcnevtot, nacc, nobds, tt, ier);
		fclose(sa);
	}

	MPI_Finalize();

	

}

void prtvec(double vector[],int ncols,char name[])
{

/* This subroutine prints the double precision vector named vector.
  Elements 1 thru ncols will be printed. name is a character variable
  that describes vector. Note that if name is given in the call to
  PRTVEC, it must be enclosed in quotes. If there are more than 10
  elements in vector, 10 elements will be printed on each line.*/

	  int lines, i;

	  FILE *sa;
		sa = fopen("sa.txt", "a");
		
      fprintf(sa,"\n%s\n", name);

     // if (ncols > 10) 
		 // lines = (int)(ncols/10.0);

      for (i = 0; i <ncols ; i++)
	  {
            fprintf(sa,"%e ",vector[i]);
			if ((i%10)==0)
				fprintf(sa,"\n");
	  }
       
fprintf(sa,"\n");
fclose(sa);

}

void prtvec2(double vector[],int ncols,char name[])
{

/* This subroutine prints the double precision vector named vector.
  Elements 1 thru ncols will be printed. name is a character variable
  that describes vector. Note that if name is given in the call to
  PRTVEC, it must be enclosed in quotes. If there are more than 10
  elements in vector, 10 elements will be printed on each line.*/

	  int lines, i;

	  FILE *sa;
		sa = fopen("sa.txt", "a");
		
      fprintf(sa,"\n%s\n", name);

     // if (ncols > 10) 
		 // lines = (int)(ncols/10.0);

      for (i = 0; i <ncols ; i++)
	  {
            fprintf(sa,"%e ",vector[i]);
			if ((i%10)==0)
				fprintf(sa,"\n");
	  }
       
fprintf(sa,"\n");
fclose(sa);

}


void write_best(double x[],int n,double best_perf,int ncpu,int nfunc_1, int nfunc_total)
{

/*     This subroutine writes the current best controls and performance,
     number of function evaluations computed on processor 1
     and the total number of function evaluations to best_perf.out.*/

    int i;
	
	FILE *best;
	best=fopen("best_perf.out","w");


    for(i=0;i<n;i++)
	{
        fprintf(best,"x[%d] = %1.14e\n",i,x[i]);
	}
        fprintf(best,"Best performance = %e\n",best_perf);
        fprintf(best,"ncpu = %d\n",ncpu);
		fprintf(best,"nfunc cpu 1 = %d\n",nfunc_1);
		fprintf(best,"nfunc total = %d\n",nfunc_total);
		fprintf(best,"\n");
       

       fclose(best);


}  


/*
SpanInit(int ierror)
{
	MPI_Init(&ierror);
	sigsetup();
	return;
}
*/

// Simulated annealing over.

// ***********************************************************************




double track_err(double track[60][6])
{
	int i;
	double eleangerr, shdeleerr, shdroterr, elbflexerr, prosuperr, jntrecforce;
	double value;
	float eleang, shdele, shdrot, elbflex, prosup, jntrec;
	FILE *cFile;

	
cFile = fopen("costfun.txt","r");
//rewind(pFile);
fscanf(cFile,"%f %f %f %f %f %f",&eleang,&shdele,&shdrot,&elbflex,&prosup, &jntrec);
fclose(cFile);


	eleangerr = 0.0;
	shdeleerr = 0.0;
	shdroterr = 0.0;
	elbflexerr = 0.0;
	prosuperr = 0.0;
	jntrecforce = 0.0;
	value = 0.0;

	for (i=0; i<60; i++){
		//if(ABS(expt[i][0]) > 0.01)
			//eleangerr += pow(ABS(expt[i][0] - track[i][0]), (double)2)/(double)pow(ABS(expt[i][0]),2);
		//else
			eleangerr += pow(ABS(expt[i][0] - track[i][0]), (double)2);//(double)pow(ABS(expt[i][0]),2);
		
		//if(ABS(expt[i][1]) > 0.01)
		//	shdeleerr += pow(ABS(expt[i][1] - track[i][1]),(double)2)/(double)pow(ABS(expt[i][1]),2);
		//else
			shdeleerr += pow(ABS(expt[i][1] - track[i][1]),(double)2);//(double)pow(ABS(expt[i][1])+0.01,2);

		//if(ABS(expt[i][2]) > 0.01)
        //    shdroterr += pow(ABS(expt[i][2] - track[i][2]),(double)2)/(double)pow(ABS(expt[i][2]),2);
		//else
			shdroterr += pow(ABS(expt[i][2] - track[i][2]),(double)2);//(double)pow(ABS(expt[i][2])+0.01,2);
		//if(ABS(expt[i][3]) > 0.01)
        //    elbflexerr += pow(ABS(expt[i][3] - track[i][3]),(double)2)/(double)pow(ABS(expt[i][3]),2);
		//else
			elbflexerr += pow(ABS(expt[i][3] - track[i][3]),(double)2);//(double)pow(ABS(expt[i][3])+0.01,2);
		
		//if(ABS(expt[i][4]) > 0.01)
        //    prosuperr += pow(ABS(expt[i][4] - track[i][4]),(double)2)/(double)pow(ABS(expt[i][4]),2);
		//else
			prosuperr += pow(ABS(expt[i][4] - track[i][4]),(double)2);//(double)pow(ABS(expt[i][4]),2);
	

			jntrecforce += ABS(track[i][5]);
		}
		
   // PRINT RESULTS TO SCREEN  
/*
	printf("\n eleangerr = %e ",eleangerr);
	printf("\n shdeleerr = %e ",shdeleerr);
	printf("\n shdroterr = %e ",shdroterr);
	printf("\n elbflexerr = %e ",elbflexerr);
	printf("\n prosuperr = %e ",prosuperr);
	printf("\n jntrecforce = %e ",jntrecforce);
*/  
 // Compute value of cost function  
  
   value = eleang*eleangerr + shdele*shdeleerr + shdrot*shdroterr + elbflex*elbflexerr + prosup*prosuperr + jntrec*jntrecforce;



   return value;


} 


