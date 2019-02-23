/* Coded by T.Takahama */

#include <stdio.h>
#include <stdlib.h>

/* Random number generator in [0, 1] */
#define Rand()		((double)rand()/RAND_MAX)

/* Structure of a particle */
typedef struct {
	double *x;
	double *v;
	double f;
	double pbest;
	double *x_star;
} ParticleRec, *Particle;

/*
	Parameters for PSO
*/

/* Number of particles */
#define Nparticles	30
/* Maximum number of iterations */
#define T_MAX		100000
/* The value of inertia weight at t=0 (W_0) and t=T_MAX (W_T) */
#define W_0		0.9
#define W_T		0.4
#define MAX_V		2.0
/* The cognitive parameter (c1) and the social parameter (c2) */
#define c1		2.0
#define c2		2.0

/*
	Definitions for a problem
*/

/* Number of variables: problem dependent */
#define Nvariables	5

/* Objective function for minimization: problem dependent */
#define better(y1, y2)	(y1<y2)

/* The following is the function of Sum_i (x_i-1)^2 */
void Evaluate(Particle P)
{
	int i, j; 
	double temp = 0;

	P->f=0.0;
	for(i=0; i<Nvariables; i++) {
		for(j = 0; j < i; j++) {
			temp += P->x[j];
		}
		P->f += temp * temp;
	}
}

/* update pbest */
void UpdateBest(Particle P)
{
	int j;

	for(j=0; j<Nvariables; j++) P->x_star[j]=P->x[j];
	P->pbest=P->f;
}

/* Initialization of particles: problem dependent */
/* The function returns the index of the best particle */
int Initialize(Particle P, int n)
{
	int i, j;
	int G;		/* the index of the best particle */

	G=0;
	for(i=0; i<n; i++) {
		for(j=0; j<Nvariables; j++) {
			P[i].x[j]= -5.12 + 10.24 * Rand();	/* problem dependent */
			P[i].v[j]=0.0;		/* problem dependent */
		}
		Evaluate(&P[i]);
		UpdateBest(&P[i]);
		if(better(P[i].f, P[G].f)) G=i;
	}
	return G;
}

/*
	Definitions for PSO
*/

/* allocate new data structures */
#define New(type, n, msg)	(type *)NewCell(sizeof(type), n, msg)

void *NewCell(int size, int n, char *msg)
{
	void *new;

	if((new=malloc(size*n))==NULL) {
		fprintf(stderr, "Cannot allocate memory for %d %s\n", n, msg);
		exit(1);
	}
	return new;
}

/* allocate "n" new particles */
Particle NewParticles(int n)
{
	int i;
	Particle P;

	P=New(ParticleRec, n, "particles");
	for(i=0; i<n; i++) {
		P[i].x=New(double, Nvariables, "x");
		P[i].v=New(double, Nvariables, "v");
		P[i].x_star=New(double, Nvariables, "x*");
	}
	return P;
}

/* Print a particle */
void Print(Particle P)
{
	int j;

	for(j=0; j<Nvariables; j++)
		printf("%f ", P->x_star[j]);
	printf(" = %g\n", P->pbest);
}

/* Particle Swarm Optimization */
main()
{
	int t, i, j;
	Particle P;
	int G;
	double w;

	P=NewParticles(Nparticles);
	G=Initialize(P, Nparticles);
	w=W_0;
	for(t=1; t<=T_MAX; t++) {
		for(i=0; i<Nparticles; i++) {
			for(j=0; j<Nvariables; j++) {
				P[i].v[j]=w*P[i].v[j]
						+c1*Rand()*(P[i].x_star[j]-P[i].x[j])
						+c2*Rand()*(P[G].x_star[j]-P[i].x[j]);
				if(P[i].v[j]<-MAX_V)
					P[i].v[j]=-MAX_V;
				else if(P[i].v[j]>MAX_V)
					P[i].v[j]=MAX_V;
				P[i].x[j]+=P[i].v[j];
			}
			Evaluate(&P[i]);
			if(better(P[i].f, P[i].pbest)) {
				if(better(P[i].f, P[G].pbest)) G=i;
				UpdateBest(&P[i]);
			}
		}
		printf("%4d: ", t); Print(&P[G]);
		w-=(W_0-W_T)/T_MAX;
	}
}