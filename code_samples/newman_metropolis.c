/*********************************************************************/
/*
  Metropolis algorithm as published in Newman and Barkema
  "Monte Carlo Methods in Statistical Physics"
  page 434.

  The algorithm uses helical boundary conditions.
 */
/*********************************************************************/

#include <math.h>

#define N  (L*L)
#define XNN 1
#define YNN L

/* GLOBAL DECLARATIONS: */
int s[N];   // The Lattice
double prob[5];  // Flip probability for given temperature
double beta;     // Inverse temperature 1/kT

/*********************************************************************/
void initialize()
{
  int i;
  for (i=2; i < 5; i += 2) prob[i] = exp(-2*beta*i);
}

/*********************************************************************/
void sweep()
{
  int i, k;
  int nn,sum,delta;
  
  for (k = 0; k < N; k++) {

    /* Choose a site */
    /* NB You should replace drandom here and below with a random
       number generator that has been tested and performs
       satisfactorily */
    i = N*drandom();
    
    /* Calculate the sum of the neighbouring spins */
    if ((nn = i+XNN) >= N)   nn -= N;
    sum = s[nn];
    if ((nn = i-XNN) <  0)   nn += N;
    sum += s[nn];
    if ((nn = i+YNN) >= N)   nn -= N;
    sum += s[nn];
    if ((nn = i-YNN) <  0)   nn += N;
    sum += s[nn];
    
    /* Calculate the change in energy */
    delta = sum*s[i];
    /* Decide whether to flip the spin */
    if (delta <= 0) {
      s[i] = -s[i];
    } else if (drandom() <  prob[delta]) {
      s[i] = -s[i];
    }
  }
}

/* END OF NEWMAN METROPOLIS ALGORITHM */
/*********************************************************************/
