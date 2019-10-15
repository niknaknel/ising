/**
 * Ising Model
 * Author: Annika Nel 19907281
**/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "ran0.h"
#include "double_ran0.h"

#define L 5
#define N  (L*L)
#define XNN 1
#define YNN L
#define T_ZERO_NEG 0  // uniform 1s initial state
#define T_ZERO_POS 1  // uniform -1s initial state
#define T_INF 2       // random initial state
#define J 1

<<<<<<< HEAD
/* Functions */
void initialize(int state, int time_seed);
=======
/* Prototypes */
void initialize(int state);
>>>>>>> d95e89720bfd20f6f9baac51b37564694771e274
void sweep();
void display_lattice();
long *gen_seed(int time_seed);

/* Globals */
int s[N];         // The Lattice
double prob[5];   // Flip probability for given temperature
double beta;      // Inverse temperature 1/kT
int E;         // Instantaneous energy
int M;         // Instantaneous magnetisation
long *SEED;       // random seed


/**
 * Main simulation loop.
 * Usage: ./ising temp init_state time_seed
 * where
 *    temp [double]     :  starting temperature
 *    init_state [int]  :  initial lattice state; 0 -> (T=0) all -1s; 1 -> (T=0) all 1s; 2 -> (T=inf) random
 *    time_seed [int]   :  1 if time_seed is to be used, else 0 
 * 
 */
int main(int argc, char *argv[])
{
    double T;
    int init_state, time_seed;

    /* Handle arguments */
    T = atof(argv[1]);
    init_state = atoi(argv[2]);
    time_seed = atoi(argv[3]);

    printf("E,M\n");

    /* Initialize lattice */
    beta = 1/T;
    initialize(init_state, time_seed);

    // display_lattice();

    int i;
    for (i = 0; i < 10; i++) sweep();

    /************ end ************/
    free(SEED);
    exit(1);

}

void initialize(int state, int time_seed)
{
  int i;

  /* Precalculate probabilities */
  for (i = 2; i < 5; i += 2) prob[i] = exp(-2*beta*i);

  /* Generate random seed */
  SEED = gen_seed(time_seed);

  /* Initialize lattice */
  switch(state) 
  {
    case T_ZERO_POS:
      for (i = 0; i < N; i++) s[i] = 1;
      break;

    case T_ZERO_NEG:
      for (i = 0; i < N; i++) s[i] = -1;
      break;

    case T_INF:
      for (i = 0; i < N; i++) {
        s[i] = double_ran0(SEED) < 0.5 ? -1 : 1;
      }
      break;

    default:
      for (i = 0; i < N; i++) s[i] = 1;
      break;
  }
  
  /* Calculate initial energy */
  int sum, nn;
  for (i = 0; i < N; i++) {
    if ((nn = i+XNN) >= N)   nn -= N;
    sum = s[nn];
    if ((nn = i+YNN) >= N)   nn -= N;
    sum += s[nn];
  }
  E = -J*sum;

  /* Calculate initial magnetisation */
  M = 0;
  for (i = 0; i < N; i++) M += s[i];
  
  printf("%d, %d\n", E, M);
}

long *gen_seed(int time_seed)
{
  long* r = malloc(sizeof(long));

  if (time_seed) {
    *r = time(0);
  } else {
    *r = 123123;
  }

  return r;
}

void sweep()
{
  int i, k, nn, sum, delta;
  
  for (k = 0; k < N; k++) {

    /* Choose a site */
    i = (int) (N * ran0(SEED));
    
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
    if (delta <= 0 || double_ran0(SEED) <  prob[delta]) {
      // printf("flip! %d\n", i);
      s[i] = -s[i];

      /* Update energy and magnetisation */
      E += delta;
      M += 2*s[i];
      
    }

    printf("%d, %d\n", E, M);
    // display_lattice();
    
  }
}

void display_lattice()
{
  int i, j;
  for (j = 0; j < L; j++) {
    for (i = 0; i < L; i++) {
      printf("%2d ", s[j*L+i]);
    }
    printf("\n");
  }
  printf("\n");
}
