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
#define SEED 123123

/* Prototypes */
void initialize(int state);
void sweep();
void display_lattice();
long *gen_seed();

/* Globals */
int TIME_SEED;    // 1 if time seed is to be used else 0
int s[N];         // The Lattice
double prob[5];   // Flip probability for given temperature
double beta;      // Inverse temperature 1/kT

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
    int init_state;

    /* Handle arguments */
    T = atof(argv[1]);
    init_state = atoi(argv[2]);
    TIME_SEED = atoi(argv[3]);

    /* Initialize lattice */
    beta = 1/T;
    initialize(init_state);

    // testing random gen
    int i;
    for (i = 0; i < 10; i++) {
      long *r = gen_seed(); // not working correctly
      printf("%.5f\n", double_ran0(r));
      free(r);
      
    }

    exit(1);

}

void initialize(int state)
{
  int i;

  /* Precalculate probabilities */
  for (i = 2; i < 5; i += 2) prob[i] = exp(-2*beta*i);

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
      s[0] = 1;
      for (i = 1; i < N; i++) s[i] = s[i-1] * -1;
      break;

      // Still need to figure out random initialization
      // for (i = 0; i < N; i++) {
      //   long p = clock();
      //   s[i] = (int) ran0(&p) < 0.5;
      // }
      // break;

    default:
      for (i = 0; i < N; i++) s[i] = 1;
      break;
  }
  
}

long *gen_seed()
{
  long* r = malloc(sizeof(long));

  if (TIME_SEED) {
    *r = time(0);
  } else {
    *r = SEED;
  }

  return r;
}

void sweep()
{
  int i, k, nn, sum, delta;
  long seed;
  
  for (k = 0; k < N; k++) {

    /* Choose a site */
    seed = 1234;
    i = ran0(&seed);
    
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
    seed = time(0);
    printf("%.5f\n", double_ran0(&seed)); // meep not working lekker
    /* Decide whether to flip the spin */
    if (delta <= 0) {
      s[i] = -s[i];
      printf("flip!\n");
    } else if (double_ran0(&seed) <  prob[delta]) {
      s[i] = -s[i];
      printf("flip!\n");
    }
  }
}

void display_lattice()
{
  int i, j;
  for (j = 0; j < L; j++) {
    for (i = 0; i < L; i++) {
      printf("%d ", s[j+i]); // align!
    }
    printf("\n");
  }
  printf("\n");
}
