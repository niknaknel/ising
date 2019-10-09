/**
 * Playing around with the Ising Model
 *
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
#define ONES 0 // uniform 1s initial state
#define NEG_ONES 1 // uniform -1s initial state
#define INF 2  // random initial state

void initialize(int state);
void sweep();
void display_lattice();

/* GLOBAL DECLARATIONS: */
int s[N];         // The Lattice
double prob[5];   // Flip probability for given temperature
double beta;      // Inverse temperature 1/kT

/**
 * Main simulation loop.
 * Usage: ./ising [temperature:float] [init_state:int=0 (ZERO), 1 (INF)]
 * 
**/
int main(int argc, char *argv[])
{
    double T;
    int init_state;

    T = atof(argv[1]);
    init_state = atoi(argv[2]);
    printf("T = %.2f\n", T);

    // Test initialization
    beta = 1/T;
    initialize(init_state);

    // Display lattice
    int i;
    for (i = 0; i < 5; i++) {
        printf("%.5f ", prob[i]);
    }
    printf("\n");
    display_lattice();
    sweep();
    display_lattice();

    exit(1);

}

void initialize(int state)
{
  int i;

  // Precalculate probabilities
  for (i = 2; i < 5; i += 2) prob[i] = exp(-2*beta*i);

  switch(state) 
  {
    case ONES:
      for (i = 0; i < N; i++) s[i] = 1;
      break;

    case NEG_ONES:
      for (i = 0; i < N; i++) s[i] = -1;
      break;

    case INF:
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