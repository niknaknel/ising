/**
 * Playing around with the Ising Model
 *
**/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "double_ran0.c"

void initialize(void);

/* GLOBAL DECLARATIONS: */
double beta;
double prob[5];

int main(int argc, char *argv[])
{
    double T;

    T = atof(argv[1]);
    printf("T = %.2f\n", T);

    // Test initialization
    beta = 1/T;
    initialize();

    // Display prob
    int i;
    for (i = 0; i < 5; i++) {
        printf("%.2f ", prob[i]);
    }
    printf("\n");

}

void initialize(void)
{
  int i;
  for (i=2; i < 5; i += 2) prob[i] = exp(-2*beta*i);
}