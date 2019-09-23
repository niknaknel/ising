/* ----------------------------------------------------------- */
/*      Introductory C program  2016 HCE
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NPOINTS          2e1      // number of points to generate; note 2e2 = 200
#define CONST_SEED   1221439L   // choose some constant seed for random number generators

#include "double_ran1.c"

int main() 
{
  long i,ir;              // long integers
  double x;

  ir = (long) CONST_SEED;
  printf("USING CONSTANT SEED  %ld\n",ir);

  // initialise random number generator:
  ir=-ir; 
  x = double_ran1(&ir); 

  for (i=0; i < NPOINTS; i++) {
    x = double_ran1(&ir);
    printf("%8ld\t%12.8f\n",i,x);
  }
  
  return 0;
}

