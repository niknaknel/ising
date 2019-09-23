/* ----------------------------------------------------------- */
/*      Introductory C program  2016 HCE
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define NPOINTS          2e7    // number of points to generate; eg 2e2 = 200
#define CONST_SEED   1221439L   // choose some constant seed for random number generators

#include "double_ran1.c"

int main() 
{
  clock_t time_start, time_end;     // special type for time calculations
  long i,ir, timediff;              // long integers
  double x;

  ir = (long) CONST_SEED;
  printf("USING CONSTANT SEED  %ld\n",ir);

  // initialise random number generator:
  ir=-ir; 
  x = double_ran1(&ir); 

  time_start = clock();
  
  if ((int) NPOINTS < 2e8) {
    for (i=0; i < NPOINTS; i++) {
      x = double_ran1(&ir);
      if(i % (int) 1e6 == 0) {
	printf("%9ld\t%12.8f\n",i,x);
      }
    }
  }
  
  time_end = clock();
  timediff = (long) time_end - (long) time_start;
  printf("\nJOB COMPLETED in %9.2f seconds\n\n", (timediff / (double) CLOCKS_PER_SEC));
  return 0;
}

