/*
   Test the speed of integer and floating point operations
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define TMAX    (long) 1e9
int main() 
{
  clock_t time_start, time_end;
  long x=0L;
  long itt;
  double y=0.0;

  time_start = clock();
  for (itt=0; itt < TMAX; itt++) { x++; }
  printf("integer total = %ld\n",x);
  time_end = clock();
  printf("%g INTEGER OPERATIONS COMPLETED in %g seconds\n\n"
	 , (float) TMAX
	 ,((float) (time_end - time_start) / (float) CLOCKS_PER_SEC));


  time_start = clock();
  for (itt=0; itt < TMAX; itt++) { 
       y += (float) 1.0; 
  }
  printf("float/double total = %g\n",y);
  time_end = clock();
  printf("%g FLOAT OPERATIONS COMPLETED in %g seconds\n\n"
	 , (float) TMAX
	 ,((float) (time_end - time_start) / (float) CLOCKS_PER_SEC));

  return 0;
}
