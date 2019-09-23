/* ----------------------------------------------------------- */
/*      Introductory C program  2014 HCE
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define NPOINTS          2e5      // number of points to generate; note 2e2 = 200
#define CHOOSE_TIMESEED    0    // = 1 for a time seed, 0 for constant seed
#define CONST_SEED   1221439L   // choose some constant seed for random number generators
#define WHICH_RNG  1            // choose a particular random number generator (see list below)

#include "double_ran0.c"
#include "double_ran1.c"
#include "double_ran2.c"
#include "double_ran3.c"
#include "double_ranbinder.c"

int main() 
{
  FILE* dataf;
  char *outname = "some_name.output";  // data file, for off-line use
  clock_t time_start, time_end;        // special type for time calculations
  int rng;
  long i,ir, timediff;              // long integers
  double x;
  double *array = (double *) malloc((long) NPOINTS*sizeof(double));
  char *name_rng[] = {"double_ran0", "double_ran1", "double_ran2", "double_ran3"
		      ,"empty","empty","empty","empty","empty","empty"
		      ,"double_ranbinder","double_randu"};

  if(array == NULL) {
    printf("ERROR: requested size of array is too large");
    exit(1);
  }

  dataf = fopen(outname,"w+");

  if (CHOOSE_TIMESEED) {
    ir = (long) time(0);  // time = builtin function which give time as an integer
    printf("USING TIME SEED  %ld\n",ir);
  } 
  else {
    printf("USING CONSTANT SEED  %ld\n",ir);
    ir = (long) CONST_SEED;
  }

  rng = (int) WHICH_RNG;
  printf("WHICH_RNG = %d,\t%s\n\ti\t  array[i]\n",rng,name_rng[rng]);
  fprintf(dataf,"WHICH_RNG = %d,\t%s\ni  array[i]\n",rng,name_rng[rng]);

  // Initialise random number generator
  switch (rng) {
  case  0: {break;}                                // ran0 needs no initialisation
  case  1: {ir=-ir; x = double_ran1(&ir); break;} // ran1 needs negative number to initialise
  case  2: {ir=-ir; x = double_ran2(&ir); break;}
  case  3: {ir=-ir; x = double_ran3(&ir); break;}
  case 10: {break;}                                // ranbinder, no initialisation needed
  default: {printf("\nrng = %d is not a valid RNG\n", rng); 
      exit(1); }
  }

  time_start = clock();
  
  for (i=0; i < NPOINTS; i++) {
     switch (rng) {
     case 0: {array[i] = double_ran0(&ir);}
     case 1: {array[i] = double_ran1(&ir);}
     case 2: {array[i] = double_ran2(&ir);}
     case 3: {array[i] = double_ran3(&ir);}
     case 10: {array[i] = double_ranbinder(&ir);}
     }
  }
  
  //  if ((int) NPOINTS < 2e3) {
    for (i=0; i < NPOINTS; i++) {
      printf("%8ld\t%12.8f\n",i,array[i]);
      fprintf(dataf,"%ld %10.8f\n",i,array[i]);
    }
    //}
  fclose(dataf);

  time_end = clock();
  timediff = (long) time_end - (long) time_start;
  printf("JOB COMPLETED in %15.2f seconds\n", (timediff / (double) CLOCKS_PER_SEC));
  free(array);  // release memory used for array
  return 0;
}

