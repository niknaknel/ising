/* ----------------------------------------------------------- */
/*      Introductory C program  2016 HCE
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define NPOINTS          2e5    // number of points to generate; eg 2e2 = 200
#define CHOOSE_TIMESEED    1    // = 1 for a time seed, 0 for constant seed
#define CONST_SEED   1221439L   // choose some constant seed for random number generators
#define WHICH_RNG  1            // choose a particular random number generator (see list below)

#include "double_ran0.c"
#include "double_ran1.c"
#include "double_ran2.c"
#include "double_ran3.c"
#include "double_ranbinder.c"

int main() 
{
  clock_t time_start, time_end;     // special type for time calculations
  long i,ir, timediff;              // long integers
  double x;
  int rng;
  double *array = (double *) malloc((long) NPOINTS*sizeof(double));
  char *name_rng[] = {"double_ran0", "double_ran1", "double_ran2", "double_ran3"
		      ,"empty","empty","empty","empty","empty","empty"
		      ,"double_ranbinder","double_randu"};
  if(array == NULL) {
    printf("ERROR: requested size of array is too large");
    exit(1);
  }

  FILE* dataf;
  char *outname = "some_name.output";  // data file, for off-line use
  dataf = fopen(outname,"w+");

  if (CHOOSE_TIMESEED) {
    // time = builtin function which give time as an integer; changes once a second
    ir = (long) time(0);  
    printf("USING TIME SEED  %ld\n",ir);
  } 
  else {
    ir = (long) CONST_SEED;
    printf("USING CONSTANT SEED  %ld\n",ir);
  }

  rng = (int) WHICH_RNG;
  printf("WHICH_RNG = %d,\t%s\n\n\ti\t  array[i]\n",rng,name_rng[rng]);
  fprintf(dataf,"WHICH_RNG = %d,\t%s\ni  array[i]\n",rng,name_rng[rng]);

  // Initialise random number generator
  switch (rng) {
  case  0: {break;}                               // ran0 needs no initialisation
  case  1: {ir=-ir; x = double_ran1(&ir); break;} // ran1 needs negative number to initialise
  case  2: {ir=-ir; x = double_ran2(&ir); break;} // ran2 needs negative number to initialise
  case  3: {ir=-ir; x = double_ran3(&ir); break;} // ran3 needs negative number to initialise
  case 10: {break;}                               // ranbinder, no initialisation needed
  default: {printf("\nrng = %d is not a valid RNG\n", rng); 
      exit(1); }
  }

  time_start = clock();

  for (i=0; i < NPOINTS; i++) {
     switch (rng) {
     case 0: {array[i] = double_ran0(&ir); break;}
     case 1: {array[i] = double_ran1(&ir); break;}
     case 2: {array[i] = double_ran2(&ir); break;}
     case 3: {array[i] = double_ran3(&ir); break;}
     case 10: {array[i] = double_ranbinder(&ir); break;}
     }
  }
  
  if ((int) NPOINTS < 2e8) {
    for (i=0; i < NPOINTS; i++) {
      if(i % (int) 1e4 == 0) {
	fprintf(dataf,"%ld %10.8f\n",i,array[i]);
      }
      if(i % (int) 1e6 == 0) {
	printf("%9ld\t%12.8f\n",i,array[i]);
      }
    }
  }
  fclose(dataf);

  time_end = clock();
  timediff = (long) time_end - (long) time_start;
  printf("\nJOB COMPLETED in %9.2f seconds\n\n", (timediff / (double) CLOCKS_PER_SEC));
  printf("Delete the large output file if necessary\n\n");
  free(array);  // release memory used for array
  return 0;
}

