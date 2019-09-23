/* ----------------------------------------------------------- */

/*   Perform various tests on various random number generators 

     HCE  2019
*/

/* ----------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "parameters_testing_generators.c"

/* UTILITIES: */
//#include "gnu_plot_2d.c"
//#include "gnu_plot_3d.c"

/*--------------------------------------------------------------------*/

/*    AVAILABLE RANDOM NUMBER GENERATORS: */

/*    rand() - built-in C generator, initialised with  srand */
/*    random() built-in C generator, initialised with  srandom */
/*    drand48,erand48,lrand48,nrand48 */
/*    For details, see basic_randomnumbers.c */
/*    or type man 3 <generator> */
   
/*    The following RNGs are provided external from C:  */
#include "double_ran0.c"
#include "double_ran1.c"
#include "double_ran2.c"
#include "double_ran3.c"
#include "ran4.c"           // double version does not work
#include "psdes.c"          // needed for ran4
#include "double_ranbinder.c"
#include "double_randu_ibm.c"
#include "double_gasdev.c"  // Gaussian variates using Box-Muller method

/* --------------------------------------------------- */
//  RANDOM NUMBER GENERATOR CONSTANTS AND VARIABLES

#define CHOOSE_TIMESEED    0  // = 1 for a time seed, 0 for constant seed
#define CONST_SEED 1221439L   // constant seed for random number generators

// globals:
int rng_initialise = 1;     // flag governing initialisation of RNGs
char *name_rng[] = {"double_ran0", "double_ran1", "double_ran2", "double_ran3", "ran4"
		    ,"empty" ,"empty" ,"empty" ,"empty" ,"empty"
		    ,"double_ranbinder", "double_randu_ibm"};
/* --------------------------------------------------- */
// prototypes:
int map2d(int which_rng);
int map3d(int which_rng);
int sum_rannos(int which_rng);
double genrand(int which_rng, long *ircall);
int gnu_plot_2d(FILE* pd, char *filename, int xcol, int ycol, int linestyle);
int gnu_plot_3d(FILE* pd, char *filename, int xrot, int zrot, int linestyle);

// Function pointers (one day this may work better than the
// existing code, but is not working yet):
//static double (*randp) (long *) = &double_ranbinder; 
//static double (*randp) (long *) = &double_randu_ibm; 

/*--------------------------------------------------------------------*/
int main() 
{
  int which_rng, which_test;
  clock_t time_start, time_end;

  time_start = clock();
  
  if (AVOID_MENU) {
    which_rng  = 10;
    which_test = 2;
  }
  else {
    printf("CHOOSE RANDOM NUMBER GENERATOR:\n  Type\n"
	   "\t-1  to stop\n"
	   "\t 0  for  double_ran0\n"
	   "\t 1  for  double_ran1\n"
	   "\t 2  for  double_ran2\n"
	   "\t 3  for  double_ran3\n"
	   "\t 4  for  ran4\n"
	   "\t10  for  double_ranbinder\t(bad generator!)\n"
	   "\t11  for  double_randu_ibm\t(bad generator!)\n"
	   );
    scanf("%d",&which_rng);
    if (which_rng < 0) {exit(0);};
    printf("CHOOSE THE TEST:\n  Type\n"
	   "\t 0  for  2-D return map\n"
	   "\t 1  for  3-D return map\n"
	   "\t 2  for  histogram of sum of variates\n"
	   );
    scanf("%d",&which_test);
  }
  switch (which_test) {
    case 0: {
      map2d(which_rng); 
      break;
    }
    case 1: {
      map3d(which_rng); 
      break;
    }
    case 2: {
      printf("PLEASE SUPPLY YOUR OWN FUNCTION sum_rannos\n\n");
      break;
    }
    default: {printf("\nwhich_test = %d is not a valid test case\n"
		     , which_test); exit(1); }
  }

  time_end = clock();
  printf("JOB COMPLETED in %g seconds\n"
	, ((float) (time_end - time_start) / (float) CLOCKS_PER_SEC));
  return 0;
}

/*--------------------------------------------------------------------*/
/*
  Draw one-dimensional histogram of sum of random numbers
 */


/*--------------------------------------------------------------------*/
/* Test RNG for unwanted structure by plotting a "return map". xp and
 yp are both random number points in 2 dimensions. A bad RNG will show
 regularities in the plot, a good one will not. While regularities
 prove that a RNG is bad, the absence of regularities does not prove a
 RNG is good.
*/
int map2d(int which_rng) {
  FILE* dataf;        // data file 
  FILE* plotcommandf; // command file, to be fed to gnuplot
  FILE *pipef;        // pipe to gnuplot
  char *outname = "2dmap.output"; // data file, for off-line use
  char *cmdname = "cdf.output";  // gnuplot command file
  // Initialise integer random value to constant seed:
  long isd = CONST_SEED;  
  double xp,yp; // random points, to be plotted against each other
  int i;

  dataf = fopen(outname,"w+");

  if (GNUPLOT_ON) {
    plotcommandf = fopen(cmdname,"w+");
    fprintf(plotcommandf,
	    "set title \"Random number generator: %s\" offset 0,-1\n"
	    "set xlabel \"x_{i}\"   offset 0,0.5\n"
	    "set ylabel \"x_{i+1}\" offset 0.5,0\n"
	    ,name_rng[which_rng]);
    fclose(plotcommandf);
    //initialize pipe, give gnuplot initial commands through cdf.output
    pipef = popen("/usr/bin/gnuplot -persist cdf.output -", "w");
    fflush(pipef);          // force it through the pipe
    // printf("pipe opened to gnuplot\n");
  }
  // Insert information into data file
  if (OUTPUT_PRINT) {
    fprintf(dataf,"# Random number generator = %s\n"
	    "# CHOOSE_TIMESEED = %d\n"
	    "# CONST_SEED = %ld\n"
	    "# Number of data points generated = %d\n#\n"
	    ,name_rng[which_rng],CHOOSE_TIMESEED,CONST_SEED,NPOINTS_MAP);
  };

  // Initialise random number generator:
  xp = genrand(which_rng, &isd);

  // Fill up first line
  xp = genrand(which_rng, &isd);
  yp = genrand(which_rng, &isd);
  if (OUTPUT_PRINT) { fprintf(dataf,"%.4f %.4f\n",xp,yp);};
  //printf("%.4f  %.4f\n",xp,yp);
  // Subsequent lines
  for(i=1; i < NPOINTS_MAP; i++) {
    xp = yp;
    yp = genrand(which_rng, &isd);
    if (OUTPUT_PRINT) { fprintf(dataf,"%.4f %.4f\n",xp,yp);};
    //printf("%.4f  %f.4\n",xp,yp);
  }

  fclose(dataf);
  if (GNUPLOT_ON) {
    gnu_plot_2d(pipef, outname, 1, 2, 2);
    fclose(pipef);
  }
  return 0;
}

/*--------------------------------------------------------------------*/
/* Test RNG for unwanted structure by plotting a "return map". xp yp
 and zp are random number points in 3 dimensions. A bad RNG will
 show regularities in the plot. While regularities prove that a RNG
 is bad, the absence of regularities does not prove a RNG is good.
*/
int map3d(int which_rng) {
  FILE* dataf;        // data file 
  FILE* plotcommandf; // command file, to be fed to gnuplot
  FILE *pipef;        // pipe to gnuplot
  // Initialise integer random value to constant seed:
  long isd = CONST_SEED;
  double xp,yp,zp; // random points, to be plotted against each other
  int i;
  int xrotation, zrotation, linestyle;
  char *outname = "3dmap.output";
  char *cmdname = "cdf.output";  // gnuplot command file

  dataf = fopen(outname,"w+");
  if (GNUPLOT_ON) {
    plotcommandf = fopen(cmdname,"w+");
    //"set title \"Random number generator: %s\" ,-1\n"
    fprintf(plotcommandf,
	    "set title \"Random number generator: %s\" offset 0,-1\n"
	    "set xlabel \"x_{i}\"   offset 0,0.5\n"
	    "set ylabel \"x_{i+1}\" offset 0.5,0\n"
	    "set zlabel \"x_{i+2}\" offset 0.5,0\n"
	    ,name_rng[which_rng]);
  fclose(plotcommandf);

    //initialize pipe, give gnuplot initial commands through cdf.output
    pipef = popen("/usr/bin/gnuplot -persist cdf.output -", "w");
    fflush(pipef);          // force it through the pipe
    /* printf("pipe opened to gnuplot\n"); */
  }

  // Insert information into output file

  if (OUTPUT_PRINT) {
    fprintf(dataf,"# Random number generator = %s\n"
	    "# CHOOSE_TIMESEED = %d\n"
	    "# CONST_SEED = %ld\n"
	    "# Number of data points generated = %d\n#\n"
	    ,name_rng[which_rng],CHOOSE_TIMESEED,CONST_SEED,NPOINTS_MAP);
  };

  // Initialise random number generator:
  xp = genrand(which_rng, &isd);

  // Fill up first line
  xp = genrand(which_rng, &isd);
  yp = genrand(which_rng, &isd);
  zp = genrand(which_rng, &isd);
  if (OUTPUT_PRINT) {fprintf(dataf,"%.4f %.4f %.4f\n",xp,yp,zp);};
  //printf("%.4f %.4f %.4f\n",xp,yp,zp);
  // Subsequent lines
  for(i=1; i < NPOINTS_MAP; i++) {
    xp = yp;
    yp = zp;
    zp = genrand(which_rng, &isd);
    if (OUTPUT_PRINT) { 
      fprintf(dataf,"%.4f %.4f %.4f\n",xp,yp,zp);
      //printf("%.4f %.4f %.4f\n",xp,yp,zp);
    };
  }

  fclose(dataf);
  if (GNUPLOT_ON) {
    xrotation = 60;   zrotation = 145;    linestyle = 2;
    //xrotation = 60;   zrotation = 130;    linestyle = 2;
    gnu_plot_3d(pipef, outname, xrotation, zrotation, linestyle);
    fclose(pipef);
  }
  return 0;
}


/*--------------------------------------------------------------------*/
/* General random number function, offering a choice of specific RNGs 
   This is obviously slow, so bypass this routine once you want 
   to do real number-crunching.
*/
/*--------------------------------------------------------------------*/

double genrand(int which_rng, long *ircall) {

  double random_number;
  long ir;

  ir = *ircall;

  /* Initialise RNG: */
  /* NB Seed must be negative for most generators! */
  if (rng_initialise) {

    rng_initialise = 0;
    if (CHOOSE_TIMESEED) {
      ir = time(0); 
      printf("USING TIME SEED  %ld\n",ir);
    } 
    else {
      printf("USING CONSTANT SEED  %ld\n",ir);
    }
    switch (which_rng) {
    case -1: {break;}
    case  0: {break;}            // ran0 needs no initialisation
      // negative numbers needed to initialise:
    case  1: {ir=-ir; random_number = double_ran1(&ir); break;}
    case  2: {ir=-ir; random_number = double_ran2(&ir); break;}
    case  3: {ir=-ir; random_number = double_ran3(&ir); break;}
    case  4: {ir=-ir; random_number = (double) ran4(&ir); break;}

      /* // How do we get the in-house RNG's into the pointer system? */
      /* //case  7: {static float (*randp) (long *) = &rand;  break;} */
      /* //case  8: {static float (*randp) (long *) = &random;  break;} */
      /* //case  9: {static float (*randp) (long *) = &drand48;  break;} */

    // No initialisation needed:
    case 10: {break;}
    case 11: {break;}
    }
  } 

  else {

    /* Get one random number from the RNG of choice: */
    switch (which_rng) {
    case -1: {break;}
    case  0: {random_number = double_ran0(&ir); break;}
    case  1: {random_number = double_ran1(&ir); break;}
    case  2: {random_number = double_ran2(&ir); break;}
    case  3: {random_number = double_ran3(&ir); break;}
    case  4: {random_number = (double) ran4(&ir); break;}
    case 10: {random_number = double_ranbinder(&ir); break;}
    case 11: {random_number = double_randu_ibm(&ir); break;}
    }

  }

  *ircall = ir;
  //printf("ir = %d\n",ir);
  return random_number;
}

/*--------------------------------------------------------------------*/
/*  
    Plot a number of plots (given by multiplot_no) in 2 dimensions.
    Column 1 is the x variable in every case,
    Column 2 is the y-variable for plot 0,
    Column 3 is the y-variable for plot 1 and so on.
*/
int gnu_plot_2d(FILE* pd, char *filename
		 , int xcol, int ycol, int linestyle)  
{
  char *the_plot;
  char *format;
  //int i;
  
  switch(linestyle) {
  case 0:			
    format = "lines\n";
    break;
  case 1:
    format = "points\n";
    break;
  case 2:
    format = "dots\n";
    break;
  case 3:
    format = "boxes\n";
    break;
  case 4:
    format = "linespoints\n";
    break;
  }
  
  
  //create mem each time since c string functions aren't friendly
  the_plot = (char*)malloc(64);
  
  // pipe each plot to gnuplot 
  sprintf(the_plot,
	  // COMMANDS PIPED TO GNUPLOT: 
	  //"set mouse; "
	  //"set size ratio -1; "
	  //"plot '%s' using 1:%d with %s\n",
	  "plot '%s' using %d:%d with %s\n",
	  filename, xcol, ycol, format);
  
  fprintf(pd, "%s", the_plot);
  fflush(pd);
  free(the_plot);
  return 0;
  
}

/*--------------------------------------------------------------------*/
/*  
    Three-dimensional plot of three columns of data
*/
//void plot3d(FILE* pd, char *filename, int linestyle)  {
int gnu_plot_3d(FILE* pd, char *filename, int xrot, int zrot, int linestyle)  {
  char *the_plot;
  char *format;
  
  switch(linestyle) {
  case 0:			
    format = "lines\n";
    break;
  case 1:
    format = "points\n";
    break;
  case 2:
    format = "dots\n";
    break;
  default:
    format = "dots\n";
    break;
  }
  
  //create mem each time since c string functions aren't friendly
  the_plot = (char*)malloc(64);
  
  sprintf(the_plot,
	  //"set mouse; "
	  "set view %d, %d; splot '%s' with %s\n",
	  xrot, zrot, filename, format);
  
  fprintf(pd, "%s", the_plot);
  fflush(pd);
  free(the_plot);
  return 0;
 
}

