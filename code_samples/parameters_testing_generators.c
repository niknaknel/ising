/* ----------------------------------------------------------- */
/*   
     PARAMETERS FOR USE IN testing_generators.c

     HCE   2014
*/
/* ----------------------------------------------------------- */

//  SET THESE SWITCHES AS DESIRED:
#define AVOID_MENU   0      // switches main menu on/off
#define GNUPLOT_ON   1       // switching direct plotting on or off
#define NPOINTS_MAP 20000   // points to be plotted in map2d, map3d
//#define NPOINTS_MAP 2e8   // points to be plotted in map2d, map3d
#define NRNOS 1e8           // no. of random variates for other tests
#define OUTPUT_PRINT 1      // generate ASCII output files (or not)

// HISTOGRAM RELATED:
#define HIST_NBINS  200
#define HIST_LOWER  0.0L
#define HIST_UPPER  1.0L
#define NUSUM  1           // number of uniform variates in sum
#define XLIMIT  0.0001     // limit for conditional distribution
#define NSHOW   10000000   // on-screen progress report

