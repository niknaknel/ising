/*--------------------------------------------------------------------*/
/* 
   Generate variates distributed according to standardised
   gaussian (mean zero, variance one) using Box-Muller method.

   Original derived from gasdev.c of Press et al.

   All important quantities have been made double or long

   double_gasdev needs double_ran1 to run.
*/
/*--------------------------------------------------------------------*/
#include <math.h>

double double_gasdev(long *idum)
{
  //double double_ran1(long *idum);
  static long iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if (*idum < 0L) iset = 0L;     // for reinitialisation
  if  (iset == 0) {             // no extra variate handy, so...
    // pick two uniform rv's in the square extending from (-1,-1) to (1,1)
    do {
      v1=2.0*double_ran1(idum)-1.0;
      v2=2.0*double_ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
      // is variate inside unit circle? If not, try again
    } while (rsq >= 1.0L || rsq == 0.0L);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    // set flag: we have a second variate ready
    iset=1;
    return v2*fac;
    // Extra variate available, use and reset flag:
  } else {
    iset=0L;
    return gset;
  }
}

