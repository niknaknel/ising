/*

  Random number generator that is deliberately bad, used for testing

  Do not EVER use this for serious work

  Source: K Binder and DW Heermann 
  "Monte Carlo Simulation in Statistical Physics"

 */  
double double_ranbinder(long *idum)
{
  long mult = 1277;
  long modulo = 131072; //  2 to the power 17
  long tmp;
  double rmodulo;
  
  tmp = (long) *idum;
  rmodulo = (double) modulo;
  tmp = mult*tmp;
  *idum = (long) (tmp) % modulo;
  
  return ((double) *idum/ (double) rmodulo);
}
	
