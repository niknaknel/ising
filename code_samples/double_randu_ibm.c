/*

  Old random number generator once used as standard by IBM Corporation
  It was found to be bad and should not be used for serious work.

  Reference: Gentle, "Random Number Generators and Monte Carlo Methods"
  page 19

 */  
double double_randu_ibm(long *idum)
{
  long mult = 65539;
  long modulo = 2147483648; /* 2 to power 31, larger than RAND_MAX*/
  long tmp;

  tmp = *idum;
  tmp = mult*tmp;
  *idum = (tmp) % modulo;
  return ((double) *idum/(double) modulo);
}
	
