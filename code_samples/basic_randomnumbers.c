/*
-----------------------------------------------------------------

     BASICS ON RANDOM NUMBER GENERATORS IN C               

There are four intrinsic sets of RNGs on our system, along with their
seed initialisation routines. In order to test RNG's, we use
however routines provided by Press et al, called

1. ran0, ran1, ran2, ran3, ran4 

and their double equivalents. In addition there are

2. purposely bad generators for the purposes of testing such as
    ranbinder.c
    randu_ibm.c

3. others to be found on the internet,

There are also built-in generators. (For details on the built-in
generators, type in the shell man 3 <randversion> or in konqueror type
man:rand or man:random or man:drand48 and select the appropriate
page.) We have:

4. rand, srand (preferred within C, yield int output)

5. random, srandom (BSD implementation, yield double integer)

6. 48-bit generators:

       drand48,erand48   double precision float with/without seeding
       lrand48,nrand48   long integers with/without seeding
       mrand48, jrand48  signed integers with/without seeding

       srand48, seed48, lcong48 = seed initialisation routines

       *ran48* are called "obsolete" in SysV but appear still to be
       useful.

The idea of this program is to try out these to see whether they
are on the system and how they work.

-----------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
#include <ctype.h>
#include "ran0.c"
#include "ran1.c"
// Note that definitions do not end with a semicolon:
#define RNOS      4
#define CONSTSEED 12321432

/*-----------------------------------------------------------*/
int main()
{
  /*  long const seedconst=12321432;*/
  long cseed;
  int i, td;
  int icase, ret = 1;
  // Types: float and double, and integers int and long
  float    ans;
  double dians;
  long    ians;
  long timeseed;
  unsigned short xsubi[3];

  icase=1;
  printf("See /usr/include/stdlib.h for details of the "
	 "built-in generators\n");
  //while (icase > 0) {
  while (ret > 0 && icase > 0) {
    printf("\n--------------------------------------------------\n"
	   "Type\t 0 to stop\n"
	   "\t 1 for seed and RAND_MAX information \n"
	   "\t 2 for usage of random()\n"
	   "\t 3 for usage of rand()\n"
	   "\t 4 for comparison of division with various types\n"
	   "\t 5 for usage of ran0 \n"
	   "\t 6 for usage of ran1\n"
	   "\t 7 for an example of a time-dependent seed\n"
	   "\t 8 for usage of srand48\n"
	   "\t 9 for usage of drand48\n"
	   "\t10 for usage of erand48 (unclear if this is correct)\n"
	   "--------------------------------------------------\n"
	   );
    ret = scanf("%d", &icase);
    //ret = scanf("%[0123456789]", &icase);
    td = isdigit(ret); 
    //printf("\n\n\n%d\n\n\n",td);
    
    if (icase < 0 || icase > 10) { 
      printf("Wrong choice of icase\n\n");
      exit(1);
    }
    

    switch (icase) {
      case 0: 
	exit(0);

      case 1:
	printf("RAND_MAX = %d\t--- "
	       "Maximum value of integer random number generators "
	       "(equal to 2^31 - 1)"
	       ,RAND_MAX);
	printf("\n\nConstant seed used = \t%d\t--- "
	       "Use only seeds with at least 8 digits!",CONSTSEED);
	i = (int) time(0);
	printf("\n\nInteger from time(0) function:\t%d\t--- "
	       "Use for time-dependent seed\n"
	       "Type 1 repeatedly to see how the time integer changes",i);
	break;
	//--------------------------------
      case 2: 
	cseed = CONSTSEED;
	printf("random() using constant seed  %ld:\n",cseed);
	// Initialisation: see /usr/include/stdlib.h for details
	srandom((unsigned)cseed);
	//srandom(1);
	for( i=0; i <= RNOS; i++) {
	  ians=random();
	  printf("%ld  ",ians);
	}
	break;
	//--------------------------------
      case 3:
	cseed = CONSTSEED;
	printf("rand() using constant seed  %ld (showing that"
	       " rand and random use the same engine):\n",cseed);
	srand((unsigned) cseed);
	for( i=0; i <= RNOS; i++) {
	  ians=rand();
	  printf("%ld  ",ians);
	}
	break;
	//--------------------------------
      case 4:
	cseed = CONSTSEED;
	printf("Obtaining floats from rand():\nusing division with INTEGER type\n");
	srand((unsigned) cseed);
	for( i=0; i <= RNOS; i++) {
	  dians = rand()/ RAND_MAX;
	  printf("%.15f  ",dians);
	}

	cseed = CONSTSEED;
	printf("\nusing division with DOUBLE type\n");
	srand((unsigned) cseed);
	for( i=0; i <= RNOS; i++) {
	  dians=(double) rand()/ (double) RAND_MAX;
	  printf("%.15f  ",dians);
	}

	cseed = CONSTSEED;
	printf("\nusing division with FLOAT type:\n"
	       "Note that the numbers below agree with the above double "
	       "only to 8 digits!\n");
	srand((unsigned) cseed);
	for( i=0; i <= RNOS; i++) {
	  ans=(float) rand()/ (float) RAND_MAX;
	  printf("%.15f  ",ans);
	}
	break;
	//--------------------------------
      case 5:
	cseed = CONSTSEED;
	printf("ran0 using constant seed   %ld:\n",cseed);
	for( i=0; i <= RNOS; i++) {
	  ans=ran0(&cseed);
	  //printf("%d %f\n",i,ans);
	  printf("%f  ",ans);
	}
	break;
	//--------------------------------
      case 6:
	cseed = - CONSTSEED;
	printf("NOTE that ran1 needs to be initialised with"
	       "a NEGATIVE seed:\n");
	printf("ran1 using constant seed   %ld:\n",cseed);
	ans = ran1(&cseed);
	for( i=0; i <= RNOS; i++) {
	  ans=ran1(&cseed);
	  //printf("%d %f\n",i,ans);
	  printf("%f  ",ans);
	}
	break;
	//--------------------------------
      case 7:
	timeseed = time(0);
	printf("ran0 time-dependent seeding with timeseed  %ld\n",timeseed);
	for( i=0; i <= RNOS; i++) {
	  ans=ran0(&timeseed);
	  printf("%f  ",ans);
	}
	break;
	//--------------------------------
      case 8:
	cseed = CONSTSEED;
	printf("lrand48() with seed %d\n",CONSTSEED);
	srand48((long) cseed);
	for( i=0; i <= RNOS; i++) {
	  ians=lrand48();
	  printf("%ld  ",ians);
	}
	break;
	//--------------------------------
      case 9:
	cseed = CONSTSEED;
	printf("\n\ndrand48() with seed %d\n",CONSTSEED);
	srand48((long) cseed);
	for( i=0; i <= RNOS; i++) {
	  dians=drand48();
	  printf("%.15f  ",dians);
	}
	break;
	//--------------------------------
      case 10:
	printf("erand48() generator (unclear whether usage is correct...):\n");
	xsubi[1] = 12;
	xsubi[2] = 98;
	xsubi[3] = 43;
	for( i=0; i <= RNOS; i++) {
	  dians=erand48(xsubi);
	  printf("%.15f  ",dians);
	}
	break;
	//--------------------------------
	// case 12: {
	//break;
	// }
      default:
	printf("ERROR: No result defined for this icase\n"); exit(1);
	break;
      }
    //printf("\n");
  }

  return 0;
}

