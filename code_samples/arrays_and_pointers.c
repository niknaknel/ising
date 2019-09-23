/* ----------------------------------------------------------- */
/*   Example to show operations of pointers and how to use them to
     pass values and array values back to the calling program */
/* ----------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// prototypes:
void swap_two_numbers(int *px, int *py);
void swap_arrays(int *apx, int *apy, int *size_of_array);

/* ----------------------------------------------------------- */
int main() 
{
  int a, b;
  int size_of_array = 2;
  int arr_a[size_of_array], arr_b[size_of_array];

  // SWAP OF TWO INTEGERS:
  a = 1;
  b = 4;
  printf("Before swap: a = %d\tb = %d\n",a,b);
  swap_two_numbers(&a,&b);
  printf("After swap:  a = %d\tb = %d\n\n",a,b);


  // SWAP TWO ARRAYS:
  arr_a[0] = 200;
  arr_a[1] = 201;
  arr_b[0] = 300;
  arr_b[1] = 301;
  printf("Before swap: arr_a = %d %d\t\tarr_b = %d %d\n"
	 ,arr_a[0],arr_a[1],arr_b[0],arr_b[1]);

  // Note that we do not need to specify &arr_a, &arr_b below:
  swap_arrays(arr_a, arr_b, &size_of_array);

  printf("After swap:  arr_a = %d %d\t\tarr_b = %d %d\n"
	 ,arr_a[0],arr_a[1],arr_b[0],arr_b[1]);

  return 0;
}
/*--------------------------------------------------------------------*/
void swap_two_numbers(int *px, int *py) 
{
  int temp;
  temp = *px;
  *px = *py;
  *py = temp;
  return;
}
/*--------------------------------------------------------------------*/
// apx and apy are both arrays

// The pointers in the declaration point to apx[0] and apy[0]
// respectively
void swap_arrays(int *apx, int *apy, int *size_of_array)
{
  int temp;
  int i;
  // Below, *(apx + i) is a pointer to the content of an element of
  // array apx, ie to what is contained in apx[i]
  for (i=0; i < *size_of_array; i++) {
    temp = *(apx + i);
    *(apx + i) = *(apy + i);
    *(apy + i) = temp;
  }

  // The above can also be written as (try it):

  //  for (i=0; i < *size_of_array; i++) {
  //    temp = apx[i];
  //    apx[i] = apy[i];
  //    apy[i] = temp;
  //  }

  return;
}
/*--------------------------------------------------------------------*/

