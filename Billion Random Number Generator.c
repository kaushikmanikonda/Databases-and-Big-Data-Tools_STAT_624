#include <stdio.h>
#include "mkl_vsl.h"
#include "mkl.h"

#define NVECTOR 800L  // # different settings
#define REPS 1250000L      
int main() // hw8.c
{
/* Initializing r and average s */
double r[NVECTOR], s=0.0;
VSLStreamStatePtr stream;
int i, j;
// vslNewStream( &stream, VSL_BRNG_R250, 777 );    // A generalized feedback shift register generator
// vslNewStream( &stream, VSL_BRNG_MCG59, 777 );   // 59-bit multiplicative congruential generator
vslNewStream( &stream, VSL_BRNG_MT19937, 777 );    // Mersenne Twister pseudorandom number generator 

/* Generating random numbers*/
for ( i=0; i<REPS; i++ ) {

vdRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, stream, NVECTOR, r, 1.5, 3.0 );

for ( j=0; j<NVECTOR; j++ ) s += r[j];
}
/* Deleting the stream*/
vslDeleteStream( &stream );
printf( "Average = %f (n=%ld)\n", s/(REPS*NVECTOR), REPS*NVECTOR );

return 0;
}