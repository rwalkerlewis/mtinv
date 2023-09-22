/*
#include <stdio.h>
#include <stdlib.h>

void scale_data( float *x, int n, float a )
{
	int i;
	
	fprintf( stdout, "%s: %s: n=%d a=%g\n",
		__FILE__, __func__, n, a );

	fflush(stdout);

	for( i = 0; i < n; i++ )
	{
	  x[i] = x[i] * a;
	}
}
*/

void scale_data( float *x, int n, float a )
{
	int i;
	for( i = 0; i < n; i++ ) x[i] *= a;
}
