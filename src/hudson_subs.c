#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*** Sample call ***/
/***
	float t,k,x,y;
	int hudson_transform(float,float,float*,float*);
	
	if( hudson_transform(t,k,&x,&y) != 0 )
        {
          fprintf(stderr, "ERROR, input out of range, skipping row.\n" );
          continue;
        }
	printf("%g %g %s", x,y,rec );
***/


int hudson_transform( float t, float k, float *xc, float *yc )
{
        float x,y;
        float sign(float,float);
                                                                                                                                                          
        if( fabs(t) > 1 || fabs(k) > 1 ) return 1;
                                                                                                                                                          
        if( fabs(t) == 1 && fabs(k) == 1 ) t = 0;
                                                                                                                                                          
        if( ( k >= 0 && t <= 0 ) || ( k <= 0 && t >= 0 ) )
        {
                x = t * ( 1 - fabs(k) );
                y = k;
        }
        else
        {
                x = sign(1,k) / ( 1/(fabs(t) * (1-fabs(k))) - 0.5 );
                y = k * ( 1 + fabs(x)/2 );
                if( (y/x) < 0.25 )
                {
                        if( k != 0 )
                                y = ( 1/(1/k - (2*sign(1,k))));
                        else
                                y = 0;
                        x = t * ( 1 + fabs(y));
                }
        }
        *xc = x;
        *yc = y;
        return 0;
}
                                                                                                                                                          
float sign(float a,float b)
{
        if( b>= 0 ) return fabs(a);
        if( b < 0 ) return -fabs(a);
	fprintf(stderr,"error in sign()\n" );
	exit(-1);
}
