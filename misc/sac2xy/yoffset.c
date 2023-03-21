#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int ac, char **av)
{
	float x, xinc, x0,x1;
	float t, t0, vel;

	int setpar(int, char **);
        int getpar(), mstpar();
        void endpar();

	setpar(ac,av);
	mstpar("x0", "f", &x0);
	mstpar("x1", "f", &x1);
	mstpar("t0", "f", &t0 );
	mstpar("xinc", "f", &xinc );
	mstpar("vkmps", "f", &vel);
	endpar();

	fprintf( stdout, "> velocity km/s = %g \n", vel );
	for( x = x0 ; x < x1; x += xinc )
	{
		t = t0 + x / vel;
		fprintf( stdout, "%g %g\n", t, x );
	}
}
