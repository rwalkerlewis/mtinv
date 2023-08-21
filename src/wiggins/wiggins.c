#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/sac.h"

char progname[256];

int main( int ac, char **av )
{
	float *data;
	Sac_Header s;
	char FileName[128];
	int i;
	float new_dt;
	int new_nt;
	float *z;
	int verbose = 0;

/*** fuction prototypes ***/

	void wrtoldsac( char *FO, Sac_Header *s, float *data );

	float *readsac( Sac_Header *s, char *filename, int verbose );

	void interpolate_wiggins( float *y, int npts, float delta, float b, float *z, int new_nt, float new_dt );

	void set_sac_minmax( Sac_Header *s, float *data );

	int setpar(int ac, char **av), mstpar(),getpar();
	void endpar();

/*** begin main ***/

	setpar(ac,av);
	mstpar("f",  "s", FileName);
	mstpar("nt", "d", &new_nt);
	mstpar("dt", "f", &new_dt);
	getpar("verbose", "b", &verbose);
	endpar();

	data = readsac( &s, FileName, verbose );

	if(verbose) fprintf( stdout, "%s: %s: %s nt=%d dt=%g\n", 
		__FILE__, __func__, FileName, s.npts, s.delta );

	z = (float *)calloc( new_nt, sizeof(float) );

	interpolate_wiggins( data, s.npts, s.delta, s.b, z, new_nt, new_dt );

	s.npts = new_nt;
	s.delta = new_dt;
	
	set_sac_minmax( &s, data );

	wrtoldsac( "xint.sac", &s, z );

	return 0;
}
