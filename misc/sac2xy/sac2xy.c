/***
	Original found in /Users/ichinose1/Work/mtinv.v3.0.5/misc/sac2xy
***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sac.h"

char progname[256];

int main( int ac, char **av)
{
	Sac_Header s;
	int i, xreset=FALSE;
	float *x, *y;
	int norm = 0, xdeg = 0, xkm = 0;
	float ymax, scale = 1.0, ydist;
	int setpar(int, char **);
	int getpar();
	void endpar();
	int verbose = 0;
	int normfac = 0;
	int ifixZeroDiv = 0;
	float xmax;
	float tshift = 0;

	sprintf( progname, "%s", av[0] );

	setpar( ac, av );
	getpar( "xreset", "b", &xreset);
	getpar( "tshift", "f", &tshift );
	getpar( "norm",   "b", &norm );
	getpar( "sc",     "f", &scale );
	getpar( "xdeg",   "b", &xdeg );
	getpar( "xkm",     "b", &xkm );
	getpar( "normfac", "b", &normfac );
	getpar( "verbose", "b", &verbose );
	getpar( "fixZeroDiv", "b", &ifixZeroDiv );
	endpar();

	if( xdeg && xkm )
	{
		fprintf( stderr, "%s: error only specify one xdeg=%d or xkm=%d\n",
			progname, xdeg, xkm );
		exit(-1);
	}

	if(verbose)
	{
	  fprintf( stderr, "%s: xreset=%d norm=%d sc=%g xdeg=%d xkm=%d ifixZeroDiv=%d\n",
		progname, xreset, norm, scale, xdeg, xkm, ifixZeroDiv );
	}

/*** read the data and allocate memory ***/

	fread(&s, sizeof(Sac_Header), 1, stdin);
	x = (float *)calloc(s.npts, sizeof(float));
        y = (float *)calloc(s.npts, sizeof(float));
        fread(&y[0], s.npts*sizeof(float), 1, stdin);

	if(verbose) fprintf( stderr, "%s: npts=%d delta=%g b=%g sta=%s\n",
		progname, s.npts, s.delta, s.b, s.kstnm );

/*** reset the x-axes to 0 and add the tshift - time shift ***/

	s.b = s.b + tshift;
	/* if( s.iftype == ITIME ) s.b=0; */
	if( xreset == TRUE ) s.b=0;
	xmax = s.b + ( s.delta * (float)s.npts );

/*** find the normalizing factor ***/
	if( norm || normfac )
	{
		ymax = fabs(y[0]);
		for( i=1; i<s.npts; i++ )
		{
			if( fabs(y[i]) > ymax ) ymax = fabs(y[i]);
		}
	}
	
	if( normfac )
	{
		fprintf( stdout, "%6.3e %e %g", ymax, 1.0/ymax, xmax );
		exit(0);
	}

/*** write the results ***/

	fprintf( stdout, ">  %-8.8s %-8.8s %-8.8s %-8.8s %g\n",
		s.knetwk, s.kstnm, s.khole, s.kcmpnm, s.dist );

        for( i=0; i<s.npts; i++)
	{
		x[i]=s.b + (float)i*s.delta;

	/*** prevent division by zero ***/
		if(ifixZeroDiv)	
		{
		 if( x[i] == 0 ) x[i]=1.0e-9; 
		 if( y[i] == 0 ) y[i]=1.0e-9; 
		}

		if( xdeg )
			ydist = s.gcarc;
		else
			ydist = 0;

		if( xkm ) 
			ydist = s.dist;
		else
			ydist = 0;

		if( norm )
			fprintf( stdout, "%g %g\n", x[i], scale*(y[i]/ymax) + ydist );
		else
                	fprintf( stdout, "%g %g\n", x[i], scale*y[i] + ydist );
	}

	if(verbose) fprintf( stderr, "%s: exiting \n", progname );
}
