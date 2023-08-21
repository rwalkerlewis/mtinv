/***
	Original found in /Users/ichinose1/Work/mtinv.v3.0.5/misc/sac2xy
***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>

#include "sac.h"
#include "mytime.h"

char progname[256];

int main( int ac, char **av)
{
	Sac_Header s;
	float *spec;
	float *spec_corrected;

	int setpar(int, char **), getpar(), mstpar();
	void endpar();
	int verbose = 0;

	int i;
	int ilog10 = 1;

	float R0 = 100; /***geometrical spreading factor***/
	float eta = 0.5; /***shear waves***/
	float GS;
	float epidistkm;
	float km2m = 1000;
	void GeometricalSpreading( float epidistkm, float R0, float eta, float *GS );

	float Q0=250, qn=1.0, U0=3500;
	float Attenuation( float epidistkm, float Q0, float qn, float U0, float freq );
	float f, atten;

/*** begin main ***/
	sprintf( progname, "%s", av[0] );

	setpar( ac, av );
/***geometrical spreading***/
	mstpar( "R0", "f", &R0 );
	mstpar( "eta", "f", &eta );

/*** attenuation ***/
	mstpar( "Q0", "f", &Q0 );
        mstpar( "qn", "f", &qn );
	mstpar( "U0", "f", &U0 );

/*** use log10 amplitude ***/
	mstpar( "log10", "b", &ilog10 );

	endpar();

/*** end getpar ***/

/*** read the data and allocate memory ***/

	fread(&s, sizeof(Sac_Header), 1, stdin);
	spec = (float *)calloc(s.npts, sizeof(float));
	spec_corrected = (float *)calloc(s.npts, sizeof(float));
        fread(&spec[0], s.npts*sizeof(float), 1, stdin);

	epidistkm = s.dist;

	GeometricalSpreading( epidistkm, R0, eta, &GS );

	fprintf( stderr, "%s: nspec=%d df=%g f0=%g sta=%s distancekm=%g GS=%e\n",
                progname,
                s.npts, s.delta, s.b, s.kstnm, epidistkm, GS );

	for( i = 0; i < s.npts; i++ )
	{
		if(ilog10)
			spec_corrected[i] = -log10(GS) + spec[i];
		else
			spec_corrected[i] =  spec[i]/GS;
	}

	for( i = 0; i < s.npts; i++ )
	{
		f = s.b + (s.delta * (float)i);
		atten = Attenuation( epidistkm, Q0, qn, U0, f );
		if(ilog10)
			spec_corrected[i] = -log10( atten ) + spec_corrected[i];
		else
			spec_corrected[i] =  spec_corrected[i]/atten;
	}

	fwrite( &s, sizeof(Sac_Header), 1, stdout );
	fwrite( &spec_corrected[0], s.npts*sizeof(float), 1, stdout );

	if(verbose) fprintf( stderr, "%s: exiting \n", progname );
}

void GeometricalSpreading( float epidistkm, float R0, float eta, float *GS )
{
	float km2m = 1000;
	if( epidistkm > R0 )
		*GS = 1.0/(R0 * km2m) * pow( ((R0 * km2m)/(epidistkm * km2m)), eta );
	else
		*GS = pow( ( epidistkm * km2m ), -eta );
}

float Attenuation( float epidistkm, float Q0, float qn, float U0, float freq )
{
	float atten, Q;
	Q = Q0 * pow( freq, qn );
	atten = exp( -M_PI * freq * epidistkm / ( U0 * Q ) );
	return atten;
}
