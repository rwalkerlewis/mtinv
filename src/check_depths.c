#include <stdio.h>
#include <stdlib.h>
#include "../include/mt.h"
char progname[128];

void check_iso_depth( FixISOZ *myfixisoz, float *z, int nz, int verbose )
{
	int is_my_z_valid( float, float *, int );

	if( ( myfixisoz->indexz = is_my_z_valid( myfixisoz->z, z, nz )) == -999 )
	{
		fprintf( stderr, "%s: %s: %s: ERROR! FixISOZ.z=%g is ***NOT*** a valid depth in ginv file\n",
	  	  progname, __FILE__, __func__, myfixisoz->z );
		exit(-1);
	}
	else
	{
		if( verbose )
		{
			fprintf( stderr,
			  "%s: %s: %s: FixISOZ.z=%g FixISOZ.index=%d is a valid depth in ginv files\n",
			    progname, __FILE__, __func__, myfixisoz->z, myfixisoz->indexz );
		}
	}
}

void check_depth( float FixMyZ, int *FixMyiz, float *z, int nz, int verbose )
{
	int is_my_z_valid( float, float *, int );

	if(verbose)
	{
		fprintf( stderr, "%s: %s: %s: searching for FixMyZ=%g nz=%d\n",
			progname, __FILE__, __func__, FixMyZ, nz );
	}

	if( ( *FixMyiz = is_my_z_valid( FixMyZ, z, nz )) == -999 )
	{
		fprintf( stderr,
			"%s: %s: %s: force_z=%g is ***NOT*** a valid depth in ginv files\n",
			progname, __FILE__, __func__, FixMyZ );
		exit(-1);
	}
	else
	{
		if(verbose)
		{
			fprintf( stderr, "%s: %s: %s: fixz=%g iz=%d is a valid depth in ginv files\n",
				progname, __FILE__, __func__, FixMyZ, *FixMyiz );
		}
	}
}

int is_my_z_valid( float myz, float *z, int nz )
{
        int iz, myiz;
        for( iz=0; iz<nz; iz++ )
        {
                if( myz == z[iz] )
                {
                        myiz = iz;
                        return myiz;
                }
        }
        return -999;
}
