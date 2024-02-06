#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/*** mkdirp ***/
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

char progname[256];

/*********************************************************************/
/*** executes system call to run gmt grd2xyz lune.grd              ***/
/*** reads lon, lat, %VR from lune.grd and computes spherical area ***/
/*** of area covered contour %VR greater than contour 80%          ***/
/*********************************************************************/

int main(int ac, char **av)
{
	char rec[256];

	float lat, lon, data, dx, dy, sum_norm = 0, sum_total = 0, total_area = 0, area = 0;
	float contour = 80;
	float contour_norm = 0.95;
	int i, n = 1, count = 0;
	FILE *fp;
	char currentdir[128];
	char grd_filename[256];
	char command_line[512];
	float *x, *y, *z;
	float xmin, xmax, ymin, ymax, zmin, zmax, znorm = 0, d2r;
	int verbose = 0;

	int setpar(int ac, char **av);
	int mstpar( char *, char *, void * );
        int getpar( char *, char *, void * );
	void endpar();
	void findminmax( float *z, int n, float *zmin, float *zmax );
	float estimate_spacing( float *x, int n );

	d2r = M_PI / 180;

	strcpy( progname, av[0] );
	strcpy( grd_filename, "lune.grd" );

	setpar(ac,av);
	getpar( "f", "s", grd_filename );
	getpar( "contour_norm", "f", &contour_norm );
	getpar( "contour", "f", &contour );
	getpar( "verbose", "b", &verbose );
	endpar();

	getcwd( currentdir, 128 );
	chdir( currentdir );

	if( (fp = fopen( grd_filename, "rb" )) == NULL )
	{
		fprintf( stderr, "%s: %s: %s: Error cannot open file %s for reading\n",
			progname, __FILE__, __func__, grd_filename );
		exit(-1);
	}
	fclose(fp);

	sprintf( command_line, "gmt grd2xyz %s > lune.out", grd_filename );
	fprintf( stderr, "%s: %s: %s: executing command %s\n",
                progname, __FILE__, __func__, command_line );
	system( command_line );

	if( (fp = fopen( "lune.out", "r" )) == NULL )
	{
		fprintf( stderr, "%s: %s: FATAL ERROR: cannot open file lune.out for reading.\n",
			__FILE__, __func__ );
		exit(-1);
	}

	x = calloc( 1, sizeof(float) );
        y = calloc( 1, sizeof(float) );
        z = calloc( 1, sizeof(float) );
	i  = 0;
	count = 0;
	while( fgets( rec, 256, fp ) != NULL )
	{
		x = realloc(x,(i+1)*sizeof(float));
                y = realloc(y,(i+1)*sizeof(float));
                z = realloc(z,(i+1)*sizeof(float));
		/***                     lon      lat      z value ***/
		sscanf( rec, "%f %f %f", &(x[i]), &(y[i]), &(z[i]) );
		i++;
		count++;
	}
	n = i;
	fclose(fp);
	system( "/bin/rm -f lune.out" );

	findminmax( z, n, &zmin, &zmax );
	findminmax( x, n, &xmin, &xmax );
	findminmax( y, n, &ymin, &ymax );

	dx = estimate_spacing( x, n );
	dy = estimate_spacing( y, n );

	fprintf( stderr, "%s: %s: %s: n=%d dx=%g dy=%g xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g %s\n",
                progname, __FILE__, __func__,
                n, dx, dy,
                xmin, xmax,
                ymin, ymax,
                zmin, zmax,
                currentdir );

	if(verbose)
	{
	  fprintf( stderr, "%s: %s: %s: n=%d dx=%g dy=%g xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g %s\n",
		progname, __FILE__, __func__, 
		n, dx, dy,
		xmin, xmax, 
		ymin, ymax, 
		zmin, zmax, 
		currentdir );
	}

	fp = fopen( "lune_norm.xyz", "w");
	sum_norm = 0;
	sum_total = 0;
	total_area = 0;
	for( i = 0; i < n; i++ )
	{
		znorm  = ( z[i] + fabs(zmin) ) / ( zmax - zmin );
		fprintf( fp, "%g %g %g\n", x[i], y[i], znorm );

		area = dx * dy * d2r * cos( y[i] * d2r );
		if( znorm >= contour_norm ) sum_norm += area;
		if( z[i]  >= contour      ) sum_total += area;
		total_area += area;
	}
	fclose(fp);

	sprintf( command_line, "gmt xyz2grd -R%g/%g/%g/%g -I%g/%g -Glune_norm.grd lune_norm.xyz", 
		xmin, xmax, ymin, ymax, dx, dy );
	
	fprintf( stderr, "%s: %s: %s: executing command %s\n", 
		progname, __FILE__, __func__, command_line );

        system( command_line );

	
	fprintf( stdout, "count= %d total_area= %g contour= %g countour_norm= %g area= %g ratio= %g norm-area= %g norm-ratio= %g %s\n",
		count, total_area, contour, contour_norm, sum_total, (sum_total/total_area), sum_norm, (sum_norm/total_area), currentdir );

	free(x);
	free(y);
	free(z);
	return 0;
}

void findminmax( float *z, int n, float *zmin, float *zmax )
{
	int i;
	*zmin = z[0]; 
	*zmax = z[0];
	for( i = 1; i < n; i++ )
	{
		if( z[i] > *zmax ) *zmax = z[i];
		if( z[i] < *zmin ) *zmin = z[i];
	}
}

float estimate_spacing( float *data, int n )
{
	float delta, max = 0; 
	int i;
	for( i = 0; i < n-1; i++ )
	{
		delta = fabs( fabs(data[i+1]) - fabs(data[i]) );
		if(delta>max) max = delta;
	}
	return max;
}
