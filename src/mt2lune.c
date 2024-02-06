/***
	reads output from CSS database using MT and Focal tables to convert mt to lune parameters
***/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#include "../include/nrutil.h"
#include "../include/mt.h"

char progname[128];

typedef struct {
	int label;
	long oridin;
	char mt_type[12]; /* FULL or DEV */
	float momag;
	float m0;
	int pdc, pclvd, piso;
	float vred;
	float mxx, myy, mzz, mxy, mxz, myz;
	float mrr, mtt, mff, mrt, mrf, mtf;
	char auth[32];
	float mtdepth;
	long evid;
	long oridout;
	long mtid;
	long fpid;
	long magid;
	float lat, lon;
	float lune_lat, lune_lon;
	char ot_str[32];
	char lddate_str[32];
	char Region[64];
	float str1, dip1, rak1, str2, dip2, rak2;

	float eval_full[4]; /* 1,2,3 */

	float Mtotal, Miso, Mdev, Mclvd, Mdc;
	float k, T, f_factor;
	float pdc2, pclvd2, pdev2, piso2;
} DBrow;

int main(int ac, char **av)
{
	DBrow *db;
	int i, nrec = 0;
	char dbfilename[256];
	DBrow  *load_dbrows( char *dbfilename, char *MT_coordinate_system, int *nrec, int verbose );
	DBrow *load_dbrows2( char *dbfilename, char *MT_coordinate_system, int *nrec, int verbose );

	void    write_dbrows( DBrow *db, int nrec );
	void   write_results( DBrow *db, int nrec );
	void  write_MTLoader( DBrow *db, int nrec );
	void       write_gmt( DBrow *db, int nrec );
	void write_dbrows_csv( DBrow *db, int nrec );

	float eig[3] = {0,0,0}; /* eig[0,1,2] */
	float lat = 0, lon = 0;
	int verbose = 0;
	void  eig2lune( float *eig, float *latitude, float *longitude, int verbose );
	char MT_coordinate_system[32]; 

	void copyMT2DB( MomentTensor *Ma, DBrow *db );
	float Mw;
	int igmt = 0;

/* c=Cartesian(x,y,z) s=spherical(r,t,f) */
/*
	mrr mtt mff mrt  mrf  mtf 
	Mzz Mxx Myy Mzx -Mzy -Mxy
*/
	MomentTensor Ma;
	void mt2eig_version2( MomentTensor Ma, float *eig, int verbose );
	void eig2iso_version2( MomentTensor Ma, float *eig, 
		float *Mtotal, float *Mdev, float *Miso, float *Mclvd, float *Mdc,
		float *k, float *T, float *f_factor, 
		float *pdc2, float *pclvd2, float *pdev2, float *piso2, 
		int verbose );

	int setpar(int ac, char **av),getpar(),mstpar();
	void endpar();

/***** begin main ******/
	strcpy(progname,av[0]);

	setpar(ac,av);
	mstpar("sys","s",MT_coordinate_system );
	mstpar("f", "s", dbfilename );
	getpar("verbose","b",&verbose);
	getpar("gmt", "b", &igmt );
	endpar();
	
	/* db = load_dbrows( dbfilename, MT_coordinate_system, &nrec, verbose ); */
	db = load_dbrows2( dbfilename, MT_coordinate_system, &nrec, verbose );

	for( i = 0; i < nrec; i++ )
	{
		copyMT2DB( &Ma, &(db[i]) );

	/*** convert eigenvalues ***/

		mt2eig_version2( Ma, eig, verbose );

		db[i].eval_full[1] = eig[0];
		db[i].eval_full[2] = eig[1];
		db[i].eval_full[3] = eig[2];
		
	/*** compute lune ***/

		eig2lune( eig, &(db[i].lune_lat), &(db[i].lune_lon), verbose );

		eig2iso_version2( Ma, eig,
			&(db[i].Mtotal), &(db[i].Mdev), &(db[i].Miso), &(db[i].Mclvd), &(db[i].Mdc),
			&(db[i].k), &(db[i].T), &(db[i].f_factor),
			&(db[i].pdc2), &(db[i].pclvd2), &(db[i].pdev2), &(db[i].piso2), verbose );

	}

/*** cross check m0, momag(Mw), and Mij ****/
/*** cross check redecompose Mij->piso,pclvd,pdc with original piso,pclvd,pdc ***/

	for( i = 0; i < nrec; i++ )
	{
		Mw = ( log10(db[i].Mtotal) + 7 )/1.5 - 10.73;

		if( fabs( Mw - db[i].momag ) > 0.01 )
		{
			db[i].mxx = db[i].mxx * 10;
			db[i].myy = db[i].myy * 10;
			db[i].mzz = db[i].mzz * 10;
			db[i].mxy = db[i].mxy * 10;
			db[i].mxz = db[i].mxz * 10;
			db[i].myz = db[i].myz * 10;

			copyMT2DB( &Ma, &(db[i]) );

			mt2eig_version2( Ma, eig, verbose );

			db[i].eval_full[1] = eig[0];
			db[i].eval_full[2] = eig[1];
			db[i].eval_full[3] = eig[2];

			eig2lune( eig, &(db[i].lune_lat), &(db[i].lune_lon), verbose );

			eig2iso_version2( Ma, eig,
                        	&(db[i].Mtotal), &(db[i].Mdev), &(db[i].Miso), &(db[i].Mclvd), &(db[i].Mdc),
                        	&(db[i].k), &(db[i].T), &(db[i].f_factor),
                        	&(db[i].pdc2), &(db[i].pclvd2), &(db[i].pdev2), &(db[i].piso2), verbose );
		}
	}

/*** switch remap and check again plotmech and plotlune ***/

	/* write_dbrows( db, nrec ); */

	write_dbrows_csv( db, nrec );

	if(igmt)
	{
	  write_gmt( db, nrec );
	}
	else
	{
	 write_MTLoader( db, nrec );
	 write_results( db, nrec );
	}
}


void copyMT2DB( MomentTensor *Ma, DBrow *db )
{
/*** set the moment tensor structure ****/
/*** x_vector[1,...,6] = 0= 1=Mxx 2=Myy 3=Mxy 4=Mxz 5=Myz 6=Mzz ****/
                
                Ma->xx = db->mxx;
                Ma->yy = db->myy;
                Ma->xy = db->mxy;
                Ma->xz = db->mxz;
                Ma->yz = db->myz;
                Ma->zz = db->mzz;
                
                Ma->yx = Ma->xy;
                Ma->zx = Ma->xz;
                Ma->zy = Ma->yz;

	/*** deviatoric Ma->zz = -( Ma->xx + Ma->yy ); ***/
                
                Ma->mt[1][1] = Ma->xx;
                Ma->mt[1][2] = Ma->xy;
                Ma->mt[1][3] = Ma->xz;
                Ma->mt[2][1] = Ma->yx;
                Ma->mt[2][2] = Ma->yy;
                Ma->mt[2][3] = Ma->yz;
                Ma->mt[3][1] = Ma->zx;
                Ma->mt[3][2] = Ma->zy;
                Ma->mt[3][3] = Ma->zz;
}

/**

column myz format 9.99EEEE
column momag format 9.99
column m0 format 9.99EEEE

SELECT
        m.mxx, m.myy, m.mzz, m.mxy, m.mxz, m.myz,
        TO_TIMESTAMP( '1970-01-01', 'YYYY-MM-DD' ) + numtodsinterval( o.time, 'SECOND' ) as OT,
        o.lat,
        o.lon,
        o.depth,
        m.Momag,
        o.evid,
        m.mt_type,
        m.m0,
        m.var_red,
        m.piso,
        m.pclvd,
        m.pdc, 
        f.strike1,
        f.dip1,
        f.rake1,  
        m.auth,
        m.lddate,
        r.name as Region
from
        llnl.origin o,
        llnl.region_name r,
        llnl.moment m,
        llnl.focal_plane f
where
        f.fpid  = m.fpid and
        o.orid  = m.oridout and
        o.grn   = r.grn and
        m.mt_type = 'FULL'
order by m.auth, o.time, m.lddate;


**/

DBrow *load_dbrows2( char *dbfilename, char *MT_coordinate_system, int *nrec, int verbose )
{
	DBrow *db;
	FILE *fp;
	char rec[512];
	int i = 0;
	long longdumb;
	float mxx, myy, mzz, mxy, mxz, myz;
	void truncate_strings( char *str, int n );
	void fix_ot_str( char *str );
	
/**Begin ***/
	if( (fp = fopen( dbfilename, "r" )) == NULL )
	{
		fprintf( stderr, "%s: %s: cannot open file %s\n",
			__FILE__, __func__, dbfilename );
		exit(-1);
	}

	db = calloc( 1, sizeof(DBrow) );
	while( fgets( rec, 512, fp ) != NULL )
	{
		db = realloc( db, (i+1)*sizeof(DBrow) );

		db[i].label = 0;
		db[i].oridin = -1;

		/*             1  2  3  4  5  6  7  8  9 10 11  12 13 14 15 16 17 18 19 20 21 22 23 24 */
		sscanf( rec, "%e %e %e %e %e %e %s %f %f %f %f %ld %s %e %f %d %d %d %f %f %f %s %s %[^\n]s",
			&(db[i].mxx),
			&(db[i].myy), 
			&(db[i].mzz),
			&(db[i].mxy),
			&(db[i].mxz), /* 5 */
			&(db[i].myz),
			db[i].ot_str,
			&(db[i].lat),
			&(db[i].lon),    
			&(db[i].mtdepth), /* 10 */
			&(db[i].momag),
			&(db[i].evid),
                        db[i].mt_type,
                        &(db[i].m0),
			&(db[i].vred), /* 15 */
			&(db[i].piso),
			&(db[i].pclvd),
                        &(db[i].pdc),
                       	&(db[i].str1),
			&(db[i].dip1),  /* 20 */
			&(db[i].rak1),
                        db[i].auth,
                        db[i].lddate_str,
                        db[i].Region );  /* 24 */

                truncate_strings( db[i].Region, 64 );
		fix_ot_str( db[i].ot_str );

		if( strcmp( MT_coordinate_system, "c" ) == 0 || strcmp( MT_coordinate_system, "cartesian" ) == 0 )
                {
                  db[i].mrr =   db[i].mzz;
                  db[i].mtt =   db[i].mxx;
                  db[i].mff =   db[i].myy;
                  db[i].mrt =   db[i].mxz;
                  db[i].mrf =  -db[i].myz;
                  db[i].mtf =  -db[i].mxy;
                }
                else if( strcmp( MT_coordinate_system, "s" ) == 0 || strcmp( MT_coordinate_system, "spherical" ) == 0 )
                {
                          db[i].mrr = db[i].mxx;
                          db[i].mtt = db[i].myy;
                          db[i].mff = db[i].mzz;
                          db[i].mrt = db[i].mxy;
                          db[i].mrf = db[i].mxz;
                          db[i].mtf = db[i].myz;

                          db[i].mxx =  db[i].mtt;
                          db[i].myy =  db[i].mff;
                          db[i].mzz =  db[i].mrr;
                          db[i].mxy = -db[i].mtf;
                          db[i].mxz =  db[i].mrt;
                          db[i].myz = -db[i].mrf;
                }
                else if( strcmp( MT_coordinate_system, "r" ) == 0 || strcmp( MT_coordinate_system, "remap" ) == 0 )
                {
                /**** save to tmp space, remap cartesian ***/
                        mzz = db[i].mxx;
                        mxx = db[i].myy;
                        myy = db[i].mzz;
                        mxz = db[i].mxy;
                        myz = db[i].mxz;
                        mxy = db[i].myz;

                /*** write tmp space back ****/
                        db[i].mxx = mxx;
                        db[i].myy = myy;
                        db[i].mzz = mzz;
                        db[i].mxy = mxy;
                        db[i].mxz = mxz;
                        db[i].myz = myz;

                /*** convert cartesian to spherical ***/
                        db[i].mrr =   db[i].mzz;
                        db[i].mtt =   db[i].mxx;
                        db[i].mff =   db[i].myy;
                        db[i].mrt =   db[i].mxz;
                        db[i].mrf =  -db[i].myz;
                        db[i].mtf =  -db[i].mxy;
                }
                else
                {
                        fprintf(stderr, "%s: sys must be either c, s, or r got %s\n",
                                progname, MT_coordinate_system );
                        exit(0);
                }

                if(verbose)
                {
                  fprintf( stderr, "%s: %s: %s: MT_coordinate_system=%s read Mxx=%g Myy=%g Mzz=%g Mxy=%g Mxz=%g Myz=%g\n",
                       progname, __FILE__, __func__, MT_coordinate_system,
                       db[i].mxx, db[i].myy, db[i].mzz, db[i].mxy, db[i].mxz, db[i].myz );

                  fprintf( stderr, "%s: %s: %s: MT_coordinate_system=%s read mrr=%g mtt=%g mff=%g mrt=%g mrf=%g mtf=%g\n",
                       progname, __FILE__, __func__, MT_coordinate_system,
                       db[i].mrr, db[i].mtt, db[i].mff, db[i].mrt, db[i].mrf, db[i].mtf);
                }

                db[i].lune_lat =  -999;
                db[i].lune_lon =  -999;
		i++;
        }
        *nrec = i;
        fclose(fp);
        fprintf(stderr, "%s: %s: loaded nrec=%d\n", __FILE__, __func__, *nrec );
        return (DBrow *)db;
}

DBrow *load_dbrows( char *dbfilename, char *MT_coordinate_system, int *nrec, int verbose )
{
	DBrow *db;
	FILE *fp;
	char rec[512];
	int i = 0;
	long longdumb;
	float mxx, myy, mzz, mxy, mxz, myz;

	void truncate_strings( char *str, int n );

	if( (fp = fopen( dbfilename, "r" )) == NULL )
	{
		fprintf( stderr, "%s: %s: cannot open file %s\n",
			__FILE__, __func__, dbfilename );
		exit(-1);
	}

	db = calloc( 1, sizeof(DBrow) );
	while( fgets( rec, 512, fp ) != NULL )
	{
		db = realloc( db, (i+1)*sizeof(DBrow) );

		/*             1  2   3 4  5  6  7  8  9  10 11 12 13 14 15 16 17  18  19  20  21  22 23 24  25  26 27     */
		sscanf( rec, "%ld %s %f %e %d %d %d %f %e %e %e %e %e %e %s %f %ld %ld %ld %ld %ld %ld %s %f %f %s %[^\n]s",
			&(db[i].oridin), /* 1 */
			db[i].mt_type,
			&(db[i].momag),
			&(db[i].m0),
			&(db[i].pdc),   /* 5 */
			&(db[i].pclvd),
			&(db[i].piso),
			&(db[i].vred),
			&(db[i].mxx),
			&(db[i].myy),   /* 10 */
			&(db[i].mzz),
			&(db[i].mxy),
			&(db[i].mxz),
			&(db[i].myz),
			db[i].auth,       /* 15 */
			&(db[i].mtdepth), 
			&(db[i].evid),
			&longdumb,
			&(db[i].oridout),
			&(db[i].mtid),    /* 20 */
			&(db[i].fpid),
			&(db[i].magid),
			db[i].ot_str,
			&(db[i].lat),
			&(db[i].lon),       /* 25 */
			db[i].lddate_str,
			db[i].Region );

		truncate_strings( db[i].Region, 64 );

/****
To explain, I first converted from cartesian to spherical coordinates
	mzz mxx myy mxz  myz  mxy
	mrr mtt mff mrt -mrf -mtf
Later switched back in the wrong order…
	mrr mtt mff mrt -mrf -mtf
	mxx myy mzz mxy  mxz  myz
****/
		if( strcmp( MT_coordinate_system, "c" ) == 0 || strcmp( MT_coordinate_system, "cartesian" ) == 0 )
		{
		  db[i].mrr =   db[i].mzz;
                  db[i].mtt =   db[i].mxx;
                  db[i].mff =   db[i].myy;
                  db[i].mrt =   db[i].mxz;
                  db[i].mrf =  -db[i].myz;
                  db[i].mtf =  -db[i].mxy;
		}
		else if( strcmp( MT_coordinate_system, "s" ) == 0 || strcmp( MT_coordinate_system, "spherical" ) == 0 )
		{
			  db[i].mrr = db[i].mxx;
                          db[i].mtt = db[i].myy;
                          db[i].mff = db[i].mzz;
                          db[i].mrt = db[i].mxy;
                          db[i].mrf = db[i].mxz;
                          db[i].mtf = db[i].myz;

			  db[i].mxx =  db[i].mtt;
			  db[i].myy =  db[i].mff;
			  db[i].mzz =  db[i].mrr;
			  db[i].mxy = -db[i].mtf;
			  db[i].mxz =  db[i].mrt;
			  db[i].myz = -db[i].mrf;
		}
		else if( strcmp( MT_coordinate_system, "r" ) == 0 || strcmp( MT_coordinate_system, "remap" ) == 0 )
		{
				
		/**** save to tmp space, remap cartesian ***/
		  	mzz = db[i].mxx;
			mxx = db[i].myy;
			myy = db[i].mzz;
			mxz = db[i].mxy;
			myz = db[i].mxz;
			mxy = db[i].myz;

		/*** write tmp space back ****/
			db[i].mxx = mxx;
			db[i].myy = myy;
			db[i].mzz = mzz;
			db[i].mxy = mxy;
			db[i].mxz = mxz;
			db[i].myz = myz;

		/*** convert cartesian to spherical ***/
			db[i].mrr =   db[i].mzz;
                       	db[i].mtt =   db[i].mxx;
                       	db[i].mff =   db[i].myy;
                       	db[i].mrt =   db[i].mxz;
                        db[i].mrf =  -db[i].myz;
                       	db[i].mtf =  -db[i].mxy;
		}
		else
		{
                       	fprintf(stderr, "%s: sys must be either c or s got %s\n",
                               	progname, MT_coordinate_system );
                       	exit(0);
		}

		if(verbose)
                {       
                           fprintf( stderr, "%s: %s: %s: MT_coordinate_system=%s read Mxx=%g Myy=%g Mzz=%g Mxy=%g Mxz=%g Myz=%g\n",
                               progname, __FILE__, __func__, MT_coordinate_system,
                               db[i].mxx, db[i].myy, db[i].mzz, db[i].mxy, db[i].mxz, db[i].myz );

                           fprintf( stderr, "%s: %s: %s: MT_coordinate_system=%s read mrr=%g mtt=%g mff=%g mrt=%g mrf=%g mtf=%g\n",
                               progname, __FILE__, __func__, MT_coordinate_system,
                               db[i].mrr, db[i].mtt, db[i].mff, db[i].mrt, db[i].mrf, db[i].mtf);
		}

		db[i].lune_lat =  -999;
		db[i].lune_lon =  -999;

		i++;		
	}
	*nrec = i;
	fclose(fp);
	fprintf(stderr, "%s: %s: loaded nrec=%d\n", __FILE__, __func__, *nrec );
	return (DBrow *)db;
}

void fix_ot_str( char *str )
{
	int i, n=25;
	for( i = 0; i < n; i++ )if( str[i] == ',' ) str[i]='T';
}

void truncate_strings( char *str, int n )
{
	int i;
	for( i = 0; i < n; i++ )
	{
		/* fprintf( stdout, "[%d (%c)] ", i, str[i] ); */
		if( str[i] == ',' ) str[i]='-';

		if( str[i] == '\0' ) break;
		if( ( str[i+1] == ' ' ) && ( str[i] == ' ' ) ) 
		{
			str[i] = '\0';
			break;
		}
	}
	/* fprintf( stdout, "\n" ); */
}

/*** ***/
void write_dbrows_csv( DBrow *db, int nrec )
{
	int i;
	FILE *fp;
	fp = fopen("ML_mt.csv", "w" );
	fprintf( fp, "lune_lon,lune_lat,pvr,pdc,pclvd,piso,Mo,mxx,myy,mzz,mxy,mxz,myz,eig1,eig2,eig3,origin_time,lat,lon,momag,depth,evid,auth,region,label\n" );
	for( i = 0; i < nrec; i++ )
	{
	  fprintf( fp, 
"%+.4f,%+.4f,%.2f,%d,%d,%d,%+.3e,%+.3e,%+.3e,%+.3e,%+.3e,%+.3e,%+.3e,%+.3e,%+.3e,%+.3e,%s,%+.3f,%+.3f,%.2f,%.2f,%ld,%s,%s,%d\n",
		db[i].lune_lon,
		db[i].lune_lat,
		db[i].vred, 
		db[i].pdc,
		db[i].pclvd,
		db[i].piso,
		db[i].Mtotal,
		db[i].mxx/db[i].Mtotal,
                db[i].myy/db[i].Mtotal,
                db[i].mzz/db[i].Mtotal,
                db[i].mxy/db[i].Mtotal,
                db[i].mxz/db[i].Mtotal,
                db[i].myz/db[i].Mtotal,
		db[i].eval_full[1]/db[i].Mtotal,
		db[i].eval_full[2]/db[i].Mtotal,
		db[i].eval_full[3]/db[i].Mtotal,
		db[i].ot_str,
		db[i].lat,
                db[i].lon,
		db[i].momag,
		db[i].mtdepth,
		db[i].evid,
		db[i].auth,
		db[i].Region,
		db[i].label );
	}
	fclose(fp);
}

void write_gmt( DBrow *db, int nrec )
{
	int i;
	for( i = 0; i < nrec; i++ )
	{
		fprintf( stdout,
"%+12.4f %+12.4f %6.2f %8s %3.2f %4.2e %3d %3d %3d %6.1f %+6.2e %+6.2e %+6.2e %+6.2e %+6.2e %+6.2e %10s %5.1f   # evid=%10ld oridin=%10ld %24s lat=%+6.3f lon=%+8.3f lune_lat=%+7.3f lune_lon=%+8.3f (%s) (%s)\n",
		db[i].lune_lon,
		db[i].lune_lat,
		db[i].vred, /* db[i].oridin, %-12ld */
                db[i].mt_type,
                db[i].momag,
                db[i].m0,
                db[i].pdc,
                db[i].pclvd,
                db[i].piso,
                db[i].vred,
                db[i].mxx,
                db[i].myy,
                db[i].mzz,
                db[i].mxy,
                db[i].mxz,
                db[i].myz,
                db[i].auth,
                db[i].mtdepth,
                db[i].evid,
                db[i].oridin,
                db[i].ot_str,
                db[i].lat,
                db[i].lon,
                db[i].lune_lat,
                db[i].lune_lon,
                db[i].lddate_str,
                db[i].Region );
	}
}

void write_MTLoader( DBrow *db, int nrec )
{
	int i;
	FILE *fp;
	float Mw;

	fp = fopen("MTLoader.fixed.txt", "w");
	
fprintf( fp, "orid     mt_type momag  m0         pdc pclvd piso var_red mxx        myy        mzz        mxy         mxz        myz         author     depth\n" );
	
	for( i = 0; i < nrec; i++ )
	{
		Mw = ( log10(db[i].Mtotal) + 7 )/1.5 - 10.73;

		if( Mw < 9 ) 
		{
fprintf( fp,
"%-12ld %8s %3.2f %4.2e %3d %3d %3d %6.1f %+6.2e %+6.2e %+6.2e %+6.2e %+6.2e %+6.2e %10s %5.1f   # evid=%10ld oridin=%10ld %24s lat=%+6.3f lon=%+8.3f lune_lat=%+7.3f lune_lon=%+8.3f (%s) (%s)\n",
		db[i].oridin,
		db[i].mt_type,
		db[i].momag,
		db[i].m0,
		db[i].pdc,
		db[i].pclvd,
		db[i].piso,
		db[i].vred,
                db[i].mxx,
                db[i].myy,
                db[i].mzz,
                db[i].mxy,
                db[i].mxz,
                db[i].myz,
                db[i].auth,
                db[i].mtdepth,
                db[i].evid,
                db[i].oridin,
                db[i].ot_str,
                db[i].lat,
                db[i].lon,
                db[i].lune_lat,
                db[i].lune_lon,
                db[i].lddate_str,
                db[i].Region );

		}
	}
	fclose(fp);
}

void write_results( DBrow *db, int nrec )
{
	float Mw;
	int i;
	for( i = 0; i < nrec; i++ )
	{
		/* Mo in dyne cm */
		Mw = ( log10(db[i].Mtotal) + 7 )/1.5 - 10.73;

		fprintf( stdout,
"%5s mw=%3.2f %3.2f m0=%4.2e %5.2e dc=%3d %3.0f clvd=%3d %3.0f iso=%3d %3.0f lat=%6.2f lon=%6.2f xx=%5.2e yy=%5.2e zz=%5.2e xy=%5.2e xz=%5.2e yz=%5.2e evid=%10ld oridin=%10ld ot=(%s) k=%5.2f T=%5.2f f=%.2f (%s) (%22s)\n",
			db[i].mt_type,
                        db[i].momag,
			Mw,
                        db[i].m0,
			db[i].Mtotal,
                        db[i].pdc,
			db[i].pdc2,
                        db[i].pclvd,
			db[i].pclvd2,
                        db[i].piso,
			db[i].piso2,
			db[i].lune_lat,
                        db[i].lune_lon,
                        db[i].mxx,
                        db[i].myy,
                        db[i].mzz,
                        db[i].mxy,
                        db[i].mxz,
                        db[i].myz,
                        db[i].evid,
                        db[i].oridin,
			db[i].ot_str,
			db[i].k,
			db[i].T,
			db[i].f_factor,
                        db[i].lddate_str,
                        db[i].Region );
	}
}

void write_dbrows( DBrow *db, int nrec )
{
	int i;
	for( i = 0; i < nrec; i++ )
	{

fprintf( stdout,
"%5s mw=%3.2f m0=%4.2e dc=%3d clvd=%3d iso=%3d vr=%6.2f xx=%5.2e yy=%5.2e zz=%5.2e xy=%5.2e xz=%5.2e yz=%5.2e auth=%10s z=%5.1f evid=%10ld oridin=%10ld oridout=%10ld mtid=%10ld fpid=%10ld magid=%10ld %8s lat=%.1g lon=%.1g lune=%.1f/%.1f (%s) (%-32.32s)\n",
                        db[i].mt_type,
                        db[i].momag,
                        db[i].m0,
                        db[i].pdc,
                        db[i].pclvd,
                        db[i].piso,
                        db[i].vred,
                        db[i].mxx,
                        db[i].myy,
                        db[i].mzz,
                        db[i].mxy,
                        db[i].mxz,
			db[i].myz,
                        db[i].auth,
                        db[i].mtdepth,
                        db[i].evid,
                        db[i].oridin,
                        db[i].oridout,
                        db[i].mtid,
                        db[i].fpid,
                        db[i].magid,
                        db[i].ot_str,
                        db[i].lat,
                        db[i].lon,
			db[i].lune_lat,
			db[i].lune_lon, 
                        db[i].lddate_str,
                        db[i].Region );
	}
}



void mt2eig_version2( MomentTensor Ma, float *eig, int verbose )
{
	float **z1, **z2;
	float *eval1, *eval2, *evp, *evd;
	int i, j, nx=3;
	void tred2( float **, int, float *, float * );
	void tqli(float d[], float e[], int n, float **z);

/*** allocate memory ***/
        eval1 = vector( 0, 4 );
        eval2 = vector( 0, 4 );
        evp = vector( 0, 4 );
        evd = vector( 0, 4 );
        z1 = matrix( 0, 4, 0, 4 );
        z2 = matrix( 0, 4, 0, 4 );

/*** initalize arrays ***/
        for( i=1; i<=nx; i++ )
        {
                for( j=1; j<=nx; j++)
                {
                        z1[i][j] = Ma.mt[i][j];
                        z2[i][j] = z1[i][j];
                }
        }
/*** do the eigenvalue/eigenvector calculations ***/

        if(verbose)
        {
          fprintf( stderr, "%s: %s: %s: z1[][] = \n %g %g %g \n %g %g %g \n %g %g %g\n",
                progname, __FILE__, __func__,
                z1[1][1], z1[1][2], z1[1][3],
                z1[2][1], z1[2][2], z1[2][3],
                z1[3][1], z1[3][2], z1[3][3] );
        }

        tred2( z2, nx, eval1, eval2 );

/*** dont print out "z2" or "eval2", its just scratch space ***/

        if(verbose)
          fprintf( stderr, "%s: %s: %s: eval1: %g %g %g\n",
                progname, __FILE__, __func__, eval1[1], eval1[2], eval1[3] );

        tqli( eval1, eval2, nx, z2 );

/*** write out ***/

        if(verbose)
        {
                for( i=1; i<=nx; i++ )
                {
                        /*** eigenvalues ***/
                        fprintf( stderr, "%s: %s: %s: eval1=%6.2e z2=",
                                progname, __FILE__, __func__, eval1[i]);

                        /*** eigenvectors ***/
                        for( j=1; j<=nx; j++)
                        {
                                fprintf( stderr, "%6.2f ", z2[i][j] );
                        }
                        fprintf(stderr, "\n" );
                }
        }
	
/*** set return values ***/
	eig[0] = eval1[1];
	eig[1] = eval1[2];
	eig[2] = eval1[3];
}


void eig2iso_version2( MomentTensor Ma, float *eig, 
                float *Mtotal, float *Mdev, float *Miso, float *Mclvd, float *Mdc,
                float *k, float *T, float *f_factor, 
                float *pdc2, float *pclvd2, float *pdev2, float *piso2, 
                int verbose )
{
	float sign[4];
	int index[4];
	float etmp[4];

	int nx=3, i, j;
	float iso;
	static float zero = 1.0E-09;
	float epsilon;
	float eval_full[4]; /* 1,2,3 */
	float eval_dev[4];
	float eval_devsort[4];

	float  *tmp_abs_dev_eigen;
        float  *tmp_dev_eigen;
        float   max_abs_dev_eigen   = 0;
        float   min_abs_dev_eigen   = 0;
	float   min_dev_eigen_corsp = 0;

	void floatsort( float *, int );
	float fsign( float );
	void indexx( int, float *, int * );
	void flt_sort_asc(int npts, float *fvec_a, float *fvec_b);

/*** compute isotropic moment ***/
	iso = ( eig[0] + eig[1] + eig[2] )/3;
	*Miso = fabs(iso);

/** set the eigenvalues of the full moment tensor ***/
	for( i = 1; i <= nx; i++ )
		eval_full[i] = eig[i-1];

/*** calculate deviatoric moment ***/
/**** DEVIATORIC TENSOR - remove iso ****/
	for( i=1; i<=3; i++ )
		eval_dev[i] = eval_full[i] - iso;

/********************************************************/
/*** k value source type Hudson et al 1989            ***/
/********************************************************/

/* ************ New way, as per Gordon Kraft (JIRN, 19 Feb 2014) ********* */
/* Sort the absolute values of the deviatoric components and sort          */
/* a copy of the deviatoric components in the corresponding manner.        */
/* Because the first arguments is ignored, pass the address of the         */
/* first index, not the zeroth index.   The logic below is from lam2Tk.m   */
/* lam2Tk.m from Carl Tape (12/2012 matlab code):                          */
/* http://compearth.googlecode.com/svn/trunk/momenttensor/matlab/lam2Tk.m  */
/* Using the lam2Tk.m convenction, T should then be defined as -2*epsilon. */
/*****************************************************************************/

	tmp_abs_dev_eigen   = (float *)calloc( 4, sizeof(float) );
	tmp_dev_eigen       = (float *)calloc( 4, sizeof(float) );
	max_abs_dev_eigen   = 0.0f;
	min_abs_dev_eigen   = 0.0f;
	min_dev_eigen_corsp = 0.0f;

	tmp_dev_eigen[0]     = 0.0f;
        tmp_dev_eigen[1]     = eval_dev[1];
        tmp_dev_eigen[2]     = eval_dev[2];
        tmp_dev_eigen[3]     = eval_dev[3];

	tmp_abs_dev_eigen[0] = 0.0f;
        tmp_abs_dev_eigen[1] = fabs( eval_dev[1] );
        tmp_abs_dev_eigen[2] = fabs( eval_dev[2] );
        tmp_abs_dev_eigen[3] = fabs( eval_dev[3] );

/**********************************************************************************/
/* Sort the abs(eigen) values, but keep a one-to-one match with the signed values */
/* Need to pass index one, since flt_sort_asc uses standard C array indexing.     */
/**********************************************************************************/

        flt_sort_asc( 3, &(tmp_abs_dev_eigen[1]), &(tmp_dev_eigen[1]) );

	min_abs_dev_eigen   = tmp_abs_dev_eigen[1];		/* min( abs(eigen[]) )                   */
        min_dev_eigen_corsp = tmp_dev_eigen[1];     		/* Signed version of the min(abs(eigen)) */
        max_abs_dev_eigen   = tmp_abs_dev_eigen[3]; 		/* max( abs(eigen[]) )                   */

	*Mdev = max_abs_dev_eigen;

/* Get the values for k and epsilon */

        *k = iso / ( fabs(iso) + max_abs_dev_eigen );

        if( fabs(min_dev_eigen_corsp) <= zero || fabs(max_abs_dev_eigen) <= zero )
        {
                epsilon  = 0;
        }
        else
        {
                epsilon  = (-1.0) * min_dev_eigen_corsp / max_abs_dev_eigen;
        }

/********************************************************/
/*** epsilon value source type ( Hudson et al 1989 ) ***/
/********************************************************/

/**********************************************************************************/
/* Gene, It was much cleaner to keep the computation of the k and epsilon   */
/* together, since they new way reuse the same sorted arrays.  For viewing  */
/* the new of computing epsilon, go up ~30 lines of code JIRN, 19 Feb 2014  */
/**********************************************************************************/

        *Mclvd  = 2.0 * fabs( epsilon ) * *Mdev;
        *Mdc    = *Mdev - *Mclvd;
        *Mtotal = *Mdev + *Miso;

/***********************************************************/
/*** tectonic release f-factor (Toksoz and Kehrner 1972) ***/
/***********************************************************/

	*f_factor = 0;
        if( fabs( *Miso ) != 0 )
        {
          *f_factor = fabs( *Mdev ) / fabs( *Miso );
        }
        else
        {
          *f_factor = fabs( *Mdev ) / zero;
        }
	*T  = 2*(epsilon);
	*pdc2   = 100*( *Mdc   / *Mtotal );
	*pdev2  = 100*( *Mdev  / *Mtotal );
        *piso2  = 100*( *Miso  / *Mtotal );
        *pclvd2 = 100*( *Mclvd / *Mtotal );
}

float fsign( float x )
{
        float y = 1;
        if( x >= 0 ) return y;
        if( x <  0 ) return -y;
        return 1;
}


void flt_sort_asc(int npts, float *fvec_a, float *fvec_b)
{
        /* Begin declaration(s) of function prototype(s) */
        void flt_swap(float *, float *);

        /* Begin declaration(s) of local variable(s)     */
        int i;
        int  j  = 1;

        /* Sort the a and b arrays; b is sorted with respect of a */
        if( npts < 1 )  {
                fprintf(stderr, "*** %s(): Numer of points less than one.  Exiting!\n", __func__);
                return;
        }

        for(i=0; i<npts; i++)  {
                j = i+1;
                while( j < npts )  {
                        if( fvec_a[i] > fvec_a[j])  {

                                flt_swap( &(fvec_a[i]), &(fvec_a[j]));

                                if( fvec_b != NULL)  {
                                        flt_swap( &(fvec_b[i]), &(fvec_b[j]));
                                }
                        }
                        j++;
                }
        }

        return;
}


void flt_swap(float *a, float *b)
{
        /* Begin declaration(s) of local variable(s) */
        float tempflt = 0.0f;

        /* Verify that the arguments are not NULL.  Otherwise, swap the values */
        if( (a == NULL) || (b == NULL) )  {
                fprintf(stderr, "*** %s(): WARNING:  One or both of the arguments is/are NULL.\n", __func__);
        }
        else  {
                tempflt = *a;
                *a      = *b;
                *b      = tempflt;
        }

        return;
}
