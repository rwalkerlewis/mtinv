#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/mt.h"         /** global datatype and structure declarations **/

extern char progname[128];

void wrtgrn2sac( Greens *g, int ista, char *wavetype )
{
	/* int ngreen = 10 */
	int ngreen = 22;
        Sac_Header sp;
        FILE *fp;
        char filename[256];

	void sac_minmax( float *x, int n, float *max, float *min, float *mean );

 /*** old ***/
 /* char cmp[][3] = { "rss", "rds", "rdd", "rep", "zss", "zds", "zdd", "zep", "tss", "tds" }; */

/**** new ***/
/*** mtinv.v4.0.0 -  Green's function now include 3C rotational motions w1,w2,w3 ***/
	char cmp[][5] = { "rss", "rds", "rdd", "rep", 
                          "zss", "zds", "zdd", "zep", 
                          "tss", "tds", 
                          "w1ss", "w1ds", "w1dd", "w1ex", 
                          "w2ss", "w2ds", "w2dd", "w2ex", 
                          "w3ss", "w3ds", "w3dd", "w3ex" };
        int i;
	int verbose = 0;

/*** misc_tools/ampshift.c ***/

	int remove_offset = 0;

        void  ampshift( float *x, int n, int verbose );

/******************************************************/
/*** for shallow sources, this shifts the amplitude ***/
/*** of the first sample to zero, and mitigates the ***/
/*** artifact of the step offset sometimes from Grn ***/
/*** functions from shallow sources, before filter  ***/
/******************************************************/

	fprintf( stdout, "%s: %s: %s: ista=%d sta=%s wavetype=%s\n",
                        progname, __FILE__, __func__, ista, g->stnm, wavetype );

	if( remove_offset )
	{
		ampshift( &(g->g.rss[0]), g->nt, verbose );
		ampshift( &(g->g.rds[0]), g->nt, verbose );
		ampshift( &(g->g.rdd[0]), g->nt, verbose );
		ampshift( &(g->g.rep[0]), g->nt, verbose );

		ampshift( &(g->g.zss[0]), g->nt, verbose );
		ampshift( &(g->g.zds[0]), g->nt, verbose );
		ampshift( &(g->g.zdd[0]), g->nt, verbose );
		ampshift( &(g->g.zep[0]), g->nt, verbose );

		ampshift( &(g->g.tss[0]), g->nt, verbose );
		ampshift( &(g->g.tds[0]), g->nt, verbose );

		if( strcmp( wavetype, "Rotational" ) == 0 )
		{
		 ampshift( &(g->g.w1ss[0]), g->nt, verbose );
		 ampshift( &(g->g.w1ds[0]), g->nt, verbose );
		 ampshift( &(g->g.w1dd[0]), g->nt, verbose );
		 ampshift( &(g->g.w1ex[0]), g->nt, verbose );

		 ampshift( &(g->g.w2ss[0]), g->nt, verbose );
		 ampshift( &(g->g.w2ds[0]), g->nt, verbose );
		 ampshift( &(g->g.w2dd[0]), g->nt, verbose );
		 ampshift( &(g->g.w2ex[0]), g->nt, verbose );

		 ampshift( &(g->g.w3ss[0]), g->nt, verbose );
		 ampshift( &(g->g.w3ds[0]), g->nt, verbose );
		 ampshift( &(g->g.w3dd[0]), g->nt, verbose );
		 ampshift( &(g->g.w3ex[0]), g->nt, verbose );
		}
	}

	if( strcmp( wavetype, "Surf/Pnl" ) == 0 ) ngreen = 10;
	if( strcmp( wavetype, "Rotational" ) == 0 ) ngreen = 22;

        for( i=0; i<ngreen; i++ )
        {

                sprintf( filename, "%s.%03d.%s.%g.grns", g->filename, ista, cmp[i], g->evdp );

                sp              = sac_null;
                sp.b            = g->t0;
                sp.delta        = g->dt;
                sp.e            = g->t0 + (g->nt * g->dt);
                sp.stla         = g->stla;
                sp.stlo         = g->stlo;
                sp.evla         = g->evla;
                sp.evlo         = g->evlo;
                sp.evdp         = g->evdp;

                sp.nvhdr = 6;
                sp.norid = 0;
                sp.nevid = 0;
                sp.iftype = ITIME; /* data type: IRLIM spec re im | IAMPH amp ph | IXY general x,y  */
                sp.idep   = IUNKN; /* not disp vel acc or volts */
                sp.iztype = IB;    /* types: IUNKN,IB,IDAY,IO,IA,ITn  */
                sp.ievtyp = IUNKN; /* type of event IQUAKE - earthquake */
                sp.npts = g->nt;

                sp.leven  = TRUE;  /* is data evenly sampled  */
                sp.lpspol = FALSE; /* LLL sets it */
                sp.lcalda = TRUE;  /* should az,baz,gcarc,dist be calculated?  */
                sp.lhdr5  = FALSE; /* LLL sets it */

		sprintf( sp.knetwk, "%s", g->net );
		sprintf( sp.kstnm,  "%s", g->stnm );
		sprintf( sp.khole,  "%s", g->loc );
                sprintf( sp.kcmpnm, "%s", cmp[i] );

		if( strcmp( wavetype, "Surf/Pnl" ) == 0 )
		{
		  if( (fp = fopen(filename,"wb")) == NULL )
		  {
			printf("%s: cannot open output file %s\n", progname, filename );
			exit(-1);
		  }
		  if( i==0 ) sac_minmax( &(g->g.rss[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==1 ) sac_minmax( &(g->g.rds[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==2 ) sac_minmax( &(g->g.rdd[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==3 ) sac_minmax( &(g->g.rep[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==4 ) sac_minmax( &(g->g.zss[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==5 ) sac_minmax( &(g->g.zds[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==6 ) sac_minmax( &(g->g.zdd[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==7 ) sac_minmax( &(g->g.zep[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==8 ) sac_minmax( &(g->g.tss[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 
                  if( i==9 ) sac_minmax( &(g->g.tds[0]), g->nt, &(sp.depmax), &(sp.depmin), &(sp.depmen) ); 

                  fwrite( &sp, sizeof(Sac_Header), 1, fp );

                  if( i==0 ) fwrite( &(g->g.rss[0]), g->nt * sizeof(float), 1, fp );
                  if( i==1 ) fwrite( &(g->g.rds[0]), g->nt * sizeof(float), 1, fp );
                  if( i==2 ) fwrite( &(g->g.rdd[0]), g->nt * sizeof(float), 1, fp );
                  if( i==3 ) fwrite( &(g->g.rep[0]), g->nt * sizeof(float), 1, fp );
                  if( i==4 ) fwrite( &(g->g.zss[0]), g->nt * sizeof(float), 1, fp );
                  if( i==5 ) fwrite( &(g->g.zds[0]), g->nt * sizeof(float), 1, fp );
                  if( i==6 ) fwrite( &(g->g.zdd[0]), g->nt * sizeof(float), 1, fp );
                  if( i==7 ) fwrite( &(g->g.zep[0]), g->nt * sizeof(float), 1, fp );
                  if( i==8 ) fwrite( &(g->g.tss[0]), g->nt * sizeof(float), 1, fp );
                  if( i==9 ) fwrite( &(g->g.tds[0]), g->nt * sizeof(float), 1, fp ); 

		  fclose(fp);
		}
		else if( strcmp( wavetype, "Rotational" ) == 0 )
		{
		  if( (fp = fopen(filename,"wb")) == NULL )
		  {
			printf("%s: cannot open output file %s\n", progname, filename );
			exit(-1);
		  }
		  fwrite( &sp, sizeof(Sac_Header), 1, fp );

		  if( i==10 ) fwrite( &(g->g.w1ss[0]), g->nt * sizeof(float), 1, fp );
                  if( i==11 ) fwrite( &(g->g.w1ds[0]), g->nt * sizeof(float), 1, fp );
                  if( i==12 ) fwrite( &(g->g.w1dd[0]), g->nt * sizeof(float), 1, fp );
                  if( i==13 ) fwrite( &(g->g.w1ex[0]), g->nt * sizeof(float), 1, fp );

                  if( i==14 ) fwrite( &(g->g.w2ss[0]), g->nt * sizeof(float), 1, fp );
                  if( i==15 ) fwrite( &(g->g.w2ds[0]), g->nt * sizeof(float), 1, fp );
                  if( i==16 ) fwrite( &(g->g.w2dd[0]), g->nt * sizeof(float), 1, fp );
                  if( i==17 ) fwrite( &(g->g.w2ex[0]), g->nt * sizeof(float), 1, fp );

                  if( i==18 ) fwrite( &(g->g.w3ss[0]), g->nt * sizeof(float), 1, fp );
                  if( i==19 ) fwrite( &(g->g.w3ds[0]), g->nt * sizeof(float), 1, fp );
		  if( i==20 ) fwrite( &(g->g.w3dd[0]), g->nt * sizeof(float), 1, fp );
                  if( i==21 ) fwrite( &(g->g.w3ex[0]), g->nt * sizeof(float), 1, fp );

		  fclose(fp);

		} /**** if rotational wavetype ***/

        } /*** loop over i-th greens function ***/

} /*** end of wrtgrn2sac.c ***/


/*** duplicated in sacio/sacio.c ***/

/*** 
void sac_minmax( float *x, int n, float *max, float *min, float *mean )
{
        int i;
        float sum=0;
        *max = x[0];
        *min = x[0];
        *mean = 0;
        for( i=0; i<n; i++ )
        {
                if( x[i] > *max ) *max=x[i];
                if( x[i] < *min ) *min=x[i];
                sum = sum + x[i];
        }
        *mean = sum/n;
}

***/
