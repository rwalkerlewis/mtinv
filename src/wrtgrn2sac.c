#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

/*** mkdirp ***/
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "../include/mt.h"         /** global datatype and structure declarations **/

extern char progname[128];

void plotgrnlib_GMT5( Greens *g, int ista, char *wavetype, int make_output_dirs )
{
	int i;
	int gstart, gend;
	FILE *fp;
	char filename[256];
	char script_filename[12] = {"plotgrn.csh"};
	float yoffset;
	float begin, end;
	int first = 1, line_color_switch = 0;
	char line_color[][8] = { "red", "black" };
	char cmp[][5] = { "rss", "rds", "rdd", "rep",
                          "zss", "zds", "zdd", "zep",
                          "tss", "tds",
                          "w1ss", "w1ds", "w1dd", "w1ex",
                          "w2ss", "w2ds", "w2dd", "w2ex",
                          "w3ss", "w3ds", "w3dd", "w3ex" }; 
	char title_string[512];
	float peak;
	float find_abspeak_value_from_float_array( float *x, int nt );

/*** mkdirp ***/
        int mkdirp2( const char *, mode_t mode );
        char outDirectory[64];

/*************************/
/*** begin subroutine ***/
/*************************/

/*** option make_output_dirs assumes Greens is either already or will be written here, more neater to do this so files not all cluttered ***/

	if( make_output_dirs )
	{
		snprintf( outDirectory, 8, "%s", g->stnm );
		mkdirp2( outDirectory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		chdir( outDirectory );

		fprintf( stderr, "%s: %s: %s: make_output_dirs=%d Writting output to directory outDirectory=%s\n",
                        progname, __FILE__, __func__, make_output_dirs, outDirectory );
	}

/*** get the start and stop pointers from glib file depending on wavetype ***/
	if( strcmp( wavetype, "Surf/Pnl" ) == 0 ) 
	{
		gstart = 0; 
		gend = 10;
	}
        if( strcmp( wavetype, "Rotational" ) == 0 ) 
	{ 
		gstart = 10;
		gend = 22;
	}

/*** open the GMT cshell script ***/

	if( (fp = fopen(script_filename, "w")) == NULL )
	{
		fprintf( stderr, "%s: %s: %s: cannot open file %s\n", 
			progname, __FILE__, __func__, script_filename );
		exit(-1);
	}
	chmod( script_filename, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH );

	fprintf( fp, "#!/bin/csh\n" );
	fprintf( fp, "############################################################################\n" );
	fprintf( fp, "## This C-shell script was automatically generated                       ###\n" );
        fprintf( fp, "## and requires Generic Mapping Tools (GMT) http://gmt.soest.hawaii.edu/ ###\n" );
        fprintf( fp, "## SCRIPT ONLY WORKS WITH GMT Version 5.x.x                              ###\n" );
        fprintf( fp, "############################################################################\n" );
        fprintf( fp, "\n" );

	fprintf( fp, "gmt set PS_PAGE_ORIENTATION        portrait  \n" );
        fprintf( fp, "gmt set MAP_ANNOT_OFFSET_PRIMARY   2p         \n" );
        fprintf( fp, "gmt set MAP_ANNOT_OFFSET_SECONDARY 2p      \n" );
        fprintf( fp, "gmt set MAP_LABEL_OFFSET           0p       \n" );
        fprintf( fp, "gmt set PS_MEDIA                   letter  \n" );
        fprintf( fp, "gmt set FONT_ANNOT_PRIMARY          9p,Helvetica,black   \n" ); 
        fprintf( fp, "gmt set FONT_ANNOT_SECONDARY        9p,Helvetica,black  \n" );
        fprintf( fp, "gmt set FONT_LABEL                  9p,Palatino-Bold,black \n" );
        fprintf( fp, "gmt set FONT_LOGO                   9p,Helvetica,black     \n" );
        fprintf( fp, "gmt set FONT_TITLE                 12p,Times-Bold,black \n" );
	fprintf( fp, "\n" );

	/* fprintf( fp, "set R=\"-R0/100/-1/1\" \n" ); */
	begin = g->t0;
	end   = g->t0 + (g->nt * g->dt);

	fprintf( fp, "set R=\"-R%g/%g/-1/+1\" \n", begin, end );

	fprintf( fp, "set J=\"-JX6i/0.8i\" \n" );

	if( strcmp( wavetype, "Surf/Pnl" ) == 0 )
	{
	 fprintf( fp, "set  PS=%s.%s.%s.%s.%g.plotgrn.ps\n", g->net, g->stnm, g->loc, g->v.modfile, g->evdp );
	 fprintf( fp, "set JPG=%s.%s.%s.%s.%g.plotgrn.jpg\n", g->net, g->stnm, g->loc, g->v.modfile, g->evdp );
	}
	if( strcmp( wavetype, "Rotational" ) == 0 )
	{
	 fprintf( fp, "set  PS=%s.%s.%s.%s.%g.plotgrn.rot.ps\n", g->net, g->stnm, g->loc, g->v.modfile, g->evdp );
         fprintf( fp, "set JPG=%s.%s.%s.%s.%g.plotgrn.rot.jpg\n", g->net, g->stnm, g->loc, g->v.modfile, g->evdp );
	}

	fprintf( fp, "\n" );

	fprintf( fp, "gmt psbasemap $R $J -Bxf1a10+l\"seconds\" -Byf1a1 -BS -P -K >! ${PS}\n" );

        for( i=gstart; i<gend; i++ )
        {
		snprintf( filename, 256, "%s.%03d.%s.%g.grns", g->filename, ista, cmp[i], g->evdp );
		if( first == 1 ) {
			first = 0;
			yoffset = 0.0;
		}else{
			yoffset = 0.75;
		}
		if( (i%2) == 0 )
			line_color_switch = 0;
		else
			line_color_switch = 1;

		if( strcmp( wavetype, "Surf/Pnl" ) == 0 )
		{
		 if( i==0 ) peak = find_abspeak_value_from_float_array( &(g->g.rss[0]), g->nt );
                 if( i==1 ) peak = find_abspeak_value_from_float_array( &(g->g.rds[0]), g->nt );
                 if( i==2 ) peak = find_abspeak_value_from_float_array( &(g->g.rdd[0]), g->nt );
                 if( i==3 ) peak = find_abspeak_value_from_float_array( &(g->g.rep[0]), g->nt );
                 if( i==4 ) peak = find_abspeak_value_from_float_array( &(g->g.zss[0]), g->nt );
                 if( i==5 ) peak = find_abspeak_value_from_float_array( &(g->g.zds[0]), g->nt );
                 if( i==6 ) peak = find_abspeak_value_from_float_array( &(g->g.zdd[0]), g->nt );
                 if( i==7 ) peak = find_abspeak_value_from_float_array( &(g->g.zep[0]), g->nt );
                 if( i==8 ) peak = find_abspeak_value_from_float_array( &(g->g.tss[0]), g->nt );
                 if( i==9 ) peak = find_abspeak_value_from_float_array( &(g->g.tds[0]), g->nt );

		 fprintf( fp, "sac2xy norm < %s | gmt psxy $R $J -W0.5p,%-8s -Y%gi -K -O >> ${PS}\n",
                        filename, line_color[line_color_switch], yoffset );

                fprintf( fp, "echo %g 0 %s | gmt pstext $R $J -F+f11p,Helvetica-Bold,black+jMR -N -D0i/0.0i -O -K >> ${PS}\n",
                        begin, cmp[i] );

                fprintf( fp, "echo %g 0 %5.2e | gmt pstext $R $J -F+f11p,Helvetica-Bold,black+jML -N -D-0.50i/0.1i -O -K >> ${PS}\n",
                        end, peak );
		}

		if( strcmp( wavetype, "Rotational" ) == 0 )
		{
		  if( i==10 ) peak = find_abspeak_value_from_float_array( &(g->g.w1ss[0]), g->nt );
                  if( i==11 ) peak = find_abspeak_value_from_float_array( &(g->g.w1ds[0]), g->nt );
                  if( i==12 ) peak = find_abspeak_value_from_float_array( &(g->g.w1dd[0]), g->nt );
                  if( i==13 ) peak = find_abspeak_value_from_float_array( &(g->g.w1ex[0]), g->nt );

                  if( i==14 ) peak = find_abspeak_value_from_float_array( &(g->g.w2ss[0]), g->nt );
                  if( i==15 ) peak = find_abspeak_value_from_float_array( &(g->g.w2ds[0]), g->nt );
                  if( i==16 ) peak = find_abspeak_value_from_float_array( &(g->g.w2dd[0]), g->nt );
                  if( i==17 ) peak = find_abspeak_value_from_float_array( &(g->g.w2ex[0]), g->nt );

                  if( i==18 ) peak = find_abspeak_value_from_float_array( &(g->g.w3ss[0]), g->nt );
                  if( i==19 ) peak = find_abspeak_value_from_float_array( &(g->g.w3ds[0]), g->nt );
                  if( i==20 ) peak = find_abspeak_value_from_float_array( &(g->g.w3dd[0]), g->nt );
                  if( i==21 ) peak = find_abspeak_value_from_float_array( &(g->g.w3ex[0]), g->nt );

		fprintf( fp, "sac2xy norm < %s | gmt psxy $R $J -W0.5p,%-8s -Y%gi -K -O >> ${PS}\n",
			filename, line_color[line_color_switch], yoffset );

		fprintf( fp, "echo %g 0 %s | gmt pstext $R $J -F+f11p,Helvetica-Bold,black+jMR -N -D0i/0.0i -O -K >> ${PS}\n",
			begin, cmp[i] );

		fprintf( fp, "echo %g 0 %5.2e | gmt pstext $R $J -F+f11p,Helvetica-Bold,black+jML -N -D-0.50i/0.1i -O -K >> ${PS}\n",
			end, peak );
		}

		fprintf( fp, "\n" );
	}
	sprintf( title_string, "Greens functions %s net.sta.loc= %s.%s.%s model= %s depth= %g",
		wavetype, g->net, g->stnm, g->loc, g->v.modfile, g->evdp  );

	fprintf( fp, "gmt psbasemap $R $J -Bxf1a10+l\"seconds\" -Byf1a1 -BN+t\"%s\" -O  >> ${PS}\n", title_string );

	fprintf( fp, "gmt psconvert -Tj -E300 ${PS}\n" );

	fprintf( fp, "### cleanup and sys dep plotting\n" );
	fprintf( fp, "/bin/rm -f ${PS}\n" );
	fprintf( fp, "# open ${JPG}\n" );
	fclose(fp);

	system( "plotgrn.csh" );
	if( make_output_dirs ) chdir( ".." );
}	

void wrtgrn2sac( Greens *g, int ista, char *wavetype, int make_output_dirs )
{
	/* int ngreen = 10 */
	int ngreen = 22;
        Sac_Header sp;
        FILE *fp;
        char filename[256];

	void sac_minmax( float *x, int n, float *max, float *min, float *mean );
	
	/***prevents the writting of empty files ***/
	void  write_it( Sac_Header *sp, float *x, int nt, char *filename );

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

/*** mkdirp ***/
	int mkdirp2( const char *, mode_t mode );
	char outDirectory[64];

/******************************************************/
/*** for shallow sources, this shifts the amplitude ***/
/*** of the first sample to zero, and mitigates the ***/
/*** artifact of the step offset sometimes from Grn ***/
/*** functions from shallow sources, before filter  ***/
/******************************************************/

/*** attempt to make the directory if it doesn't exist already ***/

	if( make_output_dirs )
	{
		snprintf( outDirectory, 8, "%s", g->stnm );
		mkdirp2( outDirectory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		chdir( outDirectory );

		fprintf( stderr, "%s: %s: %s: make_output_dirs=%d Writting output to directory outDirectory=%s\n",
                        progname, __FILE__, __func__, make_output_dirs, outDirectory );
	}

	if(verbose) 
	{ 
		fprintf( stdout, "%s: %s: %s: ista=%d sta=%s wavetype=%s\n",
                        progname, __FILE__, __func__, ista, g->stnm, wavetype );
	}

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

			/*** 256 -> 256 1 3 1 3 1 */
                snprintf( filename, 256, 
			"%s.%03d.%s.%g.grns", 	
			g->filename, ista, cmp[i], g->evdp );

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
	
                  if( i==0 ) write_it( &sp, &(g->g.rss[0]), g->nt, filename );
                  if( i==1 ) write_it( &sp, &(g->g.rds[0]), g->nt, filename );
                  if( i==2 ) write_it( &sp, &(g->g.rdd[0]), g->nt, filename );
                  if( i==3 ) write_it( &sp, &(g->g.rep[0]), g->nt, filename );

                  if( i==4 ) write_it( &sp, &(g->g.zss[0]), g->nt, filename );
                  if( i==5 ) write_it( &sp, &(g->g.zds[0]), g->nt, filename );
                  if( i==6 ) write_it( &sp, &(g->g.zdd[0]), g->nt, filename );
                  if( i==7 ) write_it( &sp, &(g->g.zep[0]), g->nt, filename );

                  if( i==8 ) write_it( &sp, &(g->g.tss[0]), g->nt, filename );
                  if( i==9 ) write_it( &sp, &(g->g.tds[0]), g->nt, filename );
		}
		else if( strcmp( wavetype, "Rotational" ) == 0 )
		{
		  if( i==10 ) write_it( &sp, &(g->g.w1ss[0]), g->nt, filename );
                  if( i==11 ) write_it( &sp, &(g->g.w1ds[0]), g->nt, filename );
                  if( i==12 ) write_it( &sp, &(g->g.w1dd[0]), g->nt, filename );
                  if( i==13 ) write_it( &sp, &(g->g.w1ex[0]), g->nt, filename );

                  if( i==14 ) write_it( &sp, &(g->g.w2ss[0]), g->nt, filename );
                  if( i==15 ) write_it( &sp, &(g->g.w2ds[0]), g->nt, filename );
                  if( i==16 ) write_it( &sp, &(g->g.w2dd[0]), g->nt, filename );
                  if( i==17 ) write_it( &sp, &(g->g.w2ex[0]), g->nt, filename );

                  if( i==18 ) write_it( &sp, &(g->g.w3ss[0]), g->nt, filename );
                  if( i==19 ) write_it( &sp, &(g->g.w3ds[0]), g->nt, filename );
		  if( i==20 ) write_it( &sp, &(g->g.w3dd[0]), g->nt, filename );
                  if( i==21 ) write_it( &sp, &(g->g.w3ex[0]), g->nt, filename );

		} /**** if rotational wavetype ***/

        } /*** loop over i-th greens function ***/

	if( make_output_dirs ) chdir( ".." );

} /*** end of wrtgrn2sac.c ***/

void  write_it( Sac_Header *sp, float *x, int nt, char *filename )
{
	FILE *fp;
	void sac_minmax( float *x, int n, float *max, float *min, float *mean );

	if( (fp = fopen(filename,"wb")) == NULL )
	{
		fprintf( stderr, "%s: %s: %s: ERROR! cannot open output file %s\n", 
			progname, __FILE__, __func__, filename );
		exit(-1);
	}

	sac_minmax( x, nt, &(sp->depmax), &(sp->depmin), &(sp->depmen) );

	fwrite( sp, sizeof(Sac_Header), 1, fp );
	fwrite( x, nt * sizeof(float), 1, fp );
	fclose(fp);
}

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

/*** used here and in grnlib2sac ***/

float find_abspeak_value_from_float_array( float *x, int nt )
{
        int i; 
        float peak;
        peak = 1.0E-29;
        for( i = 0; i < nt; i++ )
        {
                if( fabs(x[i]) > peak )
                        peak = fabs(x[i]);
        }
        return peak;
}       
