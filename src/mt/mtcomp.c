#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/nrutil.h"     /** numerical recipes **/
#include "../include/mt.h"         /** global datatype and structure declarations **/

#define VTOL 1.0E-07
#define TTOL 1.0

char progname[128];

/*** this program inputs

     Total Moment
     DC     - S/D/R   %Mo
     CLVD   - Vertical, Horizontal, %Mo
     ISO    - %Mo

     output moment tensor
            moment tensor decomposition
            GMT plot deviatoric only
****/

int main( int ac, char **av )
{
/*** local variables and function prototypes ***/

	Tensor ISO_MT, CLVD_MT, DC_MT, MT, PBT_MT, SDR_EIG_MT;
	float MoTotal = 0, str = 0, dip = 0, rak = 0, pdc = 0;
	float ppl = 0, paz = 0, tpl = 0, taz = 0;
	float eigP = 0, eigB = 0, eigT = 0; /*** in same units of Moment ***/
	char clvd_type[6];
	float pclvd = 0, piso = 0;
	float yaw_z = 0, pitch_y = 0, roll_x = 0; /*** 3D Tensor rotation - CLVD MT ***/

	void computeFullMT( Tensor isoM, Tensor dcM, Tensor clvdM, Tensor SDR_EIG_MT, Tensor PBT_MT, Tensor *M );

	void compute_PBT( Tensor *M, float taz, float tpl, float paz, float ppl, float eigP, float eigB, float eigT );

	void compute_SDR_Eig( Tensor *M, float str, float dip, float rak, float eigP, float eigB, float eigT );

	void create_CLVD( float Mo, float pclvd, char *clvdType, Tensor *M,
		float yaw_z, float pitch_y, float roll_x );

	void create_DC( float Mo, float pdc, float str, float dip, float rak, Tensor *M );
	void create_ISO( float Mo, float piso, Tensor *M );
	void writeTensor( Tensor *M, char *label );
	void initTensor( Tensor *T );

/*** variables and functions from mtinv toolkit libglib.a and mtinv_subs.o ***/

	MomentTensor Ma, Mn;
        float x_vector[7];
        int mtdegfree = 6;
        int verbose   = 0;
        Solution *sol;
        int iz;

/* see include/mt.h Mo = math.pow( 10.0, 1.5*(Mw+10.73) ) = 1.2445146117713818e+16; Reference Mw = 0.0 */

	void set_moment_tensor( MomentTensor *Ma, float *x, int idf, int verbose );
	void normalize_moment_tensor( MomentTensor *Ma, MomentTensor *Mn, float Mtotal, int verbose );
	void mt2eig( MomentTensor Mn, Solution *sol, int iz, int verbose );
	void eig2iso( Solution *sol, int iz, int verbose );
	void Eig2MajorDC( Solution *sol, int iz, int verbose );
	void Eig2MinorDC( Solution *sol, int iz, int verbose );
	void eig2lune_4mtinv( Solution *sol, int iz, int verbose );
	
/*** misc ***/
	int setpar(int,char **), mstpar(), getpar();	
	void endpar();

/**********************/
/**** begin program ***/
/**********************/
	strcpy( progname, av[0] );

	setpar(ac,av);
	mstpar("Mo", "f", &MoTotal);

/*** getpar input Double Couple ***/
	getpar("str", "f", &str );
	getpar("dip", "f", &dip );
	getpar("rak", "f", &rak );

	mstpar("pdc", "f", &pdc );

/*** getpar input Compensated Linear Vector Dipole ***/
/*** horizontal1, horizontal2 or vertical - negative ( -h1, -h2, -v ) ***/
/*** horizontal1, horizontal2 or vertical - positive ( +h1, +h2, +v ) ***/
	sprintf( clvd_type, "+v" );
	getpar("clvd_type", "s", clvd_type );
	mstpar("pclvd", "f", &pclvd );

/*** 3D Tensor rotation - CLVD MT, angles in degrees ***/
	getpar( "yaw_z",   "f", &yaw_z );
	getpar( "pitch_y", "f", &pitch_y );
	getpar( "roll_x",  "f", &roll_x );

/*** eigenvalues and PBT or SDR moment tensor ***/
	getpar( "eigP", "f", &eigP );
	getpar( "eigB", "f", &eigB );
	getpar( "eigT", "f", &eigT );
	getpar( "taz", "f", &taz );
	getpar( "tpl", "f", &tpl );
	getpar( "paz", "f", &paz );
	getpar( "ppl", "f", &ppl );
	
/*** GetPar Input Isotropic ***/
	mstpar("piso", "f", &piso );

	getpar("verbose","b",&verbose);
	endpar();

	fprintf( stderr, "%s: Mo=%e str=%g dip=%g rak=%g pdc=%g clvd_type=%s pclvd=%g piso=%g (%g,%g,%g)\n",
		progname, MoTotal, str, dip, rak, pdc, clvd_type, pclvd, piso,
		yaw_z, pitch_y, roll_x );

/*** create the moment tensor ***/
	
	create_ISO( MoTotal, piso, &ISO_MT );
	writeTensor( &ISO_MT, "isotropic" );

	create_DC( MoTotal, pdc, str, dip, rak, &DC_MT );
	writeTensor( &DC_MT, "double-couple" );

	create_CLVD( MoTotal, pclvd, clvd_type, &CLVD_MT,
		yaw_z, pitch_y, roll_x );
	writeTensor( &CLVD_MT, "compensated linear vector dipole" );

	initTensor( &SDR_EIG_MT );
/*
	compute_SDR_Eig( &SDR_EIG_MT, str, dip, rak, eigP, eigB, eigT );
	writeTensor( &SDR_EIG_MT, "eigenvalues and SDR MT" );
*/
	initTensor( &PBT_MT );
	compute_PBT( &PBT_MT, taz, tpl, paz, ppl, eigP, eigB, eigT );
	writeTensor( &PBT_MT, "eigenvalues and PT-axes MT" );

	initTensor( &MT );
	computeFullMT( ISO_MT, DC_MT, CLVD_MT, SDR_EIG_MT, PBT_MT, &MT );
	writeTensor( &MT, "full - moment tensor" );


/*** set and normalize the moment tensor ***/

	x_vector[1] = MT.xx;
	x_vector[2] = MT.yy;
	x_vector[3] = MT.xy;
	x_vector[4] = MT.xz;
	x_vector[5] = MT.yz;
	x_vector[6] = MT.zz;

	fprintf( stdout, "%+6.2e %+6.2e %+6.2e %+6.2e %+6.2e %+6.2e\n",
		MT.xx, MT.xy, MT.xz, MT.yy, MT.yz, MT.zz );

	mtdegfree = 6;
	set_moment_tensor( &Ma, x_vector, mtdegfree, verbose );
	
	sol = (Solution *)malloc(2*sizeof(Solution));
	iz = 0;
	sol[iz].mt_type  = FULL_MT;

/*** decompose the moment tensor ***/

	mt2eig( Ma, sol, iz, verbose );
	eig2iso( sol, iz, verbose );

/*** normalize all moment tensor and eigenvalues ***/

	normalize_moment_tensor( &Ma, &Mn, sol[iz].Mtotal, verbose );
	sol[iz].dmoment  = Ma.moment;
        sol[iz].mw       = Ma.Mw;
        sol[iz].exponent = Ma.expon;
        sol[iz].abcassa  = Ma.abcassa;

	sol[iz].FullMT.eig[1].val /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.eig[2].val /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.eig[3].val /= pow(10.0, sol[iz].exponent );

        sol[iz].Dev.eig[1].val /= pow(10.0, sol[iz].exponent );
        sol[iz].Dev.eig[2].val /= pow(10.0, sol[iz].exponent );
        sol[iz].Dev.eig[3].val /= pow(10.0, sol[iz].exponent );

        sol[iz].FullMT.T.ev /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.B.ev /= pow(10.0, sol[iz].exponent );
        sol[iz].FullMT.P.ev /= pow(10.0, sol[iz].exponent );

        sol[iz].mrr = Mn.rr;
        sol[iz].mtt = Mn.tt;
        sol[iz].mff = Mn.ff;
        sol[iz].mrt = Mn.rt;
        sol[iz].mrf = Mn.rf;
        sol[iz].mtf = Mn.tf;

        sol[iz].mxx = Mn.xx;
        sol[iz].mxy = Mn.xy;
        sol[iz].mxz = Mn.xz;
        sol[iz].myy = Mn.yy;
        sol[iz].myz = Mn.yz;
        sol[iz].mzz = Mn.zz;

	Eig2MajorDC( sol, iz, verbose );
	Eig2MinorDC( sol, iz, verbose );
	eig2lune_4mtinv( sol, iz, verbose );
}

void initTensor( Tensor *T ) { 
	T->xx = 0; T->xy = 0; T->xz = 0;
	T->yx = 0; T->yy = 0; T->yz = 0;
	T->zx = 0; T->zy = 0; T->zz = 0;
}

void symetricTensor( Tensor *T ) { T->yx = T->xy; T->zx = T->xz; T->zy = T->yz; }

void compute_PBT( Tensor *M, float taz, float tpl, float paz, float ppl, float eigP, float eigB, float eigT )
{
	float pplr, pazr, tazr, tplr;
	float baz=0, bpl=0;
	Vector3 t, b, p;
	Tensor Evec, Eval, EvecT, Mtmp;
	char label[256];

	void productTensor( Tensor *M1, Tensor *M2, Tensor *M3 );
        void copyTensor( Tensor *Mi, Tensor *Mo );
        void transposeTensor( Tensor *Mi, Tensor *Mo );
        void initTensor( Tensor *T );
        void writeTensor( Tensor *M, char *label );
        void writeVector3( Vector3 *v, char *label );
	
/*** convert deg to radians ***/
	tazr = taz * (M_PI/180);
	tplr = tpl * (M_PI/180);
	pazr = paz * (M_PI/180);
	pplr = ppl * (M_PI/180);
	
/*** P-axes Eigenvectors from azimuth and plunge to directional cosines ***/
	p.x = cos( pazr ) * cos( pplr );
        p.y = sin( pazr ) * cos( pplr );
        p.z = sin( pplr );
        sprintf( label, "PBT EIG: P-axes vector from eigP=%+7.2e paz=%3.0f ppl=%3.0f", eigP, paz, ppl );
        writeVector3( &p, label );

/*** T-axes Eigenvectors from azimuth and plunge to directional cosines ***/
	t.x = cos( tazr ) * cos( tplr );
	t.y = sin( tazr ) * cos( tplr );
	t.z = sin( tplr );
	sprintf( label, "PBT EIG: T-axes vector from eigT=%+7.2e taz=%3.0f tpl=%3.0f", eigT, taz, tpl );
	writeVector3( &t, label );

/*** B axes from t curl p ***/
	b.x =  ( t.y * p.z - t.z * p.y );
	b.y =  ( t.z * p.x - t.x * p.z );
	b.z =  ( t.x * p.y - t.y * p.x );
	if( b.x == 0 )
		baz = 0;
	else
		baz = (180/M_PI) * atan( b.y / b.x );
	bpl = (180/M_PI) * asin( b.z );
	sprintf( label, "PBT EIG: B-axes vector from eigB=%+7.2e baz=%3.0f bpl=%3.0f", eigB, baz, bpl );
        writeVector3( &b, label );

/*** eigenvalues ***/
	initTensor( &Eval );
	Eval.xx = eigT;
	Eval.yy = eigB;
	Eval.zz = eigP;
	writeTensor( &Eval, "PBT EIG: Eigenvalues T  B  P" );

/*** eigenvectors ***/
	initTensor( &Evec );
	Evec.xx = t.x;  Evec.xy = b.x; Evec.xz = p.x;
	Evec.yx = t.y;  Evec.yy = b.y; Evec.yz = p.y;
	Evec.zx = t.z;  Evec.zy = b.z; Evec.zz = p.z;
	writeTensor( &Evec, "PBT EIG: Eigenvectors  T   B   P" );

/*** generate moment tensors  Evec * Eval * EvecT = MT ***/
/*** For orthonganol symmetric square matrix EvecT = inv(Evec)  ***/

	initTensor( &EvecT );
	transposeTensor( &Evec, &EvecT );
	/* writeTensor( &EvecT, "PBT Eigenvectors Transpose" ); */

	productTensor( &Evec, &Eval, &Mtmp );
	productTensor( &Mtmp, &EvecT, M );
	/* writeTensor( M, "Moment Tensor - from PBT eigval and eigvec" ); */
}

void compute_SDR_Eig( Tensor *M, float str, float dip, float rak, float eigP, float eigB, float eigT )
{
	float strr,dipr,rakr;
	Vector3 u, v, t, b, p;
	Tensor Evec, Eval, EvecT, Mtmp;

	void productTensor( Tensor *M1, Tensor *M2, Tensor *M3 );
        void copyTensor( Tensor *Mi, Tensor *Mo );
        void transposeTensor( Tensor *Mi, Tensor *Mo );
	void initTensor( Tensor *T );
	void writeTensor( Tensor *M, char *label );
	void writeVector3( Vector3 *v, char *label );

/*** convert deg to radians ***/
	strr = str * (M_PI / 180);	
	dipr = dip * (M_PI / 180);
	rakr = rak * (M_PI / 180);

/*** slip vector ***/
	u.x =  cos(rakr)*cos(strr) + cos(dipr)*sin(rakr)*sin(strr);
	u.y =  cos(rakr)*sin(strr) - cos(dipr)*sin(rakr)*cos(strr);
	u.z = -sin(dipr) * sin(rakr);
	writeVector3( &u, "U vector" );

/*** fault normal vector ***/
	v.x = -sin(dipr)*sin(strr);
	v.y = sin(dipr)*cos(strr);
	v.z = -cos(dipr);
	writeVector3( &v, "V vector" );

/*** t axes vector ***/
	t.x = (1/M_SQRT2) * ( v.x + u.x );
	t.y = (1/M_SQRT2) * ( v.y + u.y );
	t.z = (1/M_SQRT2) * ( v.z + u.z );
	writeVector3( &t, "T vector" );

/**** v curl u = b ***/
	b.x = ( v.y * u.z - u.y * v.z );
	b.y = ( v.z * u.x - u.z * v.x );
	b.z = ( v.x * u.y - u.x * v.y );
	writeVector3( &t, "B vector = v curl u" );

/*** p axes vector ***/
	p.x = (1/M_SQRT2) * ( v.x - u.x );
	p.y = (1/M_SQRT2) * ( v.y - u.y );
	p.z = (1/M_SQRT2) * ( v.z - u.z );
	writeVector3( &t, "P vector" );

/**** t curl p = b ***/
/*** test, produces same result as v curl u =b ***/
/***
	b.x = ( t.y * p.z - p.y * t.z );
        b.y = ( t.z * p.x - p.z * t.x );
        b.z = ( t.x * p.y - p.x * t.y );
        writeVector3( &t, "B vector = T curl P" );
***/

/*** eigenvalues ***/
	initTensor( &Eval );
/*	float str = 0, dip = 45, rak = 90; */
/*
	Eval.xx = +1.0E+23;
        Eval.yy =  0.0E+00;
        Eval.zz = -1.0E+23;
*/
	Eval.xx = eigT;
	Eval.yy = eigB;
	Eval.zz = eigP;
	writeTensor( &Eval, "SDR Eigenvalues T  B  P" );

/*** eigenvectors ***/
	initTensor( &Evec );
	Evec.xx = t.x;  Evec.xy = b.x; Evec.xz = p.x;
        Evec.yx = t.y;  Evec.yy = b.y; Evec.yz = p.y;
        Evec.zx = t.z;  Evec.zy = b.z; Evec.zz = p.z;
	writeTensor( &Evec, "SDR Eigenvectors  T   B   P" );

/*** generate moment tensors  Evec * Eval * EvecT = MT ***/
/*** For orthonganol symmetric square matrix EvecT = inv(Evec)  ***/

	initTensor( &EvecT );
	transposeTensor( &Evec, &EvecT );
	/* writeTensor( &EvecT, "Eigenvectors Transpose" ); */

	productTensor( &Evec, &Eval, &Mtmp );
	productTensor( &Mtmp, &EvecT, M );
	/* writeTensor( M, "Moment Tensor - from eigval and eigvec" ); */
}

void computeFullMT( Tensor isoM, Tensor dcM, Tensor clvdM, Tensor sdr, Tensor pbt, Tensor *M )
{
	M->xx = isoM.xx + dcM.xx + clvdM.xx + sdr.xx + pbt.xx;
	M->yy = isoM.yy + dcM.yy + clvdM.yy + sdr.yy + pbt.yy;
	M->zz = isoM.zz + dcM.zz + clvdM.zz + sdr.zz + pbt.zz;
	M->xy = isoM.xy + dcM.xy + clvdM.xy + sdr.xy + pbt.xy;
	M->xz = isoM.xz + dcM.xz + clvdM.xz + sdr.xz + pbt.xz;
	M->yz = isoM.yz + dcM.yz + clvdM.yz + sdr.yz + pbt.yz;

/*** forces symmetry ***/
	M->yx = M->xy;
	M->zx = M->xz;
	M->zy = M->yz;
}

/*** horizontal1, horizontal2 or vertical - negative ( -h1, -h2, -v ) ***/
/*** -2  0  0      +1  0  0      +1  0  0
      0 +1  0       0 -2  0       0 +1  0
      0  0 +1       0  0 +1       0  0 -2  ***/
/*** horizontal1, horizontal2 or vertical - positive ( +h1, +h2, +v ) ***/
/***
     +2  0  0      -1  0  0      -1  0  0
      0 -1  0       0 +2  0       0 -1  0
      0  0 -1       0  0 -1       0  0 +2  ***/

void create_CLVD( float Mo, float pclvd, char *clvdType, Tensor *MCLVD,
		float yaw_z, float pitch_y, float roll_x )
{
	float factor = 0.5;
	Tensor M;

	void rotateTensor( Tensor *Mi, Tensor *Mo, float yaw_z, float pitch_y, float roll_x );
	void writeTensor( Tensor *M, char *label );
	void copyTensor( Tensor *Mi, Tensor *Mo );

/*** start ***/

	if( strcmp(clvdType,"+v") == 0 )
	{
		M.xx = -1.0 * Mo * factor * pclvd;
		M.yy = -1.0 * Mo * factor * pclvd;
		M.zz = +2.0 * Mo * factor * pclvd;
	}
	else if( strcmp(clvdType,"+h2") == 0 )	
	{
		M.xx = -1.0 * Mo * factor * pclvd;
                M.yy = +2.0 * Mo * factor * pclvd;
                M.zz = -1.0 * Mo * factor * pclvd;
	}
	else if( strcmp(clvdType,"+h1") == 0 )
	{
		M.xx = +2.0 * Mo * factor * pclvd;
                M.yy = -1.0 * Mo * factor * pclvd;
                M.zz = -1.0 * Mo * factor * pclvd;
	}
	else if( strcmp(clvdType,"-v") == 0 )
        {
		M.xx = +1.0 * Mo * factor * pclvd;
                M.yy = +1.0 * Mo * factor * pclvd;
                M.zz = -2.0 * Mo * factor * pclvd;
        }
        else if( strcmp(clvdType,"-h2") == 0 )  
        {
		M.xx = +1.0 * Mo * factor * pclvd;
                M.yy = -2.0 * Mo * factor * pclvd;
                M.zz = +1.0 * Mo * factor * pclvd;
        }
        else if( strcmp(clvdType,"-h1") == 0 )
        {
		M.xx = -2.0 * Mo * factor * pclvd;
                M.yy = +1.0 * Mo * factor * pclvd;
                M.zz = +1.0 * Mo * factor * pclvd;
        }

	M.xy = 0;
        M.xz = 0;
        M.yz = 0;
        M.yx = 0;
        M.zx = 0;
        M.zy = 0;
	
/*
	M.xx = Mo;
	M.yy = Mo;
	M.zz = Mo;
*/
	writeTensor( &M, "basis compensated linear vector dipole" );

	/*** debug ***/
	/* copyTensor( &M, MCLVD ); */

	rotateTensor( &M, MCLVD, yaw_z, pitch_y, roll_x );

	writeTensor( MCLVD, "rotated compensated linear vector dipole" ); 
}

void create_DC( float Mo, float pdc, float str, float dip, float rak, Tensor *M )
{
	float d2r, strr, dipr, rakr, tol = 1.0E-07;
	void scaleTensor( Tensor *M, float scale );
	void writeTensor( Tensor *M, char *label );
	void floorTensor( Tensor *M, float tol );

	d2r   = M_PI / 180.0;
	strr  = str * d2r;
	dipr  = dip * d2r;
	rakr  = rak * d2r;

	M->xx = -(sin(dipr)*cos(rakr)*sin(2*strr)+sin(2*dipr)*sin(rakr)*sin(strr)*sin(strr));
        M->yy =  (sin(dipr)*cos(rakr)*sin(2*strr)-sin(2*dipr)*sin(rakr)*cos(strr)*cos(strr));
        M->zz =  (sin(2*dipr)*sin(rakr));
        M->xy =  (sin(dipr)*cos(rakr)*cos(2*strr)+0.5*sin(2*dipr)*sin(rakr)*sin(2*strr));
        M->xz = -(cos(dipr)*cos(rakr)*cos(strr)+cos(2*dipr)*sin(rakr)*sin(strr));
        M->yz = -(cos(dipr)*cos(rakr)*sin(strr)-cos(2*dipr)*sin(rakr)*cos(strr));
        M->yx = M->xy;
        M->zx = M->xz;
        M->zy = M->yz;

	/* writeTensor( M, "DC before scaling" ); */
	floorTensor( M, tol );
	scaleTensor( M, (Mo*pdc) );
}

void floorTensor( Tensor *M, float tol )
{
	if( fabs(M->xx) < tol ) M->xx = 0;
	if( fabs(M->xy) < tol ) M->xy = 0;
	if( fabs(M->xz) < tol ) M->xz = 0;
	if( fabs(M->yy) < tol ) M->yy = 0;
	if( fabs(M->yz) < tol ) M->yz = 0;
	if( fabs(M->zz) < tol ) M->zz = 0;
	M->yx = M->xy;
	M->zx = M->xz;
	M->zy = M->yz;
}

void floorVector( Vector3 *v, float tol )
{
	if( fabs(v->x) < tol ) v->x = 0;
	if( fabs(v->y) < tol ) v->y = 0;
	if( fabs(v->z) < tol ) v->z = 0;
}

void scaleTensor( Tensor *M, float scale )
{
	M->xx *= scale;
	M->yy *= scale;
        M->zz *= scale;
        M->xy *= scale;
        M->xz *= scale;
        M->yz *= scale;
        M->yx = M->xy;
        M->zx = M->xz;
        M->zy = M->yz;
}

void create_ISO( float Mo, float piso, Tensor *M )
{
	/* float factor = 0.333333333333333333; */
	float factor = 1.0;

	M->xx = Mo * factor * piso;
	M->yy = Mo * factor * piso;
	M->zz = Mo * factor * piso;
	M->xy = 0;
	M->xz = 0;
	M->yz = 0;
	M->yx = 0;
	M->zx = 0;
	M->zy = 0;
}

void writeVector3( Vector3 *v, char *label ) 
{
	void floorVector( Vector3 *v, float tol );
	floorVector( v, VTOL );
	fprintf( stderr, "%s: x=%+6.2e y=%+6.2e z=%+6.2e\n", label, v->x, v->y, v->z );
}

void writeTensor( Tensor *M, char *label )
{
	void floorTensor( Tensor *T, float tol );

	/* floorTensor( M, TTOL ); */

	fprintf( stderr, "%s\n", label );
	fprintf( stderr, "\t %+6.2e %+6.2e %+6.2e\n", M->xx, M->xy, M->xz );
	fprintf( stderr, "\t %+6.2e %+6.2e %+6.2e\n", M->yx, M->yy, M->yz );
	fprintf( stderr, "\t %+6.2e %+6.2e %+6.2e\n", M->zx, M->zy, M->zz );
	fprintf( stderr, "\n" );
}

void transposeTensor( Tensor *Mi, Tensor *Mo )
{
	Tensor T;
/*
	void copyTensor( Tensor *Mi, Tensor *Mo );
	copyTensor( Mi, &T );
*/
/*** Transpose second order rank tensor Aij = Aji ***/
	Mo->xx = Mi->xx;
	Mo->yy = Mi->yy;
	Mo->zz = Mi->zz;

	Mo->xy = Mi->yx;
	Mo->xz = Mi->zx;
	Mo->yx = Mi->xy;
	Mo->yz = Mi->zy;
	Mo->zx = Mi->xz;
	Mo->zy = Mi->yz;
}

void rotateTensor( Tensor *Mi, Tensor *Mo, float yaw_z, float pitch_y, float roll_x )
{
	Tensor Rx, Ry, Rz, R;
	float alpha, beta, gamma, d2r;
	Tensor RxT, RyT, RzT;

	void productTensor( Tensor *M1, Tensor *M2, Tensor *M3 );
	void copyTensor( Tensor *Mi, Tensor *Mo );
	void transposeTensor( Tensor *Mi, Tensor *Mo );

/***********************/
/*** rotation matrix ***/
/***********************/

	d2r = M_PI / 180.0;
	alpha = yaw_z   * d2r;
	beta  = pitch_y * d2r;
	gamma = roll_x  * d2r;

	Rz.xx = cos( alpha ); Rz.xy = -sin( alpha ); Rz.xz = 0;
	Rz.yx = sin( alpha ); Rz.yy = cos( alpha );  Rz.xz = 0;
	Rz.zx = 0;            Rz.zy = 0;             Rz.zz = 1;

	Ry.xx = cos( beta );  Ry.xy = 0;             Ry.xz = sin( beta );
	Ry.yx = 0;            Ry.yy = 1;             Ry.yz = 0;
	Ry.zx = -sin( beta ); Ry.zy = 0;             Ry.zz = cos( beta );

	Rx.xx = 1;            Rx.xy = 0;             Rx.xz = 0;
	Rx.yx = 0;            Rx.yy = cos( gamma );  Rx.yz = -sin( gamma );
	Rx.zx = 0;            Rx.zy = sin( gamma );  Rx.zz = cos( gamma );

	productTensor( &Rx, Mi, Mo );
	copyTensor( Mo, Mi );

	productTensor( &Ry, Mi, Mo );
	copyTensor( Mo, Mi );

	productTensor( &Rz, Mi, Mo );
	
/*** identity matrix R * R' = I ***/
/***
	productTensor( &Rx, Mi, Mo );
	transposeTensor( &Rx, &RxT );
	productTensor( Mo, &RxT, Mi );

	productTensor( &Ry, Mi, Mo );
	transposeTensor( &Ry, &RyT );
        productTensor( Mo, &RyT, Mi );

	productTensor( &Rz, Mi, Mo );
	transposeTensor( &Rz, &RzT );
        productTensor( Mo, &RzT, Mi );

	copyTensor( Mi, Mo );
***/

}

void productTensor( Tensor *M1, Tensor *M2, Tensor *M3 )
{
	M3->xx = M1->xx*M2->xx + M1->xy*M2->yx + M1->xz*M2->zx;
	M3->xy = M1->xx*M2->xy + M1->xy*M2->yy + M1->xz*M2->zy;
	M3->xz = M1->xx*M2->xz + M1->xy*M2->yz + M1->xz*M2->zz;
	M3->yx = M1->yx*M2->xx + M1->yy*M2->yx + M1->yz*M2->zx;
        M3->yy = M1->yx*M2->xy + M1->yy*M2->yy + M1->yz*M2->zy;
        M3->yz = M1->yx*M2->xz + M1->yy*M2->yz + M1->yz*M2->zz;
	M3->zx = M1->zx*M2->xx + M1->zy*M2->yx + M1->zz*M2->zx;
        M3->zy = M1->zx*M2->xy + M1->zy*M2->yy + M1->zz*M2->zy;
        M3->zz = M1->zx*M2->xz + M1->zy*M2->yz + M1->zz*M2->zz;
}

void copyTensor( Tensor *Mi, Tensor *Mo )
{
	Mo->xx = Mi->xx;
	Mo->xy = Mi->xy;
	Mo->xz = Mi->xz;
	Mo->yx = Mi->yx;
	Mo->yy = Mi->yy;
	Mo->yz = Mi->yz;
	Mo->zx = Mi->zx;
	Mo->zy = Mi->zy;
	Mo->zz = Mi->zz;
}
