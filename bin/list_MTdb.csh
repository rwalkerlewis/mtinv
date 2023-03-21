#!/bin/csh

if( $#argv == 0 ) then
	set DATABASE=/Users/ichinose/Work/mtinv.v3.0.6/data/mt.db
	echo "Default Database = ${DATABASE} "
else if ( $#argv == 1 ) then
	set DATABASE=$argv[1]
	echo "Database = ${DATABASE} "
else
        echo "$0 needs 1 argument"
	echo " Usage : $0 ${DATABASE}"
	echo "Default Database = ${DATABASE} "
        exit
endif

sqlite3 ${DATABASE} << EOF
.headers on
.mode col
.stats off
.timer off

.width -10
select sqlite_version();

.print "\n" 

.width  -7 -6 -19 -7 -8 -5 -5 -6 -5 -12 -32 -19
select
	a.orid,
	b.mtid,
        datetime( a.time, 'unixepoch' ) as OriginTime,
        a.lat,
        a.lon,
        a.depth,
	b.momag as Mw, 
	rtrim(b.pdc,0) as PDC,
	round(b.var_red,1) as VR,
	/* rtrim(c.strike1,'0') as S1, c.dip1 as D1, rtrim(c.rake1,'0') as R1, */
        a.algorithm,
	g.grname,
	a.lddate
from 
	MT_ORIGIN_STAGE a,
	MOMENT_STAGE b,
	FOCAL_PLANE_STAGE c,
	GEOREGION g
where
	( ( a.grn = g.grn ) or ( a.grn is NULL ) ) and
	a.orid = b.orid and
	b.fpid = c.fpid
order by time;

.print "\n"

.width 16
select count(*) as rows_counted from MT_ORIGIN_STAGE;

.print "\n"
.quit
EOF
