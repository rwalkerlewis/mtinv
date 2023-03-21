#!/bin/csh

#set DATABASE=/Users/ichinos1/Work/mtinv.v3.0.6/data/mt.db
set DATABASE=${MT_DATABASE_FILE}

if( $#argv == 1 ) then
        set MTID=$argv[1]
	echo DATABASE=${DATABASE}
else if ( $#argv == 2 ) then
        set MTID=$argv[1]
        set DATABASE=$argv[2]
	echo DATABASE=${DATABASE}
else
        echo "$0 needs 1 argument"
        echo " Usage : $0 MTID (required) DATABASE (optional,default=${DATABASE})"
        exit
endif

sqlite3 ${DATABASE} << EOF
delete from MT_EARTHMODEL_STAGE where mtvmodid in 
  ( select mtvmodid from MT_WAVEFORM_SEGMENT_STAGE where mtdataid in ( select mtdataid from MT_DATA_STAGE where mtid = ${MTID} ) );
delete from MT_FILTER_STAGE where mtfilterid in 
  ( select mtfilterid from MT_WAVEFORM_SEGMENT_STAGE where mtdataid in ( select mtdataid from MT_DATA_STAGE where mtid = ${MTID} ) );
delete from MT_WAVEFORM_SEGMENT_STAGE where mtdataid in ( select mtdataid from MT_DATA_STAGE where mtid = ${MTID} );
delete from MT_DATA_STAGE where mtid = ${MTID};
delete from FOCAL_PLANE_STAGE where orid = (select orid from MOMENT_STAGE where mtid = ${MTID} );
delete from MT_ORIGIN_STAGE where orid = ( select orid from MOMENT_STAGE where mtid = ${MTID} );
delete from MOMENT_STAGE where mtid = ${MTID};
.quit
EOF

##
## G. Ichinose 2019March11
## How do we orphan MOMENT_STAGE rows with no corresponding MT_ORIGIN_STAGE rows by orid link?
##
