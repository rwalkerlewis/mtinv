#!/bin/csh
sqlite3 /Users/ichinose/Work/mtinv.v3.0.6/data/mt.db << EOF
.output dump.out
.mode column
.header off
select lon, lat, orid from MT_ORIGIN_STAGE;
.quit
EOF
