#!/bin/csh
sqlite3 /Users/ichinose/Work/mtinv.v3.0.6/data/mt.db << EOF
.read update.sql
.quit
EOF
