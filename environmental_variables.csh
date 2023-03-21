#!/bin/tcsh

set workingdir = `pwd`
echo ${workingdir}
# echo $PATH

### MacOS
###
setenv PATH    ${PATH}:"${workingdir}/bin"
#setenv MANPATH ${MANPATH}:"${workingdir}/man"

### Linux
###
# append_path PATH ${workingdir}/bin
# append_path MANPATH ${workingdir}/man

echo $PATH 
#echo $MANPATH

# setenv MTINV_PATH          /Users/ichinose1/Work/mtinv.v3.0.6
# setenv MT_DATABASE_FILE    /Users/ichinose1/Work/mtinv.v3.0.6/data/mt.db
#
setenv MTINV_PATH        "${workingdir}"
setenv MT_DATABASE_FILE  "${workingdir}/data/mt.db"

echo $MTINV_PATH
echo $MT_DATABASE_FILE

### Global Topography, E-topo5 world grids 
### ask me if you want these.
###
setenv MTINV_GMT_GRID_FILE /Users/ichinose1/Work/topogmt/etopo5.grd
setenv MTINV_GMT_INT_FILE  /Users/ichinose1/Work/topogmt/etopo5.int
setenv MTINV_GMT_CPT_FILE  /Users/ichinose1/Work/topogmt/etopo5.cpt

echo $MTINV_GMT_GRID_FILE
echo $MTINV_GMT_INT_FILE
echo $MTINV_GMT_CPT_FILE
