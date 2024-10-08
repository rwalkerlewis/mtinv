#!/bin/csh
######################################################################
### Autogenerated C-shell script by setupMT.c                      ###
###  This C-shell script sets up mkgrnlin to create                ###
###  Green's functions from all the SAC data hdrs                  ###
######################################################################
### DataDir = (../Data)
### RespDir = (../Resp)
######################################################################

### don't forget to add ${MTINV_PATH}/bin to your executable path in the shell startup file ###

setenv MTINV_PATH /Users/ichinose1/Work/mtinv.v4.0.1

echo " MTINV_PATH = ${MTINV_PATH} " 

### 2021/11/18T02h53m04.00 36.9077 -90.543 New Madrid, MO

cat >! cus.par << EOF
velmod=cus
#zrange=3,3,33
zrange=1,1,25
evla=36.9077
evlo=-90.543
dt=0.15
nt=2048
fmax=0.4
t0=0.
redv=18.
damp=1.
kmax=20000
eps=0.0005
smin=0.0005
modeldb=${MTINV_PATH}/data/modeldb/
stadb=../Data/rdseed.stations
noverbose
nodump
EOF

cat >! mkgrnlib.par << EOF
### station-code network-code location-code mkgrnlib.parfile dt(sec/sample) ### comments
CGM3     NM     "00" cus.par     0.050 ### R=    90 Az=061 NM.CGM3.00   nchan=1 nseg=1 (HHZ)
PENM     NM     "00" cus.par     0.050 ### R=    96 Az=122 NM.PENM.00   nchan=1 nseg=1 (HHZ)
HENM     NM     "00" cus.par     0.050 ### R=    98 Az=102 NM.HENM.00   nchan=1 nseg=1 (HHZ)
GNAR     NM     "00" cus.par     0.050 ### R=   115 Az=156 NM.GNAR.00   nchan=1 nseg=1 (HHZ)
CCM      IU     "00" cus.par     0.070 ### R=   142 Az=334 IU.CCM.00    nchan=1 nseg=1 (BHZ)
CCM      IU     "10" cus.par     0.070 ### R=   142 Az=334 IU.CCM.10    nchan=2 nseg=2 (BHZ HHZ)
SIUC     NM     "00" cus.par     0.070 ### R=   148 Az=052 NM.SIUC.00   nchan=1 nseg=1 (HHZ)
SLM      NM     "00" cus.par     0.080 ### R=   194 Az=008 NM.SLM.00    nchan=2 nseg=2 (BHZ HHZ)
WVT      IU     "10" cus.par     0.110 ### R=   257 Az=109 IU.WVT.10    nchan=2 nseg=2 (BHZ HHZ)
WVT      IU     "00" cus.par     0.110 ### R=   257 Az=109 IU.WVT.00    nchan=2 nseg=2 (BHZ HHZ)
EOF

multithread_mkgrnlib \
     parfile=mkgrnlib.par \
     executable_pathname=${MTINV_PATH}/bin/mkgrnlib > multithread_mkgrnlib.out

makepar com="New Madrid, MO" \
    date="2021-11-18T02:53:04.00" \
    DataDir=../Data \
    RespDir=../Resp \
    gmt5 nooracle nomysql sqlite \
    maxsta=8 mindist=70 maxdist=800 \
    lf=0.075 hf=0.15 \
    minsnr=3.0 ctol=0.85 maxshift=10 realtime nolocal *.glib
###
### dump and plot the greens functions for single depth
###
# grnlib2sac glib=NM.CGM3."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=NM.PENM."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=NM.HENM."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=NM.GNAR."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=IU.CCM."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=IU.CCM."10".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=NM.SIUC."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=NM.SLM."00".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=IU.WVT."10".cus.glib dumpgrn plotgrn z=1
# grnlib2sac glib=IU.WVT."00".cus.glib dumpgrn plotgrn z=1
