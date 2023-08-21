#!/bin/csh

rm -f x.sac x.sac.int xint.sac

# create random timeseries, low pass filter below nyquist before deciation using SAC
# 
# fnyq = 1/(2*0.01) = 50 Hz
# fnyq = 1/(2*0.20) = 2.5 Hz
#
sac << EOF
fg random 1 1 npts 2048 delta 0.01
lowpass butter npole 3 corner 2.5
write x.sac
EOF

## test the SAC interpolate command
##
sac << EOF
read x.sac
lh npts delta
interpolate delta 0.2
write append .int
quit
EOF

## test the wiggins interpolate 
##
wiggins f=x.sac nt=102 dt=0.2

### plot all the results 
###
cat >! sac.mac << EOF
read x.sac xint.sac x.sac.int
qdp off
color on inc on
lh npts delta
p1
pause
p2
pause
quit
EOF

sac sac.mac

# clean up
/bin/rm -f sac.mac x.sac xint.sac x.sac.int
