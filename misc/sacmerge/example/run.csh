#!/bin/csh

### Demostration of sacmerge
### merges these 3 sac files into 1 single trace, fills gaps with zeros
# 2007.185.00.00.00.0000.IU.KBL..BHZ.Q.SAC
# 2007.185.00.35.50.0000.IU.KBL..BHZ.Q.SAC
# 2007.185.14.28.43.0250.IU.KBL..BHZ.R.SAC

##clean up outputs
##
/bin/rm -f 2007.185.00.00.00.0000.IU.KBL..BHZ.D.SAC

### merge these three sac files, read input SAC filenames from command line, wildcards accepted
### NOTE: I've decimated the data to 0.2 sec/sample just so the tar file of this example can be emailed under 6MB.
###
# sacmerge 2007*KBL*BHZ*SAC
../../../bin/sacmerge 2007*KBL*BHZ*SAC

### output written to sac file, gaps are filled with zeros and overlap is filled using later file
###

### Plot result, assumes you have SAC installed in your executable PATH
###
cat >! sac.plot << EOF
bd x
qdp off
read *.SAC
p1
pause
color on inc on
p2
pause
quit
EOF

sac sac.plot
