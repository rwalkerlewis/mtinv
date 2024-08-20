#!/bin/csh
#####################################################################
### this is an autogenerated script for mteig                     ###
### mteig compute %VR NSS on the lune following Ford et al., 2012 ###
#####################################################################
###
### Created by mtbestfit fitType=BEST VARIANCE-REDUCTION
### ot=-0 fsec=4 z=14 vred=86.5165 fit=86.5165 pdc=88 piso=10 pclvd=3
### fit_max=86.5165 fit_max_threshold=20 vred_diff=0 vred_diff_threshold=30
###

#set DEGFREE=6  ### mteig assumes degree-of-freedom is 6 for full-MT to compute %VR on lune

cat >! mtinv.par << EOF
#### REGION COMMENT ############################
CM New Madrid, MO
#### Date and Origin Time ######################
OT 2021-11-18T02:53:04.000
#### Forward Calculations ######################
##    stk    dip    rak   Mw  evlo  evla   Z ##########
EV   26   66   -169 4.09    -90.543    36.9077  14
#####################################################################################################
# sta net loc model  np pas lf  hf  nt  dt   tr  tt v/d  mulfac used(Y/N)  ts0  weight ###              #
CGM3     NM   00   cus  3  2 0.075 0.150  512 0.1 0 0 d 1 y +0.000 1 Surf/Pnl ## R=89.6 Az=61
PENM     NM   00   cus  3  2 0.075 0.150  512 0.1 0 0 d 1 y +0.000 1 Surf/Pnl ## R=96.1 Az=122
HENM     NM   00   cus  3  2 0.075 0.150  512 0.1 0 0 d 1 y +0.000 1 Surf/Pnl ## R=97.6 Az=102
GNAR     NM   00   cus  3  2 0.075 0.150  512 0.1 0 0 d 1 y +0.000 1 Surf/Pnl ## R=114.8 Az=156
CCM      IU   00   cus  3  2 0.075 0.150  512 0.14 0 0 d 1 y +0.000 1 Surf/Pnl ## R=141.8 Az=334
CCM      IU   10   cus  3  2 0.075 0.150  512 0.14 0 0 d 1 y +0.000 1 Surf/Pnl ## R=141.8 Az=334
SIUC     NM   00   cus  3  2 0.075 0.150  512 0.14 0 0 d 1 y +0.000 1 Surf/Pnl ## R=147.5 Az=52
SLM      NM   00   cus  3  2 0.075 0.150  512 0.16 0 0 d 1 y +0.000 1 Surf/Pnl ## R=193.9 Az=8
### NOTE! stas not used, are commented out and not loaded by mteig ###
### because there is no prediction only millions of forward calcs ###
EOF

### PROCESS GREENS FUNCTIONS ###
glib2inv par=mtinv.par noverbose parallel

### PROCESS DATA ###
sacdata2inv par=mtinv.par path=../Data respdir=../Resp noverbose nodumpsac parallel

### 
### changed nsim_eig=2000 to 500 to speed up runtime
### 
time mteig par=mtinv.par nthreads=8 \
                         nsim_eig=500 nsim_evec=4000 eigvec_fac=17000 \
                         Mo=1.670965e+22 fixz=14 \
                         color doplt seed=1 gmt5 parallel \
                         title="2021-11-18 New Madrid, MO" Add_user_Eig e0=+1.67096 e1=+0.144478 e2=-1.31764
