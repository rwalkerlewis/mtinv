#!/bin/csh

# This is to get the ZSS, RSS, TSS 

mtdecomp2 Mxx=1 Myy=-1 Mzz=0 Mxy=1 Mxz=0 Myz=0 Mo=1e20 verbose
mtdecomp str=158 dip=90 rak=0 Mo=1.0E+20 pdc=1 clvd_type=None pclvd=0 piso=0 verbose

