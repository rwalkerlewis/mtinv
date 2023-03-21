#!/bin/csh

#### example mixed
mtcomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.1 clvd_type=+v pclvd=0.2 piso=0.7

#### pure dc
#mtcomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=1.0 clvd_type=+v pclvd=0.0 piso=0.0

### pure clvd
#mtcomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.0 clvd_type=+v pclvd=1.0 piso=0.0

### pure ex
#mtcomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.0 clvd_type=+v pclvd=0.0 piso=1.0

### test sdr eig
# mtcomp verbose Mo=1.0E+23 str=0 dip=45 rak=90 pdc=0.0 clvd_type=+v pclvd=0.0 piso=0.0 eigP=-1.0E+23 eigB=0.0 eigT=+1.0E+23
# mtcomp verbose Mo=1.0E+23 str=45 dip=30 rak=90 pdc=0.0 clvd_type=+v pclvd=0.0 piso=1.0 eigP=-0.5E+23 eigB=-0.5E+23 eigT=+1.0E+23 >! out.1
# mtcomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.0 clvd_type=+v pclvd=0.0 piso=1.0 eigP=-0.5E+23 eigB=-0.5E+23 eigT=+1.0E+23 >! out.2
# mtcomp verbose Mo=1.0E+23 str=45 dip=60 rak=90 pdc=0.0 clvd_type=+v pclvd=0.0 piso=1.0 eigP=-0.5E+23 eigB=-0.5E+23 eigT=+1.0E+23 >! out.3

### test pbt eig

#mtcomp verbose Mo=0.0E+0 pdc=0.0 clvd_type=+v pclvd=0.0 piso=0.0 eigP=-1.0E+23 eigB=0.0 eigT=+1.0E+23 taz=180 tpl=90 paz=270 ppl=0
#mtcomp verbose Mo=0.0E+0 pdc=0.0 clvd_type=+v pclvd=0.0 piso=0.0 eigP=-1.0E+23 eigB=0.0 eigT=+1.0E+23 taz=45 tpl=0 paz=135 ppl=0
#mtcomp verbose Mo=0.0E+0 pdc=0.0 clvd_type=+v pclvd=0.0 piso=0.0 eigP=-1.0E+23 eigB=0.0 eigT=+1.0E+23 taz=270 tpl=45 paz=90 ppl=45
#mtcomp verbose Mo=0.0E+0 pdc=0.0 clvd_type=+v pclvd=0.0 piso=0.0 eigP=-6.74E+23 eigB=+3.85E+23 eigT=+5.89E+23 taz=219 tpl=18 paz=128 ppl=4

#mtcomp verbose Mo=2.0E+23 pdc=0.0 clvd_type=+v pclvd=0.0 piso=1.0 eigP=-0.5E+23 eigB=-0.5E+23 eigT=+1.0E+23 taz=160 tpl=70 paz=70 ppl=0

### example mixed % isotropic opening crack
###
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.03 pclvd=0.07 piso=0.90 >! out.90
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.07 pclvd=0.13 piso=0.80 >! out.80
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.10 pclvd=0.20 piso=0.70 >! out.70
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.13 pclvd=0.27 piso=0.60 >! out.60
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.15 pclvd=0.30 piso=0.55 >! out.55
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.17 pclvd=0.33 piso=0.50 >! out.50
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.18 pclvd=0.37 piso=0.45 >! out.45
# mtcomp noverbose Mo=1.0E+23 str=45 dip=45 rak=90 clvd_type=+v pdc=0.20 pclvd=0.40 piso=0.40 >! out.40
