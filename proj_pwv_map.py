import Read_BICEP_ts as bts
import plotAzscan as pA
import numpy as np
import pylab as pl
import os, sys
import pickle as pk
import AM_functions as am
import math
import scipy.interpolate as interp
#import importlib #importlib.reload(packagename)

raz=pA.ReadAzscan()
ets=bts.extract_ts()
x_am=am.AM_functions()


#wvr_scan = '20200418_190135_scanAz_fast.txt'
#wvr_scan = '20190103_150135_scanAz_fast.txt'
#wvr_scan = '20200418_150134_scanAz_fast.txt'
wvr_scan = '20200418_140135_scanAz_fast.txt'


pf='../../Postings/WVR_postings/?/plots'


#x_am.plot_Trj_az(wvr_scan, rewrite_txt=0, out_pk=wvr_scan[:-4]+'_Trj_pk_BAK_0-360.txt')

# fn='tod/ba/'+wvr_scan+'_pwv_zoom_on.png'
# fn_posting=pf+'/'+wvr_scan+'_pwv_zoom_on.png'
#x

waz, az_pwv, fs, idx, pwv, D_pwv = raz.read_pwvatmo(wvr_scan, show=1)


scan_v=360./30.
t=waz/scan_v

#w_y=0.0005
w_y=0
w_x=0.1
#assuming wind just along x

el=55


x=np.zeros(len(waz))
y=np.zeros(len(waz))

for i in range (len(waz)):
    waz_i=waz[i]
    t_i=waz_i/scan_v
    (scannum_i, az_i) = divmod(waz_i,360)
    y[i]=(1./np.tan(math.radians(el)))*np.sin(math.radians(az_i))+ w_y*t_i
    x[i]=(1./np.tan(math.radians(el)))*np.cos(math.radians(az_i))+ w_x*t_i


pl.scatter(x, y, s=4, c=pwv)#, cmap='plasma')
pl.colorbar()
pl.suptitle('PWV atmogram projected')
pl.title(wvr_scan[:-4])
pl.show()


#sys.exit()

#for BK
v_y=0.0005 #redefining bcs diff pointing
v_x=0.1

tag='20200418B01_dk203'
tod = ets.load_tod(tag)
rx=210
det=10

el_bk=tod.pointing.hor.el
x_az, T_rx_psum, T_rx_pdiff, D_sum, D_diff = ets.pl_tod_atmo(tag, tod, rx)

T_rx_psum=T_rx_psum[tod.mapind]
T_rx_pdiff=T_rx_pdiff[tod.mapind]
x_az=x_az[tod.mapind]
el_bk=el_bk[tod.mapind]
t=np.array(tod.std)
t_cut=t[tod.mapind]
t_s=[(now-t_cut[0]).total_seconds() for now in t_cut]


x_bk=np.zeros(len(x_az))
y_bk=np.zeros(len(x_az))

el_avg=np.mean(el_bk)

for i in range (len(x_az)):
    y_bk[i]=(1./np.tan(math.radians(el_avg)))*np.sin(math.radians(x_az[i]))+ v_y*t_s[i]
    x_bk[i]=(1./np.tan(math.radians(el_avg)))*np.cos(math.radians(x_az[i]))+ v_x*t_s[i]



pl.scatter(x_bk, y_bk, s=4, c=T_rx_pdiff, cmap='plasma')
pl.colorbar()
pl.suptitle('T_bk '+str(rx)+' GHz projected - det_'+str(det))
pl.title(tag)
pl.show()
