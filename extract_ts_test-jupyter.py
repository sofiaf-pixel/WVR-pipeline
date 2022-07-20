import Read_BICEP_ts as bts
#import extract_ts as et
import plotAzscan as pA
import numpy as np
import pylab as pl
import os, sys
import pickle as pk
import AM_functions as am

from scipy.interpolate import interp1d
from scipy import interpolate
#import importlib #importlib.reload(packagename)

class struct(object):
    pass


wvr_scan = '20200418_190135_scanAz_fast.txt'
#tag='20200408E05_dk068'
bk_tag='20200418B01_dk203'
rx=210
det=10
rx_list = [30, 40, 210, 270]
res_150 = 1.5 #deg
rx_res = []
rx_list = [210]


raz=pA.ReadAzscan()
ets=bts.extract_ts(bk_tag, wvr_scan)
x_am=am.AM_functions()


# ets.plot_bk_atmo(rx_list, wvr_fn=wvr_scan)
# ets.plot_ts(rx_list, det)
# #ets.plot_PS(rx_list, det)
#
#wvr_struct=struct();
#wvr_struct.waz, wvr_struct.az, wvr_struct.calib_az, wvr_struct.fs, wvr_struct.idx, wvr_struct.pwv_ts, wvr_struct.D_pwv, wvr_struct.tot = raz.read_pwvatmo(wvr_scan)

#bk_struct=struct();
#tod = ets.load_tod(bk_tag)
#bk_struct.az, bk_struct.T_rx_psum, bk_struct.T_rx_pdiff, bk_struct.D_sum, bk_struct.D_diff = ets.pl_tod_atmo(bk_tag, tod, rx)

ets.correlate_az_tod(rx)
