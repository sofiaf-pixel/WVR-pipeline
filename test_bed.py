import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import scipy.interpolate as sint
import numpy as np
import pickle
import pylab as pl
from matplotlib import dates
import datetime
import scipy.optimize as sp
from pathlib import Path
import plotElnod as pE
import plotAzscan as pA
from time import perf_counter
import glob
import re as r
import BK_analysis as BK_an
import AM_functions as am
import Read_BICEP_ts as bts

re=pE.ReadElnod()
raz=pA.ReadAzscan()
x_am=am.AM_functions()
bk=BK_an.extract_ts()

wvr_scan='20200418_190135_scanAz_fast.txt'
bk_tag='20200418B01_dk203'

thesis_posting = '/Volumes/LaCie/BICEP/Postings/WVR_postings/20220720_thesis_plots/'


D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=1, out_pk='Trj_pk_BAK.txt', posting_folder=thesis_posting)
sys.exit()
ets=bts.extract_ts(bk_tag, wvr_scan)

#bk.keck_bandpasses()

# x_am.plot_Trj_az('20200418_190135_scanAz_fast.txt', Az_min=125, Az_max=185, remake_post_fig=0, rewrite_txt=1, datafolder='BAElnod_data/', posting_folder='')
#x_am.plot_Trj_az('20200418_190135_scanAz_fast.txt', Az_min=0, Az_max=360, remake_post_fig=0, rewrite_txt=1, rewrite_amc=0, datafolder='BAElnod_data/', posting_folder='')

x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol, el_offset= ets.det_fpu_location(270, fn_save = thesis_posting+bk_tag+'det_fpu_location_270.png', show_plots=1)
