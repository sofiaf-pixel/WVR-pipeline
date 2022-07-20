import plotElnod as pE
import plotAzscan as pA
import AM_functions as am
import pickle as pk
import pylab as pl
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time
import datetime
import matplotlib.dates as mdates
from math import *
import re


x_am=am.AM_functions()
raz=pA.ReadAzscan()

#scan=str(sys.argv[1])
#scan='20200418_150002_skyDip_fast'
# if scan[-7:-5]=='Az':
#     path='am_datfiles_Az/SPole_annual_50/'
#     type='azscan'
# elif scan[-7:-5]=='ip':
#     path='am_datfiles/SPole_annual_50/'
#     type='elnod'
# else:
#     print('Error in the given folder. Exiting.')
#     sys.exit()


day=str(sys.argv[1])
type=str(sys.argv[2])


# x_am.read_amc('20200418_190002_skyDip_fast_El55.025_am.dat.amc', pathtofn='am_datfiles/SPole_annual_50/20200418_190002_skyDip_fast/')
# sys.exit()
#Trial_folder

# if test==1:
#
#     os.makedirs(folder+'_test')
#     os.system('cp folder/*.amc')
#

remake_post_fig=0 #To overwrite/save plots automatically on the Posting plot folder

scan1='20200418_140135_scanAz_fast.txt'
scan2='20200405_090135_scanAz_fast.txt'
scan3='20200418_150134_scanAz_fast.txt'
scan4='20200418_190135_scanAz_fast.txt'
# Az_max=280
# Az_min=240

scan_list=[scan1]#, scan2, scan3]


Az_max=360
Az_min=0

if type=='elnod':
    x_am.plot_Trj_el(day, remake_post_fig=remake_post_fig, rewrite_txt=1)

elif type=='azscan':
    for scan in scan_list:
        print('out_pk=', scan[:-4]+'_Trj_pk_BAK_0-360.txt')
        x_am.plot_Trj_az(scan, Az_min=Az_min, Az_max=Az_max, remake_post_fig=remake_post_fig, rewrite_txt=1, out_pk=scan[:-4]+'_Trj_pk_BAK_0-360.txt')
        #x_am.plot_Trj_az(scan, Az_min=240, Az_max=241, remake_post_fig=remake_post_fig, rewrite_txt=0)
