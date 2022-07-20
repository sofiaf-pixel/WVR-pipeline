import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import numpy as np
import pickle as pk
import pylab as pl
import scipy.optimize as sp
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec
from itertools import groupby
from operator import itemgetter
from math import atan, atan2
import math
import Read_BICEP_ts as bts
import plotAzscan as pA
import AM_functions as am
import plotElnod as pE
import pvlib
import pandas
import datetime
from astropy.modeling import models, fitting, functional_models
import pointing as pnt

xp=pnt.pointing()

speed='slower' #'slower' #change just this

if speed=='slower':
    day='06'
elif speed =='nominal':
    day='05'

#test_file_fast='20181209_210002_scanAz_beamMap_fast.txt'
#test_file_slow='20181209_210002_scanAz_beamMap_slow.txt'
#test_file_slow='20181210_160002_scanAz_beamMap_slow.txt'
#test_file_slow='20181210_140002_scanAz_beamMap_slow.txt'
#test_file_slow='20181210_120002_scanAz_beamMap_slow.txt'
#test_file_slow='20181209_140002_scanAz_beamMap_slow.txt'
# test_file_slow='20181209_180002_scanAz_beamMap_slow.txt'

#path_to_test='pointing_data/'+test_file_slow[:-9]+'/'


#pf='../../../Postings/WVR_postings/20210803_first_pointing_model/plots'
# pf2='../../../Postings/WVR_postings/20210925_pointing_model2_old/plots'
# pf2_fast='../../../Postings/WVR_postings/20210925_pointing_model2/plots'

#2022 data
pf='../../../Postings/WVR_postings/20220106_beam_maps_2022/plots'
pf2='../../../Postings/WVR_postings/20220106_beam_maps_2022/plots'
pf2_fast='../../../Postings/WVR_postings/20220106_beam_maps_2022/plots'


#
# test_file_slow_list=['20181210_160002_scanAz_beamMap_slow.txt', '20181210_140002_scanAz_beamMap_slow.txt',
#             '20181210_120002_scanAz_beamMap_slow.txt', '20181209_140002_scanAz_beamMap_slow.txt', '20181209_180002_scanAz_beamMap_slow.txt', '20181209_210002_scanAz_beamMap_slow.txt']
#
# test_file_fast_list=['20181210_160002_scanAz_beamMap_fast.txt', '20181210_140002_scanAz_beamMap_fast.txt',
#             '20181210_120002_scanAz_beamMap_fast.txt', '20181209_140002_scanAz_beamMap_fast.txt', '20181209_180002_scanAz_beamMap_fast.txt', '20181209_210002_scanAz_beamMap_fast.txt']
#
#
# files_that_failed= ['20181208_120002_scanAz_beamMap_fast.txt', '20181208_130002_scanAz_beamMap_fast.txt', '20181208_140002_scanAz_beamMap_fast.txt', '20181209_090002_scanAz_beamMap_fast.txt']
#


#test_file_slow_list=['20181209_210002_scanAz_beamMap_slow.txt']

#path_to_test='wvr1_data/pointing_data/'#+test_file[:-9]+'/'
path_to_test='../../wvr1_data_local/pointing_data/'#+test_file[:-9]+'/'

test_file_slow_list_all=[]
test_file_fast_list_all=[]
fn_core_all_list=[]

for flist in os.listdir(path_to_test):
    print(flist[:8])
    if flist[:8]=='202201'+day:
        if flist[-7:]=='beamMap':
            test_file_slow_list_all.append(flist+'_slow.txt')
            test_file_fast_list_all.append(flist+'_fast.txt')
            fn_core_all_list.append(flist)
            print(flist)

#print(test_file_slow_list_all)
#
# #
xp.extract_avg_par(figfolder=pf2_fast, fn='pointing_data_2022_fast_'+speed+'.txt', showplots=1)

sys.exit()

print('processing beam files:', fn_core_all_list)

eloffs_list_slow=[]
azoffs_list_slow=[]
pointfile_list_slow=[]
sigmaX_list_slow=[]
sigmaY_list_slow=[]

eloffs_list_fast=[]
azoffs_list_fast=[]
pointfile_list_fast=[]
sigmaX_list_fast=[]
sigmaY_list_fast=[]

fn_that_failed=[]


# for map in test_file_fast_list:
#     az_offs_axis, el_list, D_map, D_cut = xp.make_beam_map(map, pf2_fast, sampling='fast', showplots=0) #tfast
#
# sys.exit()

#
#for test_file in files_that_failed[3:]:
#
# test_file_fast='20220105_000002_scanAz_beamMap_fast.txt'
#
# az_offs_axis, el_list, D_map, D_cut = xp.make_beam_map(test_file_fast, pf2_fast, sampling='fast', showplots=1)
#
#
# sys.exit()






for fn_core in fn_core_all_list:
    try:
        #fn_core=test_file[:-9]
        print('Making beam map for '+fn_core)
        test_file_slow=fn_core+'_slow.txt'
        test_file_fast=fn_core+'_fast.txt'

        az_offs_axis, el_list, D_map, D_cut = xp.make_beam_map(test_file_fast, pf2_fast, sampling='fast', scan_speed=speed, showplots=0) #tfast

        print('az_offs_axis, el_list=',az_offs_axis, el_list)

        eloffs_center, azoffs_center, sigmaX, sigmaY= xp.fit_sun_2D(az_offs_axis, el_list, D_map, test_file_fast, pf2_fast, scan_speed=speed, showplots=0) #fast
        eloffs_list_fast.append(eloffs_center)
        azoffs_list_fast.append(azoffs_center)
        sigmaX_list_fast.append(sigmaX)
        sigmaY_list_fast.append(sigmaY)
        pointfile_list_fast.append(test_file_fast[:-24])
        print('fast output parameters:', eloffs_center, azoffs_center, sigmaX, sigmaY)

        # az_offs_axis, el_list, D_map, D_cut = xp.make_beam_map(test_file_slow, pf2_fast, showplots=0) #slow
        # eloff_center, azoffs_center, sigmaX, sigmaY= xp.fit_sun_2D(az_offs_axis, el_list, D_map, test_file_slow, pf2_fast, showplots=1) #slow
        # eloffs_list_slow.append(eloffs_center)
        # azoffs_list_slow.append(azoffs_center)
        # sigmaX_list_slow.append(sigmaX)
        # sigmaY_list_slow.append(sigmaY)
        # pointfile_list_slow.append(test_file_slow[:-24])
        # print('slow output parameters:', eloffs_center, azoffs_center, sigmaX, sigmaY)

        print('beam map '+fn_core+' done.')
    except:
        print('beam map '+fn_core+' failed.')
        fn_that_failed.append(fn_core)



D_p_fast = {'fn':[], 'az_off': [], 'el_off':[], 'az_sigma':[], 'el_sigma':[]}
D_p_slow = {'fn':[], 'az_off': [], 'el_off':[], 'az_sigma':[], 'el_sigma':[]}

D_p_fast['fn']=pointfile_list_fast
D_p_fast['az_offs']=azoffs_list_fast
D_p_fast['el_offs']=eloffs_list_fast
D_p_fast['az_sigma']=sigmaY_list_fast
D_p_fast['el_sigma']=sigmaX_list_fast

D_p_slow['fn']=pointfile_list_slow
D_p_slow['az_offs']=azoffs_list_slow
D_p_slow['el_offs']=eloffs_list_slow
D_p_slow['az_sigma']=sigmaY_list_slow
D_p_slow['el_sigma']=sigmaX_list_slow

print('parameters:', D_p_fast)

f = open('pointing_data_2020_fast_'+speed+'.txt','wb')
pk.dump(D_p_fast, f)
f.close()

f = open('pointing_data_2020_slow_'+speed+'.txt','wb')
pk.dump(D_p_slow, f)
f.close()


print('List of files that failed:', fn_that_failed)
