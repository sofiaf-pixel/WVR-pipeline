import plotElnod as pE
import plotAzscan as pA
import AM_functions as am
import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time
import datetime
import matplotlib.dates as mdates
from math import *
import extract_mod_par as emp
import Read_BICEP_ts as bts
from dateutil import parser


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()
x_emp = emp.mod_parameters()

#WVR options

wvr_scan='20200921_150134_scanAz_fast.txt'

template='SPole_annual_50.amc'

month = wvr_scan[:6]


#0. Unzip data
#x_emp.unzip(month)

#1. Extract Modulation Parameters

# if not os.path.exists('mod_param/tilt_angle_fullmonth'+month+'.txt'):
#     x_emp.ModAmptoTilt(month)
# else:
#     print('mod_param/tilt_angle_fullmonth'+month+'.txt already exists.\n')
#     answer = input("Do you want to overwrite it? ")
#     if answer == "y":
#         x_emp.ModAmptoTilt(month)
#     elif answer == "n":
#         print('Using existing tilt model.')
#     else:
#         print("Please enter y or n.")


#2. Extract PWv atmogram from Temperatures
path_to_pwv = 'am_datfiles_Az/SPole_annual_50/'+wvr_scan[:-4]+'/'+wvr_scan[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'
if not os.path.exists(path_to_pwv):
    #x_am.fit_w_am_Az(wvr_scan, clean_method='import_model')
    data=x_am.fit_w_am_Az(wvr_scan, clean_method='fit')
else:
    print(path_to_pwv+' already exists.\n')
    answer = input("Do you want to overwrite it? ")
    if answer == "y":
        #x_am.fit_w_am_Az(wvr_scan, clean_method='import_model')
        data=x_am.fit_w_am_Az(wvr_scan, clean_method='fit')
    elif answer == "n":
        print('Converting PWV to Trj.')
    else:
        print("Please enter y or n.")
