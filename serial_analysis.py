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
from zipfile import ZipFile
import tarfile
import math


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()


def extract_PWV_atmo_folder(folder, template='SPole_annual_50.amc', out_folder='serial_plots', posting_folder='None'):

    data_folder='../../wvr1_data_local/'
    for flist in os.listdir(data_folder+folder):
        try:
            print('file:', flist)
            path_to_test=data_folder+folder+'/'+flist
            #path_to_pickle='am_datfiles_Az/'+template[:-4]+'/'+flist[:-4]
            #pickle_fn_temps_impmod=path_to_pickle+'/'+flist[:-4]+'_clean_mod3_clean_method_import_model_pickle_temps.txt'
            #pickle_fn_temps_fit=path_to_pickle+'/'+flist[:-4]+'_clean_mod3_clean_method_fit_pickle_temps.txt'

            print('Starting fit.')
            x_am.fit_w_am_Az(flist, clean_method='import_model', out_path=out_folder, template=template)
            #x_am.fit_w_am_Az(flist, clean_method='fit')
        except:
            print(flist + 'failed.')




folder='20200418_tag'
extract_PWV_atmo_folder(folder)
