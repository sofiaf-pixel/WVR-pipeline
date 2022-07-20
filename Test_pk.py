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
from matplotlib.dates import DateFormatter
from scipy.signal import butter, filtfilt


month_list=['202004', '202001', '201909',  '201904']
#month_list=['201909']

for month in month_list:


    double_pkfn='doubmemod_data/doublemod_'+month+'.txt'

    print(double_pkfn)

    f = open(double_pkfn,'rb')
    doublemod_dict=pk.load(f)
    f.close()

    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.title('Double Mod Amp')
    pl.show()

    pl.plot(doublemod_dict['date'], doublemod_dict['offset0'])
    pl.plot(doublemod_dict['date'], doublemod_dict['offset1'])
    pl.plot(doublemod_dict['date'], doublemod_dict['offset2'])
    pl.plot(doublemod_dict['date'], doublemod_dict['offset3'])
    pl.title('Atm T')
    pl.show()

    print(doublemod_dict)
