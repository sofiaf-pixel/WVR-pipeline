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
from scipy import stats
from scipy import signal
import Read_BICEP_ts as bts



lag=np.arange(2*nt-1)-nt
pl.figure(figsize=(15,8))

npts=49
delay=np.zeros(npts)
phase=np.arange(npts)*2

for k in range(npts):
    rowcor1=np.correlate(atmogram[2*k,:], atmogram[2*k+2,:])
    delay[k]=rowcor1.argmax()


pl.plot(phase, nt=delay, '-bs', markersize=10, linewidth=0.7)
pl.xlim()
pl.xlabel('Azimuthal phase')
pl.ylabel('Delay')
