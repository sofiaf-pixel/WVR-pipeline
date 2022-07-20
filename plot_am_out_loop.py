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
day_list=[]
month=str(sys.argv[1])
type='elnod'

for i in range(29, 30):
    if 0<i<10:
        day_list.append(month+'0'+str(i))
    elif i>10:
        day_list.append(month+str(i))

for day in day_list:
	string='python3 plot_am_out.py '+str(day)+' '+type
	os.system(string)
