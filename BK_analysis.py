import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import matplotlib.dates as mdates
import datetime
import scipy.io
import h5py
import mat73
import plotAzscan as pA
import pickle as pk
import ephem
from astropy.time import Time



raz=pA.ReadAzscan()

class extract_ts(object):

    def __init__(self, unit=None, verb=True):

        '''
        '''
        self.lat = -90 #deg
        self.lon = 0 #deg
        self.alt = 2835 #m

    def keck_bandpasses(self, bp_file='K270_frequency_spectrum_20170710.txt', path_to_file='../txtfiles/'):
        file_variable = open(path_to_file+bp_file)
        lines = file_variable.readlines()
        freqs = []
        response_unform = []
        response_rj = []

        for i in range (9, len(lines)):
            var = lines[i].split(',')
            freqs.append(float(var[0]))
            response_unform.append(float(var[1]))
            response_rj.append(float(var[2]))

        # pl.plot(freqs, response_rj)
        # pl.plot(freqs, response_unform)
        # pl.title(bp_file)
        # pl.show()

        return freqs, response_unform, response_rj


    def make_fp_map(self, tod, rx):

        #ut = 2455822.20000367 #julian date

        # Which Julian Date does Ephem start its own count at?
        J0 = ephem.julian_date(0)
        t = Time(tod.t, format='datetime')

        observer = ephem.Observer()
        observer.lon = str(self.lon)  # str() forces deg -> rad conversion
        observer.lat = str(self.lat)  # deg -> rad
        observer.elevation = self.alt
        #observer.date = ut - J0
        observer.date = t.jd - J0

        x_az=tod.pointing.hor.az

        ra, dec = observer.radec_of(az, el)

        az = 3.30084818 #rad
        el = 0.94610742 #rad


        ts_cal, ts_cal_p3 = self.calib_tod_rx(tod, rx)
        a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)
        t=np.array(tod.std)
        #t=t[tod.mapind] #not using mapind as fs.sf;fs.se should take care of this

        #x_az=x_az[tod.mapind]
        #T_rx=ts_cal.psum[tod.mapind,det]

        if i_det == 'all':
            i_det_list = [j for j in range (0, len(ts_cal.pdiff[1,:]))]
        else:
            i_det_list = [i_det]

        for i_det in i_det_list:

            if p3==False:
                T_rx_pdiff=ts_cal.pdiff[:,i_det]
                T_rx_psum=ts_cal.psum[:,i_det]
            elif p3==True:
                T_rx_pdiff=ts_cal_p3.pdiff[:,i_det]
                T_rx_psum=ts_cal_p3.psum[:,i_det]
            fs=tod.fs

            x_sum, D_sum=raz.interpToImage_BK(x_az, T_rx_psum, fs)
            x_diff, D_diff=raz.interpToImage_BK(x_az, T_rx_pdiff, fs)
