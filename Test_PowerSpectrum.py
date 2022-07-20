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


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()
ets=bts.extract_ts()


def butter_highpass_filter(data, cutoff, fs, order=5):

    def butter_highpass(cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
        return b, a

    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y



def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):

    def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = signal.butter(order, [low, high], btype='band')
        return b, a

    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = signal.lfilter(b, a, data)
    return y


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def plot_WVR_ps(test_file, pf):

    path_to_test='BAElnod_data/'+test_file[:-9]+'/'

    #The data order  in FH is: TIMEWVR CH0C CH0A CH0H CH0B CH1C CH1A CH1H CH1B CH2C CH2A CH2H CH2B CH3C CH3A CH3H CH3B EL AZ
    D_clean, waz, mod_removed_data_fit, p_double, p_err_double, p_single, p_err_single, calib_data_Gavg, FH =x.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=3, clean_method='fit')
    D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single, calib_data_Gavg, FH =x.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=3, clean_method='import_model')
    D, waz_1, calib_data_Gavg, FH = x.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=0)

    TIMEWVR=np.array(FH[:,0])
    CH0C=np.array(FH[:,1])
    CH0A=np.array(FH[:,2])
    CH0H=np.array(FH[:,3])
    CH0B=FH[:,4]
    T0=np.array(FH[:,0])
    az=FH[:,18] #unwrapped

    dt_wvr=TIMEWVR[10]-TIMEWVR[9]
    daz_wvr=az[10]-az[9]
    samp_rate=1./dt_wvr
    tot_time=TIMEWVR[len(TIMEWVR)-1]-TIMEWVR[0]
    data_len=len(CH0A)

    print(test_file)
    print('time_len=', tot_time/dt_wvr)
    print('data_len=', data_len)

    maska=np.where((TIMEWVR > 0) & (CH0A > 0))
    maskb=np.where((TIMEWVR > 0) & (CH0B > 0))
    mask0=np.where((TIMEWVR > 0) & (T0 > 0))

    time_ordered_time=TIMEWVR#[mask]
    time_ordered_az=az#[mask]

    ps_cha=x.extract_PowerSpectrum(time_ordered_time, dt_wvr, time_ordered_az, daz_wvr, CH0A, data='WVR')
    ps_chb=x.extract_PowerSpectrum(time_ordered_time, dt_wvr, time_ordered_az, daz_wvr, CH0B, data='WVR')
    ps_T0=x.extract_PowerSpectrum(time_ordered_time, dt_wvr, time_ordered_az, daz_wvr, T0, data='WVR')

    # time_ordered_data0=calib_data_Gavg[:,1]
    #
    # l_scale=0.9813
    # #l_scale=0.95
    # #l_scale=1.
    #
    time_ordered_time=calib_data_Gavg[:,0]
    #time_ordered_time_shift=np.array([j/l_scale for j in time_ordered_time])

    samp_freq=1/(calib_data_Gavg[1,0]-calib_data_Gavg[0,0])

    #time_ordered_az=calib_data_Gavg[:,6]
    #time_ordered_az_shift=time_ordered_time_shift*(360./30.)
    time_ordered_az=time_ordered_time*(360./30.)

    pickle_fn='am_datfiles_Az/SPole_annual_50/'+test_file[:-4]+'/'+test_file[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'



    #
    # if not os.path.exists(pickle_fn):
    #     print('Fitting single and double mod to '+test_file)
    #     x_am.fit_w_am_Az(test_file, path_to_data=path_to_test)

    f = open(pickle_fn,'rb')
    fit_output = pk.load(f)
    f.close()
    pwv_ts=fit_output['pwv']

    az_pk=fit_output['Az'] #wrapped


    # convolve=True
    #
    #time_ordered_data=pwv_ts

    mask = np.isfinite(pwv_ts)
    print('mask =', mask)
    print('~mask=', ~mask)
    #mask = np.isfinite(pwv_ts)
    pwv_ts_masked = pwv_ts[mask]
    #pwv_ts_masked[~mask]=np.mean(pwv_ts)
    time_ordered_az_masked = time_ordered_az[mask]
    time_ordered_time_masked = time_ordered_time[mask]
    #time_ordered_az_shift_masked = time_ordered_az_shift[mask]
    #time_ordered_time_shift_masked = time_ordered_time_shift[mask]
    # pl.plot(np.gradient(pwv_ts_masked))
    # pl.show()


    ind_max=len(pwv_ts_masked)-1000
    ind_min=20
    #
    #
    ps_pwv=x.extract_PowerSpectrum(time_ordered_time_masked, dt_wvr, time_ordered_az_masked, daz_wvr, pwv_ts_masked)
    # ps_pwv_shift=x.extract_PowerSpectrum(time_ordered_time_shift_masked[ind_min:ind_max], time_ordered_az_shift_masked[ind_min:ind_max], pwv_ts_masked[ind_min:ind_max])
    #


    pl.axvline(x=1./360., c='g', alpha=0.5, label='dipole')
    pl.axvline(x=1./180., c='b', alpha=0.5, label='quadrupole')
    # pl.plot(ps_cha.az.freq, ps_cha.az.ps, c='blue', label='chA')
    # pl.plot(ps_chb.az.freq, ps_chb.az.ps, c='red', label='chB')
    #pl.plot(ps_T0.az.freq, ps_T0.az.ps, c='y', label='T0')
    pl.plot(ps_pwv.az.freq, ps_pwv.az.ps, c='green', label='pwv')
    pl.suptitle('Power Spectrum - T_ch0/PWV')
    pl.title(test_file[:-4])
    #pl.xlabel('1/t[s^(-1)]')
    pl.xlabel('1/az[deg^(-1)]')
    pl.ylabel('PS(T_uncal)[ADU^2]')
    pl.legend()
    pl.loglog()
    pl.savefig(pf+'/'+test_file+'_T0-PWV_PS.png')
    #pl.show()
    pl.close()



    #to plot filtered data
    #
    #
    # dt=time_ordered_time[10]/time_ordered_time[9]
    # daz=time_ordered_az[10]/time_ordered_az[9]
    # print('dt=', dt)
    # samp_rate=1./dt
    # print('samp_rate=', samp_rate)
    # fc=(1./90.)
    #
    # idx = np.argsort(ps_pwv.az.freq[np.where(ps_pwv.az.freq>0)])
    # ps_cut=ps_pwv.az.ps[idx]
    # freq_cut=ps_pwv.az.freq[idx]
    #
    # b=0.8*fc
    #
    # lowcut=0.1
    # highcut=0.9

    #ps_pwv_filtered = butter_highpass_filter(ps_pwv.az.ps, fc, fs=samp_rate)
    # pwv_hpfiltered = butter_highpass_filter(pwv_ts[mask], fc, fs=samp_rate)
    # pwv_bpfiltered = butter_bandpass_filter(ps_pwv.az.ps, lowcut, highcut, fs=samp_rate)
    #
    # ps_pwv_hpfiltered = x.extract_PowerSpectrum(time_ordered_time[mask], time_ordered_az[mask], pwv_hpfiltered)
    # ps_pwv_bpfiltered = x.extract_PowerSpectrum(time_ordered_time[mask], time_ordered_az[mask], pwv_bpfiltered)

    #
    # kernel = np.ones(120)
    # kernel_gauss=gaussian(np.linspace(0,len(kernel),len(kernel)), len(kernel)/2, len(kernel)/8)
    # ps_pwv_convolved_az = np.convolve(ps_pwv.az.ps, kernel_gauss, mode='same')
    # ps_pwv_convolved_az_shift = np.convolve(ps_pwv_shift.az.ps, kernel_gauss, mode='same')
    # ps_pwv_convolved_t = np.convolve(ps_pwv.t.ps, kernel_gauss, mode='same')
    #




def plot_BK_ps(tag, rx, det, pf, p3=True):

    tod = ets.load_tod(tag)
    x_az=tod.pointing.hor.az


    ts_cal, ts_cal_p3 = ets.calib_tod_rx(tod, rx)
    t=np.array(tod.std)

    time_ordered_az=x_az#[tod.mapind]

    if p3==True:
        time_ordered_data_pdiff=ts_cal_p3.pdiff[:, det]#[tod.mapind, det]
        time_ordered_data_psum=ts_cal_p3.psum[:, det]#[tod.mapind, det]
    else:
        time_ordered_data_pdiff=ts_cal.pdiff[:, det]#[tod.mapind, det]
        time_ordered_data_psum=ts_cal.psum[:, det]#[tod.mapind, det]

    time_ordered_time=t#[tod.mapind]
    #

    t_cut=t[tod.mapind]
    t_totsec=[(now-t_cut[0]).total_seconds() for now in t_cut]


    #time_ordered_az=t_totsec[len(t_totsec)-1]*scan_speed_BK
    fs=tod.fs
    x_int, waz = x.return_waz(time_ordered_az, fs)


    # x_az_map=x_az#[tod.mapind]
    # no_zero=np.where(waz != 0)

    dt=t_totsec[1]-t_totsec[0]
    az_map=x_az[tod.mapind]
    daz=az_map[1]-az_map[0]

    #non masked
    ps_BK_pdiff=x.extract_PowerSpectrum(t, dt, waz, daz, time_ordered_data_pdiff, data='BK')
    ps_BK_psum=x.extract_PowerSpectrum(t, dt, waz, daz, time_ordered_data_psum, data='BK')
    #masked
    # ps_BK_pdiff=x.extract_PowerSpectrum(time_ordered_time_masked, dt, waz_masked, daz, pdiff_masked, data='BK')
    # ps_BK_psum=x.extract_PowerSpectrum(time_ordered_time, dt, waz, daz, time_ordered_data_psum, data='BK')
    s0 = int(fs.sf[0]);e0=int(fs.ef[0])

    #tag_size = x_az[e0]-x_az[s0]
    tag_size = x_az[e0]-x_az[s0]


    print('tag_size = ', tag_size)

    fig, ax1 = pl.subplots()#, sharex=True)
    ax1.axvline(x=1./tag_size, c='k', ls='--', alpha=0.5, label='tag size')
    ax1.plot(ps_BK_pdiff.az.freq, ps_BK_pdiff.az.ps, c='blue', label='pdiff')
    ax1.plot(ps_BK_psum.az.freq, ps_BK_psum.az.ps, c='red', label='psum')
    ax1.set_xlabel('1/az[deg^(-1)]')
    ax1.set_ylabel('T_PS[K^2]')
    ax1.loglog()
    ax1.legend(loc=0)
    # ax2 = ax1.twiny()
    # ax2.plot(ps_BK_pdiff.t.freq, ps_BK_pdiff.az.ps, c='blue', label='pdiff')
    # ax2.plot(ps_BK_psum.t.freq, ps_BK_psum.az.ps, c='red', label='psum')
    # ax2=ax1.twinx()
    # ax2.plot(ps_pwv.az.freq, ps_pwv.az.ps, c='green', label='pwv')
    #ax.set_ylabel('PWV_PS[um^2]')
    # ax2.axvline(x=1./360., c='g', ls='--', alpha=0.5, label='dipole')
    pl.suptitle('Power Spectrum')
    pl.title(tag)
    # ax2.legend(loc=0)
    # ax2.loglog()
    if p3==True:
        pl.savefig(pf+'/'+tag+'_rx_'+str(rx)+'_ab-idx_'+str(det)+'_PS.png')
    else:
        pl.savefig(pf+'/'+tag+'_rx_'+str(rx)+'_ab-idx_'+str(det)+'_PS_nop3.png')
    #pl.show()
    pl.close()






pf='../../Postings/WVR_postings/20210720_BK_PS/plots'


#test_file='20190103_150135_scanAz_fast.txt'
test_file='20200418_140135_scanAz_fast.txt'
#test_file='20200101_230134_scanAz_fast.txt'
#test_file='20190101_010135_scanAz_fast.txt'

test_file_list=['20200418_140135_scanAz_fast.txt', '20200418_130135_scanAz_fast.txt', '20200418_150134_scanAz_fast.txt', '20200418_190135_scanAz_fast.txt']


#BK scan
tag='20200418B01_dk203'

rx_list = [210, 270]
#rx=210
det_list = [6,7,8,9,10]
#det=10


#plot_WVR_ps(test_file, pf)

for rx in rx_list:
    for det in det_list:
        plot_BK_ps(tag, rx, det, pf)




# fig=pl.figure()
# #pl.axvline(x=1./tag_dt, c='r', alpha=0.8, label='tag size')
# pl.plot(ps_BK_pdiff.t.freq, ps_BK_pdiff.t.ps, c='k', alpha=0.6, label='PS')
# pl.suptitle('Power Spectrum')
# pl.title(tag)
# pl.xlabel('1/t[Hz]')
# pl.ylabel('PS(T)[K^2]')
# pl.legend()
# pl.loglog()
#
# pl.show()
#
#
# fig1, ax1 = pl.subplots(4,1)
# fig2, ax2 = pl.subplots(4,1)
#
# for i in range(4):
#
#     time_ordered_data=calib_data_Gavg[ind_min:ind_max,i+1]
#     time_ordered_data_modremoved=mod_removed_data[ind_min:ind_max,i]
#     time_ordered_data_modremoved_fit=mod_removed_data_fit[ind_min:ind_max,i]
#
#     ps_T=x.extract_PowerSpectrum(time_ordered_time[ind_min:ind_max], time_ordered_az[ind_min:ind_max], time_ordered_data)
#     ps_T_modremoved=x.extract_PowerSpectrum(time_ordered_time[ind_min:ind_max], time_ordered_az[ind_min:ind_max], time_ordered_data_modremoved)
#     ps_T_modremoved_fit=x.extract_PowerSpectrum(time_ordered_time[ind_min:ind_max], time_ordered_az[ind_min:ind_max], time_ordered_data_modremoved_fit)
#
#     kernel = np.ones(120)
#     kernel_gauss=gaussian(np.linspace(0,len(kernel),len(kernel)), len(kernel)/2, len(kernel)/8)
#     ps_T_convolved_az = np.convolve(ps_T_modremoved.az.ps, kernel_gauss, mode='same')
#
#     ax1[i].axvline(x=1./360., c='g', alpha=0.5, label='single mod')
#     ax1[i].axvline(x=1./180., c='b', alpha=0.5, label='double mod')
#     #ax1[i].plot(ps_T0.az.freq, ps_T0.az.ps, c='k', label='raw')
#     ax1[i].plot(ps_T_modremoved.az.freq, ps_T_modremoved.az.ps, c='k', label='data')
#     ax1[i].plot(ps_T_modremoved.az.freq, ps_T_convolved_az/len(kernel_gauss), c='y', label='data convolved')
#     #ax1[i].plot(ps_T0_modremoved_fit.az.freq, ps_T0_modremoved_fit.az.ps, c='b', label='mod_removed_scanfit')
#     ax1[0].set_title('Power Spectrum\n'+test_file[:-4])
#     ax1[i].set_xlabel('1/az[deg^(-1)]')
#     ax1[i].set_ylabel('PS(T_ch'+str(i)+')[K^2]')
#     ax1[3].legend()
#     ax1[i].loglog()
#
#     ax2[i].axvline(x=1./30., c='g', alpha=0.5, label='single mod')
#     ax2[i].axvline(x=1./15., c='b', alpha=0.5, label='double mod')
#     ax2[i].plot(ps_T.t.freq, ps_T.t.ps, c='k', label='raw')
#     ax2[i].plot(ps_T_modremoved.t.freq, ps_T_modremoved.t.ps, c='r', label='mod_removed_monthlymodel')
#     ax2[i].plot(ps_T_modremoved_fit.t.freq, ps_T_modremoved_fit.t.ps, c='b', label='mod_removed_scanfit')
#     ax2[0].set_title('Power Spectrum\n'+test_file[:-4])
#     ax2[i].set_xlabel('1/t[s^(-1)]')
#     ax2[i].set_ylabel('PS(T_ch'+str(i)+')[K^2]')
#     ax2[3].legend()
#     ax2[i].loglog()



#just for T0

# i=0
#
# time_ordered_data=calib_data_Gavg[ind_min:ind_max,i+1]
# time_ordered_data_modremoved=mod_removed_data[ind_min:ind_max,i]
# time_ordered_data_modremoved_fit=mod_removed_data_fit[ind_min:ind_max,i]
#
# ps_T=x.extract_PowerSpectrum(time_ordered_time[ind_min:ind_max], time_ordered_az[ind_min:ind_max], time_ordered_data)
# ps_T_modremoved=x.extract_PowerSpectrum(time_ordered_time[ind_min:ind_max], time_ordered_az[ind_min:ind_max], time_ordered_data_modremoved)
# ps_T_modremoved_fit=x.extract_PowerSpectrum(time_ordered_time[ind_min:ind_max], time_ordered_az[ind_min:ind_max], time_ordered_data_modremoved_fit)
#
# kernel = np.ones(120)
# kernel_gauss=gaussian(np.linspace(0,len(kernel),len(kernel)), len(kernel)/2, len(kernel)/8)
# ps_T_convolved_az = np.convolve(ps_T_modremoved.az.ps, kernel_gauss, mode='same')
#
#
# fig=pl.figure()
#
# pl.axvline(x=1./360., c='g', alpha=0.8, label='l=1')
# pl.axvline(x=1./180., c='b', alpha=0.8, label='l=2')
# pl.axvline(x=1./120., c='c', alpha=0.8, label='l=3')
# pl.axvline(x=1./90.,  alpha=0.8, label='l=4')
# pl.plot(ps_T.az.freq, ps_T.az.ps, c='k', label='raw')
# pl.plot(ps_T_modremoved.az.freq, ps_T_modremoved.az.ps, c='r', label='mod_removed_monthlymodel')
# pl.plot(ps_T_modremoved_fit.az.freq, ps_T_modremoved_fit.az.ps, c='b', label='mod_removed_scanfit')
# pl.suptitle('Power Spectrum')
# pl.title(test_file[:-4])
# pl.xlabel('1/az[deg^(-1)]')
# pl.ylabel('T0[K^2]')
# pl.legend()
# pl.loglog()
#
# pl.show()
