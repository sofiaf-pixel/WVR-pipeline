import Read_BICEP_ts as bts
import plotAzscan as pA
import numpy as np
import pylab as pl
import os, sys
import pickle as pk
import AM_functions as am
from scipy.interpolate import interp1d
from scipy import interpolate
from datetime import *
from dateutil import parser
import mat73


raz=pA.ReadAzscan()
x_am=am.AM_functions()

class struct(object):
    pass


def make_clems_plots(tag, pf):
    data = mat73.loadmat('tod/ba/bicep_array_'+tag+'_alt_map_weight.mat')

    keys_array = data.keys()

    for key_i in keys_array:
        if key_i[0]=='d':
            data_i =  data[key_i]
            pl.figure(figsize=(14,10))
            pl.imshow(data_i.T, aspect='auto', interpolation='nearest', origin='lower', vmin=-10, vmax=10)
            pl.colorbar()
            pl.xlabel('tag number (time)')
            pl.ylabel('good pair number')
            pl.title(tag)
            pl.savefig(pf+'all_det_ts_matrix'+key_i+'.png')
            pl.close()

def ps_az(az_ordered_az, az_ordered_data):

    az_ordered_az=az_ordered_az[np.where(np.isfinite(az_ordered_data))[0]]
    az_ordered_data=az_ordered_data[np.where(np.isfinite(az_ordered_data))[0]]

    az_ = az_ordered_az
    daz = az_ordered_az[10] - az_ordered_az[9]
    az = len(az_)
    df=1./daz
    fft_wvr = np.fft.fft(az_ordered_data)/az
    ps    = np.square(np.abs(fft_wvr))#*np.hanning(len(az_ordered_data))
    freq_ = np.fft.fftfreq(az,daz)

    ps_scaled=np.array([(ps[i]/freq_[i]) for i in range(len(ps))])

    idx_pos=np.where(freq_>=0)

    return freq_[idx_pos], ps_scaled[idx_pos]

def find_planck_to_rj_coeff():
    h=6.626176e-27 #erg*s
    k=1.380662e-16 #erg/K
    c=2.99792458e8 #m/s
    T_cmb=2.725 #K

    alpha = h/(k*T_cmb)

    freq_GHz=np.linspace(0, 300, 100000)
    freq_Hz=freq_GHz*(10**9)

    x = alpha*freq_Hz

    f_v = (x**2.)*np.exp(x)/((np.exp(x)-1)**2.)

    pl.figure(figsize=(12,8))
    pl.plot(freq_Hz, f_v)
    pl.title('T_rj/T_p')

    i210=np.where(freq_Hz>=210e9)[0][0]
    i270=np.where(freq_Hz>=270e9)[0][0]

    f_210=f_v[i210]
    f_270=f_v[i270]
    pl.axvline(x=270e9, c='r', label='270GHz')
    pl.axhline(y=f_270, c='r', label='f_270 ='+str(round(f_270,2)))
    pl.axvline(x=210e9, c='y', label='210GHz')
    pl.axhline(y=f_210, c='y', label='f_210 ='+str(round(f_210,2)))
    pl.legend()

    f_v = {'210':f_210, '270':f_270}

    return f_v




def wvr_PS(ets, rx_string, pf):

    if rx_string == 'PWV':
        waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan)
        wvr_atmo = D_pwv
    elif rx_string == 'T210':
        f_v = find_planck_to_rj_coeff()
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
        wvr_atmo = D210/f_v[rx_string[1:]] #from Trj to Tcmb
    elif rx_string == 'T270':
        f_v = find_planck_to_rj_coeff()
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
        wvr_atmo = D270/f_v[rx_string[1:]]


    az_wvr = ets.wvr_struct.az_real #az_calib
    fs = ets.wvr_struct.fs

    wvr_tod=[]
    wvr_ps=[]
    nscans=len(wvr_atmo[10,:])


    pl.figure(figsize=(12,8))
    for i in range(nscans-1):
        wvr_tod.append(wvr_atmo[:,i])
        az_ordered_az = az_wvr[fs.s[i]:fs.e[i]]
        az_ordered_data = wvr_atmo[:,i]
        freq, ps = ps_az(az_ordered_az, np.array(az_ordered_data))
        pl.plot(freq, ps, c='k', alpha=0.4)
        pl.scatter(freq, ps, c='blue')


    pl.axvline(x=1/20., c='y', label='20deg')
    pl.axvline(x=1/6., c='c', label='6deg')
    pl.axvline(x=1/3., ls='--', c='orange', label='3deg beam')
    pl.xlabel('1/az[deg^-1]')
    pl.ylabel('PS')
    pl.xlim((1/180.), 1)
    #pl.xscale('log')
    pl.loglog()
    pl.legend()
    pl.suptitle(rx_string+' - Power Spectrum')
    pl.title(wvr_scan[:-4])
    pl.savefig(pf+'PS_'+wvr_scan[:-4]+'-'+rx_string+'.png')
    pl.close()




def BK_PS(ets, pf, rx=210, p3_filt=False, groundsub=False):

    if groundsub==True:
        tod = ets.bk_struct
    elif groundsub == False:
        tod = ets.bk_struct_nogroundsub

    fs_bk = tod.fs

    x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')
    x_b = x_pair[np.where(det_pol=='b')]
    x_a = x_pair[np.where(det_pol=='a')]
    y_b = y_pair[np.where(det_pol=='b')]
    y_a = y_pair[np.where(det_pol=='a')]


    for det in range(len(det_a_list)-1):

        a_det = det_a_list[det]
        i_det_b=np.where(x_b==x_a[det])[0]
        b_det = det_b_list[i_det_b[0]]


        col_p3=['cyan', 'yellow']
        col=['blue', 'orange']

        az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=False, i_det_savefig=det, posting_folder='None')
        az_bk, T_rx_psum_bk_p3, T_rx_pdiff_bk_p3, D_sum_bk_p3, D_diff_bk_p3 = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=True, i_det_savefig=det, posting_folder='None')


        az_bk=np.array(az_bk)


        bk_sum_ps=[]
        bk_diff_ps=[]

        bk_sum_p3_ps=[]
        bk_diff_p3_ps=[]

        bk_freq=[]


        pl.figure(figsize=(12,8))
        for i in range(len(fs_bk.sf)):
            try:
                az_ordered_data_sum = T_rx_psum_bk[int(fs_bk.sf[i]):int(fs_bk.ef[i])]
                az_ordered_data_diff = T_rx_pdiff_bk[int(fs_bk.sf[i]):int(fs_bk.ef[i])]
                az_ordered_data_sum_p3 = T_rx_psum_bk_p3[int(fs_bk.sf[i]):int(fs_bk.ef[i])]
                az_ordered_data_diff_p3 = T_rx_pdiff_bk_p3[int(fs_bk.sf[i]):int(fs_bk.ef[i])]
                az_ordered_az = az_bk[int(fs_bk.sf[i]):int(fs_bk.ef[i])]

                freq_sum, ps_sum = ps_az(az_ordered_az, np.array(az_ordered_data_sum))
                freq_diff, ps_diff = ps_az(az_ordered_az, np.array(az_ordered_data_diff))
                freq_sum_p3, ps_sum_p3 = ps_az(az_ordered_az, np.array(az_ordered_data_sum_p3))
                freq_diff_p3, ps_diff_p3 = ps_az(az_ordered_az, np.array(az_ordered_data_diff_p3))

                bk_sum_ps.append(ps_sum)
                bk_diff_ps.append(ps_diff)
                bk_sum_p3_ps.append(ps_sum_p3)
                bk_diff_p3_ps.append(ps_diff_p3)

                bk_freq.append(freq_sum)


                pl.plot(freq_diff[1:], ps_diff[1:], c=col[1], alpha=0.2)
                pl.plot(freq_sum[1:], ps_sum[1:], c=col[0], alpha=0.2)
                if i==0:
                    pl.scatter(freq_sum, ps_sum, c=col[0], label='psum - p3filt = False')
                    pl.scatter(freq_diff, ps_diff, c=col[1], label='pdiff- p3filt = False')
                else:
                    pl.scatter(freq_sum, ps_sum, c=col[0])
                    pl.scatter(freq_diff, ps_diff, c=col[1])

                pl.plot(freq_diff_p3[1:], ps_diff_p3[1:], c=col_p3[1], alpha=0.2)
                pl.plot(freq_sum_p3[1:], ps_sum_p3[1:], c=col_p3[0], alpha=0.2)
                if i==0:
                    pl.scatter(freq_sum_p3, ps_sum_p3, c=col_p3[0], label='psum - p3filt = True')
                    pl.scatter(freq_diff_p3, ps_diff_p3, c=col_p3[1], label='pdiff- p3filt = True')
                else:
                    pl.scatter(freq_sum_p3, ps_sum_p3, c=col_p3[0])
                    pl.scatter(freq_diff_p3, ps_diff_p3, c=col_p3[1])

            except Exception as e:
                print(e)


        #pl.axvline(x=1/7., c='y', label='7deg')
        #pl.axvline(x=1/30., c='c', label='30deg')
        # pl.axvline(x=1/3., ls='--', c='orange', label='3deg beam')
        pl.xlabel('1/az[deg^-1]')
        pl.ylabel('PS')
        pl.loglog()
        pl.legend()
        pl.suptitle(rx_string+' -  gcp_idx:'+str(det_a_list[det])+'-'+str(det_b_list[det])+'\nPower Spectrum')
        pl.title(bk_tag)
        pl.savefig(pf+'PS_'+bk_tag+'-freq'+rx_string+'_det'+str(det)+'_p3'+str(p3_filt)+'_groundsub_'+str(groundsub)+'.png')
        pl.close()









pf = '../../../Postings/WVR_postings/20220310_BK_WVR_correlations_PS/plots/'


wvr_scan='20200921_150134_scanAz_fast.txt'
bk_tag='20200921B09_dk293'

ets=bts.extract_ts(bk_tag, wvr_scan)

make_clems_plots(bk_tag, pf)

rx_string_list = ['PWV', 'T210', 'T270']

for rx_string in rx_string_list:
    wvr_PS(ets, rx_string, pf)
    for p3 in [True, False]:
        for groundfilt in [True, False]:
            if rx_string[0]=='T':
                BK_PS(ets, pf, rx=float(rx_string[1:]), p3_filt=p3, groundsub=groundfilt)
