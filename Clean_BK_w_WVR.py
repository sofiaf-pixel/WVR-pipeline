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
import scipy.stats

raz=pA.ReadAzscan()
x_am=am.AM_functions()


class struct(object):
    pass


pf='../../../Postings/WVR_postings/20220310_BK_WVR_datacleaning/plots/'

wvr_scan='20200921_150134_scanAz_fast.txt'
bk_tag='20200921B09_dk293'


wvr_scan1='20200415_140134_scanAz_fast.txt'
bk_tag1='20200415D03_dk293'

wvr_scan2='20200418_190135_scanAz_fast.txt'
bk_tag2='20200418B01_dk203'

wvr_scan3='20200723_190134_scanAz_fast.txt'
bk_tag3='20200723B09_dk203'

wvr_scan4='20200921_150134_scanAz_fast.txt'
bk_tag4='20200921B09_dk293'

wvr_scan_list=[wvr_scan2, wvr_scan3, wvr_scan4]
bk_tag_list=[bk_tag2, bk_tag3, bk_tag4]

# wvr_scan_list=[wvr_scan3, wvr_scan2]
# bk_tag_list=[bk_tag3, bk_tag2]
#
# for i_scan in range(len(wvr_scan_list)):
#     wvr_scan=wvr_scan_list[i_scan]
#     bk_tag=bk_tag_list[i_scan]
#
#     #Extracting Trj atmograms from WVR PWV
#     D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=1, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
#
# sys.exit()


for i_scan in range(1):
    wvr_scan=wvr_scan_list[i_scan]
    bk_tag=bk_tag_list[i_scan]

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
        pl.savefig(pf+'PlanckCMBu_to_RJu_coeff.png')

        pl.close()

        f_v = {'210':f_210, '270':f_270}

        return f_v



    def scale(im, nR, nC):
        number_rows = len(im)     # source number of rows
        number_columns = len(im[0])  # source number of columns
        return [[ im[int(number_rows * r / nR)][int(number_columns * c / nC)]
                     for c in range(nC)] for r in range(nR)]



    ets=bts.extract_ts(bk_tag, wvr_scan)
    f_v = find_planck_to_rj_coeff()


    #x_az=tod.pointing.hor.az
    waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan, show=0)#, az_lim=(np.min(x_az), np.max(x_az)))
    time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

    #Extracting Trj atmograms from WVR PWV
    D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

    wvr_atmo = ets.wvr_struct.D_pwv
    az_wvr = ets.wvr_struct.az_real #az_calib
    time_ordered_az=ets.wvr_struct.az_wvr
    fs = ets.wvr_struct.fs
    t_wvr=ets.wvr_struct.tot

    tod = ets.bk_struct
    fs_bk = tod.fs

    for rx in [270]:
        for gs_filt in [False]:
            for p3_filt in [True, False]:

                wvr_atmo_Trx = wvr_atmo_Trj[str(rx)]#-np.nanmean(wvr_atmo_Trj[str(rx)])

                x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol, el_fpu_center= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png', el_lim=10)
                x_b = x_pair[np.where(det_pol=='b')]
                x_a = x_pair[np.where(det_pol=='a')]
                y_b = y_pair[np.where(det_pol=='b')]
                y_a = y_pair[np.where(det_pol=='a')]


                for det in range(len(det_a_list)[:30]):

                    #try:

                    a_det = det_a_list[det]
                    i_det_b=np.where(x_b==x_a[det])[0]
                    b_det = det_b_list[i_det_b[0]]

                    #az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det)

                    corr_list_p3True=[]
                    corr_list_p3False=[]

                    az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, gs=gs_filt, i_det_savefig=det, posting_folder='None')
                    az_bk_pairs, T_rx_polA_bk, T_rx_polB_bk, D_polA_bk, D_polB_bk = ets.pl_tod_atmo_pairs(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], i_det_savefig=det, posting_folder=pf)

                    t_bk=np.array(tod.std)

                    time_UTC=[]
                    time_UTC=[parser.parse(time_str) for time_str in time_UTC_str]

                    fs_wvr=fs
                    full_time_wvr=[time_UTC[int(fs_s_i)] for fs_s_i in fs_wvr.s]
                    xlabels_wvr=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time_wvr])
                    nscans_wvr = len(fs_wvr.s)
                    xticks_wvr = np.arange(nscans_wvr)
                    x_sum, D_sum = raz.interpToImage_BK(az_bk, T_rx_psum_bk, fs_bk)
                    x_diff, D_diff = raz.interpToImage_BK(az_bk, T_rx_pdiff_bk, fs_bk)

                    full_time=[t_bk[int(fs_s_i)] for fs_s_i in fs_bk.sf]
                    xlabels=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time])
                    nscans_bk = len(fs_bk.sf)
                    xticks = np.arange(nscans_bk)

                    full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])
                    full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr])

                    t_min_wvr=datetime.time(full_time_wvr[0])
                    t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
                    t_min_bk=datetime.time(full_time[0])
                    t_max_bk=datetime.time(full_time[len(full_time)-1])

                    calib_az=np.array(calib_az)
                    #
                    # if  any(caz <0 for caz in calib_az):
                    #     print('adding 360')
                    #     calib_az=[i+360 for i in calib_az]


                    x_min_wvr_idx = np.where(calib_az>= np.min(x_diff))[0][0]
                    x_max_wvr_idx = np.where(calib_az>= np.max(x_diff))[0][0]


                    bk_time_mask = np.where((full_time_dt<=t_max_wvr)&(full_time_dt>=t_min_wvr))[0]
                    wvr_time_mask = np.where((full_time_wvr_dt>=t_min_bk)&(full_time_wvr_dt<=t_max_bk))[0]


                    xlabels=np.array(xlabels)
                    xlabels_wvr=np.array(xlabels_wvr)

                    xticks_mask=xticks[bk_time_mask]
                    xlabels_mask=xlabels[bk_time_mask]
                    xticks_wvr_mask=xticks_wvr[wvr_time_mask]
                    xlabels_wvr_mask=xlabels_wvr[wvr_time_mask]

                    D_sum_bk_tcut = D_sum_bk[:,bk_time_mask]
                    D_diff_bk_tcut = D_diff_bk[:,bk_time_mask]

                    D_polA_bk_tcut = D_polA_bk[:,bk_time_mask]
                    D_polB_bk_tcut = D_polB_bk[:,bk_time_mask]

                    wvr_atmo_tcut = wvr_atmo[:,wvr_time_mask]
                    wvr_atmo_Trx_tcut = wvr_atmo_Trx[:,wvr_time_mask]

                    wvr_atmo_Tplanck_tcut = np.array(wvr_atmo_Trx_tcut)/f_v[str(rx)]


                    #cut in Az range
                    az_wvr_max = np.max(calib_az)
                    az_wvr_min = np.min(calib_az)

                    az_wvr = np.linspace(az_wvr_min, az_wvr_max, np.shape(wvr_atmo_Tplanck_tcut)[0])


                    idx_s = np.where(az_wvr>=np.min(x_diff))[0][0]
                    idx_e = np.where(az_wvr>=np.max(x_diff))[0][0]

                    wvr_atmo_Tplanck_tcut_azcut = wvr_atmo_Tplanck_tcut[idx_s:idx_e, :]

                    print('wvr_time_mask = ', wvr_time_mask)
                    print('idx_s, idx_e = ', idx_s, idx_e)
                    print('az_wvr = ', az_wvr)
                    print('x_diff = ', x_diff)


                    if (t_min_bk <= t_min_wvr):
                        extent_wvr_cut = [0, len(wvr_time_mask), np.min(az_wvr[idx_s:idx_e]), np.max(az_wvr[idx_s:idx_e])]
                    else:
                        extent_wvr_cut = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(az_wvr[idx_s:idx_e]), np.max(az_wvr[idx_s:idx_e])]


                    D_sum = struct();

                    D_sum.D_sum_bk_tcut = D_sum_bk_tcut
                    D_sum.bk_time_mask = bk_time_mask
                    D_sum.x_sum = x_sum
                    D_sum.xticks_mask = xticks_mask
                    D_sum.xlabels_mask = xlabels_mask

                    wvr_atmo_struct = struct();

                    wvr_atmo_struct.extent_wvr_cut = extent_wvr_cut
                    wvr_atmo_struct.xticks_wvr_mask = xticks_wvr_mask
                    wvr_atmo_struct.xlabels_wvr_mask = xlabels_wvr_mask

                    #for psum/pdiff

                    wvr_atmo_Tplanck_tcut_azcut_p0 = wvr_atmo_Tplanck_tcut_azcut-np.nanmean(wvr_atmo_Tplanck_tcut_azcut)
                    wvr_atmo_struct.wvr_atmo_Tplanck_tcut_azcut_p0 = wvr_atmo_Tplanck_tcut_azcut_p0

                    wvr_atmo_Tplanck_tcut_azcut_resized = scale(wvr_atmo_Tplanck_tcut_azcut_p0, np.shape(D_sum_bk_tcut)[0], np.shape(D_sum_bk_tcut)[1])

                    wvr_atmo_struct.wvr_atmo_Tplanck_tcut_azcut_resized = wvr_atmo_Tplanck_tcut_azcut_resized

                    D_diff = struct();

                    D_diff.D_diff_bk_tcut = D_diff_bk_tcut
                    D_diff.bk_time_mask = bk_time_mask
                    D_diff.x_diff = x_diff
                    D_diff.xticks_mask = xticks_mask
                    D_diff.xlabels_mask = xlabels_mask


                    D_polB = struct();
                    D_polA = struct();

                    D_polA.D_polA_bk_tcut = D_polA_bk_tcut
                    D_polB.D_polB_bk_tcut = D_polB_bk_tcut

                    D_dict = {'wvr_atmo':wvr_atmo_struct,'D_sum': D_sum, 'D_diff': D_diff, 'D_polA': D_polA, 'D_polB': D_polB, 'a_det': a_det, 'b_det': b_det}

                    D_folder = 'bk_wvr_clean_data/'
                    D_fn = D_folder+'BK-WVR_data_cleaning_date_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filt'+str(p3_filt)+'_groundsub_'+str(gs_filt)+'.txt'

                    f = open(D_fn,'wb')
                    pk.dump(D_dict, f)
                    f.close()
                    #
                    #
                    # except Exception as e:
                    #     print(e)
                    #     print('rx = ', rx)
                    #     print('det = ', det)
                    #     print('p3_filt = ', p3_filt)
                    #     print('gs_filt = ', gs_filt)
