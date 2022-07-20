
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

class mod_parameters(object):

    def __init__(self, bk_tag, wvr_scan, unit=None, verb=True):

        '''


        '''

        data_struct=struct();
        data_struct=bts.extract_ts(bk_tag, wvr_scan)
        self.data_struct=data_struct


    def start_up(self, wvr_scan, bk_tag, pf):
        f = open('pointing_parameters_2018_fast.txt','rb')
        point_par = pk.load(f)
        f.close()
        az_offs = point_par['az_offs']

        #x_az=tod.pointing.hor.az
        waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan)#, az_lim=(np.min(x_az), np.max(x_az)))
        time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

        #Extracting Trj atmograms from WVR PWV
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=1, out_pk='Trj_pk_BAK.txt', posting_folder=pf)



    def ps_az(az_ordered_az, az_ordered_data):

        az_ordered_az=az_ordered_az[np.where(np.isfinite(az_ordered_data))]
        az_ordered_data=az_ordered_data[np.where(np.isfinite(az_ordered_data))]

        az_ = az_ordered_az
        daz = az_ordered_az[10] - az_ordered_az[9]
        az = len(az_)
        df=1./daz
        fft_wvr = np.fft.fft(az_ordered_data)/az
        ps    = np.square(np.abs(fft_wvr))*np.hanning(len(az_ordered_data))
        freq_ = np.fft.fftfreq(az,daz)

        idx_pos=np.where(freq_>=0)

        return freq_[idx_pos], ps[idx_pos]


    def plot_atmogram(self, wvr_scan, bk_tag, pf):


        #pf='../../../Postings/WVR_postings/20220210_BK_WVR_correlations/plots/'


        ets=self.data_struct

        #x_az=tod.pointing.hor.az
        waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan, show=0)#, az_lim=(np.min(x_az), np.max(x_az)))
        time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

        #Extracting Trj atmograms from WVR PWV
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
        wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

        #det=21
        #rx_list = [270, 210]
        rx_list = [210, 270] #270 is done

        wvr_atmo = ets.wvr_struct.D_pwv
        az_wvr = ets.wvr_struct.az_real #az_calib
        time_ordered_az=ets.wvr_struct.az_wvr
        fs = ets.wvr_struct.fs
        t_wvr=ets.wvr_struct.tot

        tod = ets.bk_struct
        fs_bk = tod.fs


        for rx in rx_list:

            wvr_atmo_Trx = wvr_atmo_Trj[str(rx)]#-np.nanmean(wvr_atmo_Trj[str(rx)])

            x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')
            x_b = x_pair[np.where(det_pol=='b')]
            x_a = x_pair[np.where(det_pol=='a')]
            y_b = y_pair[np.where(det_pol=='b')]
            y_a = y_pair[np.where(det_pol=='a')]

            for det in range(len(det_a_list)):

                try:

                    a_det = det_a_list[det]
                    i_det_b=np.where(x_b==x_a[det])[0]
                    b_det = det_b_list[i_det_b[0]]

                    pl.figure(figsize=(10,6))
                    pl.scatter(x,y, c='k')
                    pl.scatter(x_a[det],y_a[det], s=150, marker='o', c='r', label='det A')
                    pl.scatter(x_b[i_det_b],y_b[i_det_b], s=150, marker='*', c='y', label='det B')
                    pl.xlabel('Az[deg]')
                    pl.ylabel('El[deg]')
                    pl.suptitle('FPU Map')
                    pl.legend(title='Selected Pair')
                    pl.title('Az offs = '+str(round(x_a[det],2))+' - El Offs = '+str(round(y_a[det],2)))
                    pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(det)+'.png')
                    pl.savefig(pf+bk_tag[:8]+'det_fpu_location_'+str(rx)+'_det'+str(det)+'.png')
                    pl.close()

                    #az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det)

                    corr_list_p3True=[]
                    corr_list_p3False=[]

                    for p3_filt in [False, True]:

                        print('p3_filt = ', p3_filt)

                        az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det, posting_folder='None')

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
                        print(calib_az)
                        if  any(caz <0 for caz in calib_az):
                            print('adding 360')
                            calib_az=[i+360 for i in calib_az]
                        print(calib_az)

                        x_min_wvr_idx = np.where(calib_az>= np.min(x_diff))[0][0]
                        x_max_wvr_idx = np.where(calib_az>= np.max(x_diff))[0][0]


                        print(t_min_wvr, t_max_wvr)
                        print(t_min_bk, t_max_bk)

                        bk_time_mask = np.where((full_time_dt<=t_max_wvr)&(full_time_dt>=t_min_wvr))[0]
                        wvr_time_mask = np.where((full_time_wvr_dt>=t_min_bk)&(full_time_wvr_dt<=t_max_bk))[0]

                        if (t_min_bk <= t_min_wvr):
                            extent_wvr = [0, len(wvr_time_mask), np.min(calib_az), np.max(calib_az)]
                        else:
                            extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]

                        #
                        # full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])
                        # full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr])
                        #
                        # t_min_wvr=datetime.time(full_time_wvr[0])
                        # t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
                        # t_min_bk=datetime.time(full_time[0])
                        # t_max_bk=datetime.time(full_time[len(full_time)-1])
                        #
                        # x_min_wvr_idx = np.where(calib_az>= np.min(x_diff))[0][0]
                        #
                        # x_max_wvr_idx = np.where(calib_az>= np.max(x_diff))[0][0]
                        #
                        # if (t_min_bk <= t_min_wvr):
                        #
                        #     bk_time_mask = np.where(full_time_dt>=t_min_wvr)[0]
                        #     wvr_time_mask = np.where(full_time_wvr_dt<=t_max_bk)[0]
                        #     extent_wvr = [0, len(wvr_time_mask), np.min(calib_az), np.max(calib_az)]
                        #
                        # else:
                        #
                        #     bk_time_mask = np.where(full_time_dt<=t_max_wvr)[0]
                        #     wvr_time_mask = np.where(full_time_wvr_dt>=t_min_bk)[0]
                        #
                        #     extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]
                        #



                        fig,(ax1, ax2, ax3, ax4) = pl.subplots(4,1, figsize=(12,8))

                        xlabels=np.array(xlabels)
                        xlabels_wvr=np.array(xlabels_wvr)

                        xticks_mask=xticks[bk_time_mask]
                        xlabels_mask=xlabels[bk_time_mask]
                        xticks_wvr_mask=xticks_wvr[wvr_time_mask]
                        xlabels_wvr_mask=xlabels_wvr[wvr_time_mask]

                        D_sum_bk_tcut = D_sum_bk[:,bk_time_mask]
                        D_diff_bk_tcut = D_diff_bk[:,bk_time_mask]
                        wvr_atmo_tcut = wvr_atmo[:,wvr_time_mask]
                        wvr_atmo_Trx_tcut = wvr_atmo_Trx[:,wvr_time_mask]

                        pos1=ax1.imshow(D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_sum), np.max(x_sum)], origin='lower')
                        ax1.set_title('BK'+str(rx)+' Pair Sum')
                        cbar1 = pl.colorbar(pos1, ax=ax1)
                        cbar1.set_label('T[K]')
                        ax1.set_ylabel('Az[deg]')
                        ax1.set_xticks(xticks_mask[::10])
                        #ax1.set_xticklabels(xlabels_mask[::10])
                        ax1.set_xticklabels([])
                        #pl.yticks(fontsize=fs_ticks)

                        pos2=ax2.imshow(D_diff_bk_tcut, aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_diff), np.max(x_diff)], origin='lower')
                        ax2.set_title('BK'+str(rx)+' Pair Diff')
                        cbar2 = pl.colorbar(pos2, ax=ax2)
                        cbar2.set_label('T[K]')
                        ax2.set_ylabel('Az[deg]')
                        ax2.set_xticks(xticks_mask[::10])
                        #ax2.set_xticklabels(xlabels_mask[::10])
                        ax2.set_xticklabels([])
                        #pl.yticks(fontsize=fs_ticks)
                        #az_labels=[int(j) for j in calib_az]

                        pos3=ax3.imshow(wvr_atmo_tcut, aspect='auto', interpolation='nearest', extent=extent_wvr, origin='lower')
                        ax3.set_title('PWV')
                        cbar3 = pl.colorbar(pos3, ax=ax3)
                        cbar3.set_label('PWV[um]')
                        ax3.set_ylabel('Az[deg]')
                        ax3.set_ylim(np.min(x_diff), np.max(x_diff))
                        ax3.set_xticks(xticks_wvr_mask[::6])
                        ax3.set_xticklabels([])

                        pos4=ax4.imshow(wvr_atmo_Trx_tcut, aspect='auto', interpolation='nearest', extent=extent_wvr, origin='lower')
                        ax4.set_title('WVR Trj '+str(rx))
                        cbar4 = pl.colorbar(pos4, ax=ax4)
                        cbar4.set_label('Trj[K]')
                        ax4.set_ylabel('Az[deg]')
                        ax4.set_ylim(np.min(x_diff), np.max(x_diff))

                        ax4.set_xticks(xticks_wvr_mask[::6])
                        ax4.set_xticklabels(xlabels_wvr_mask[::6])

                        pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'.png')
                        pl.savefig(pf+'ScanPairDate_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3_'+str(p3_filt)+'.png')
                        pl.close()


                        nscans_bk = len(D_sum_bk[10,:])

                        az_bk_tag = az_bk[int(fs_bk.sf[0]):int(fs_bk.ef[nscans_bk-1])]
                        #print('az_bk_tag=', az_bk_tag)
                        az_bk_shifted = [(az_bk_i + az_offs) for az_bk_i in az_bk_tag]
                        #print('az_bk_shifted=', az_bk_shifted)

                        if  any(caz <0 for caz in az_bk_shifted):
                            print('adding 360')
                            az_bk_shifted=[i+360 for i in az_bk_shifted]
                        if  any(caz >360. for caz in az_bk_shifted):
                            print('subtracting 360')
                            az_bk_shifted=[i-360 for i in az_bk_shifted]


                        nscans=len(wvr_atmo[10,:])

                        tod_wvr = np.zeros(nscans)
                        tod_wvr_list = []
                        tod_wvr_Trx_list = []

                        x_az_atmo=np.arange(len(wvr_atmo[:,10]))
                        #print('x_az_atmo=', x_az_atmo)

                        mask_match_bk = np.where((x_az_atmo>np.nanmin(az_bk_shifted)) & (x_az_atmo<np.nanmax(az_bk_shifted)))
                        wvr_atmo_matchbk = wvr_atmo_tcut[mask_match_bk[0],:]
                        wvr_atmo_Trx_matchbk = wvr_atmo_Trx_tcut[mask_match_bk[0],:]

                        x_az_atmo_matchbk = [az_i-az_offs for az_i in x_az_atmo[mask_match_bk]]
                        x_az_atmo_matchbk_labels = [int(az_i) for az_i in x_az_atmo_matchbk]
                        #x_az_atmo_matchbk = x_az_atmo[mask_match_bk]


                        fig,ax1 = pl.subplots(1,1)
                        pos1=ax1.imshow(wvr_atmo_matchbk, aspect='auto', interpolation='nearest', origin='lower')
                        ax1.set_title('WVR PWV\n'+wvr_scan[:-4])
                        cbar1 = pl.colorbar(pos1, ax=ax1)
                        cbar1.set_label('PWV[um]')
                        ax1.set_ylabel('Az')
                        #ax1.set_xlim(0,int(nscans/2.))
                        ax1.set_yticks(np.arange(len(wvr_atmo_matchbk[:,10]))[::10])
                        ax1.set_yticklabels(x_az_atmo_matchbk_labels[::10])
                        #pl.yticks(fontsize=fs_ticks)

                        pl.savefig(pf+'wvr_'+wvr_scan[:-4]+'_bk_cut.png')
                        pl.savefig(pf+'wvr_'+wvr_scan[:8]+'_bk_cut.png')
                        pl.close()




        def plot_atmogram(self, wvr_scan, bk_tag, pf):


        #pf='../../../Postings/WVR_postings/20220210_BK_WVR_correlations/plots/'


        ets=self.data_struct

        #x_az=tod.pointing.hor.az
        waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan, show=0)#, az_lim=(np.min(x_az), np.max(x_az)))
        time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

        #Extracting Trj atmograms from WVR PWV
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
        wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

        #det=21
        #rx_list = [270, 210]
        rx_list = [210, 270] #270 is done

        wvr_atmo = ets.wvr_struct.D_pwv
        az_wvr = ets.wvr_struct.az_real #az_calib
        time_ordered_az=ets.wvr_struct.az_wvr
        fs = ets.wvr_struct.fs
        t_wvr=ets.wvr_struct.tot

        tod = ets.bk_struct
        fs_bk = tod.fs


        for rx in rx_list:

            wvr_atmo_Trx = wvr_atmo_Trj[str(rx)]#-np.nanmean(wvr_atmo_Trj[str(rx)])

            x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')
            x_b = x_pair[np.where(det_pol=='b')]
            x_a = x_pair[np.where(det_pol=='a')]
            y_b = y_pair[np.where(det_pol=='b')]
            y_a = y_pair[np.where(det_pol=='a')]

            for det in range(len(det_a_list)):

                try:

                    a_det = det_a_list[det]
                    i_det_b=np.where(x_b==x_a[det])[0]
                    b_det = det_b_list[i_det_b[0]]

                    pl.figure(figsize=(10,6))
                    pl.scatter(x,y, c='k')
                    pl.scatter(x_a[det],y_a[det], s=150, marker='o', c='r', label='det A')
                    pl.scatter(x_b[i_det_b],y_b[i_det_b], s=150, marker='*', c='y', label='det B')
                    pl.xlabel('Az[deg]')
                    pl.ylabel('El[deg]')
                    pl.suptitle('FPU Map')
                    pl.legend(title='Selected Pair')
                    pl.title('Az offs = '+str(round(x_a[det],2))+' - El Offs = '+str(round(y_a[det],2)))
                    pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(det)+'.png')
                    pl.savefig(pf+bk_tag[:8]+'det_fpu_location_'+str(rx)+'_det'+str(det)+'.png')
                    pl.close()

                    #az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det)

                    corr_list_p3True=[]
                    corr_list_p3False=[]

                    for p3_filt in [False, True]:

                    nscans_wvr_cut=len(wvr_atmo_matchbk[10,:])

                    for iscan in range(nscans_wvr_cut-1):

                        az_wvr_onescan = az_wvr[fs.s[iscan]:fs.e[iscan]]
                        time_ordered_az_onescan = time_ordered_az[fs.s[iscan]:fs.e[iscan]]

                        scan_tod = wvr_atmo_matchbk[:,iscan]
                        scan_tod_Trx = wvr_atmo_Trx_matchbk[:,iscan]
                        tod_wvr_list.append(scan_tod)
                        tod_wvr_Trx_list.append(scan_tod_Trx)

                    #for BK data

                    nscans_bk = len(D_sum_bk_tcut[10,:])
                    tod_bk_sum_list = []
                    tod_bk_diff_list = []
                    tod_bk_sum_full_list = []
                    tod_bk_diff_full_list = []

                    for iscan in range(nscans_bk-1):
                        az_bk_onescan = az_bk[int(fs_bk.sf[iscan]):int(fs_bk.ef[iscan])]
                        x_az_atmo_bk=np.arange(np.min(az_bk_onescan), np.max(az_bk_onescan))
                        scan_tod_bk_sum = D_sum_bk_tcut[:,iscan]
                        scan_tod_bk_diff = D_diff_bk_tcut[:,iscan]

                        scan_tod_bk_sum_full = D_sum_bk[:,iscan]
                        scan_tod_bk_diff_full = D_diff_bk[:,iscan]

                        tod_bk_sum_list.append(scan_tod_bk_sum)
                        tod_bk_diff_list.append(scan_tod_bk_diff)

                        tod_bk_sum_full_list.append(scan_tod_bk_sum_full)
                        tod_bk_diff_full_list.append(scan_tod_bk_diff_full)



                    #BK sum-diff correlation

                    corr_coeff_list=[]

                    for i in range(len(tod_bk_diff_list)):
                        pl.figure()
                        pl.scatter(tod_bk_diff_full_list[i], tod_bk_sum_full_list[i])
                        corr_coeff = scipy.stats.pearsonr(tod_bk_diff_full_list[i], tod_bk_sum_full_list[i])
                        corr_coeff_list.append(corr_coeff[0])


                        pl.title('r = '+str(corr_coeff))
                        #pl.show()
                        pl.close()

                    xaxis=xticks_mask[:-1]
                    xlabels=xlabels_mask[:-1]

                    pl.figure()

                    fig,(ax1, ax2, ax0) = pl.subplots(3,1, figsize=(12,8))

                    axes = [ax0, ax1, ax2]

                    pos0=ax0.scatter(xaxis, corr_coeff_list, c='r')
                    ax0.plot(xaxis, corr_coeff_list, c='k', alpha=0.5)
                    ax0.set_title('PSum-PDiff Corr Coeff')
                    ax0.set_ylabel('r')
                    ax0.set_xticks(xaxis[::10])
                    ax0.set_xticklabels(xlabels_mask[::10])

                    pos1=ax1.imshow(D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_sum), np.max(x_sum)], origin='lower')
                    ax1.set_title('BK'+str(rx)+' Pair Sum')
                    # cbar1 = pl.colorbar(pos1, ax=ax1)
                    # cbar1.set_label('T[K]')
                    ax1.set_ylabel('Az[deg]')
                    ax1.set_xticks(xticks_mask[::10])
                    #ax1.set_xticklabels(xlabels_mask[::10])
                    ax1.set_xticklabels([])
                    #pl.yticks(fontsize=fs_ticks)

                    pos2=ax2.imshow(D_diff_bk_tcut, aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_diff), np.max(x_diff)], origin='lower')
                    ax2.set_title('BK'+str(rx)+' Pair Diff')
                    cbar2 = pl.colorbar(pos2, ax=axes)
                    cbar2.set_label('T[K]')
                    ax2.set_ylabel('Az[deg]')
                    ax2.set_xticks(xticks_mask[::10])
                    ax2.set_xticklabels([])


                    pl.suptitle(bk_tag+'\nrx = '+str(rx)+', GCP idx a/b = '+str(a_det)+'/'+str(b_det))

                    pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_psum-pdiff_correlation_p3_'+str(p3_filt)+'.png')
                    pl.savefig(pf+'ScanPairDate_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_psum-pdiff_correlation_p3_'+str(p3_filt)+'.png')
                    pl.close()

                    if p3_filt==True:
                        corr_list_p3True.append(np.nanmean(corr_coeff_list))
                    elif p3_filt==False:
                        corr_list_p3False.append(np.nanmean(corr_coeff_list))

                    #WVR sky avg over time

                    sky_avg=[]
                    sky_std=[]

                    for i in range(len(tod_wvr_Trx_list)):
                        sky_avg.append(np.nanmean(tod_wvr_Trx_list[i]))
                        sky_std.append(np.nanstd(tod_wvr_Trx_list[i]))

                    xaxis_wvr = xticks_wvr_mask[:-1]

                    fig,(ax1, ax2, ax3) = pl.subplots(3,1, figsize=(12,8))

                    axes = [ax1, ax2, ax3]

                    pos1=ax1.scatter(xaxis_wvr, sky_avg, c='blue')
                    ax1.plot(xaxis_wvr, sky_avg, c='k', alpha=0.5)
                    ax1.set_xticks(xaxis_wvr[::10])
                    ax1.set_xticklabels([])
                    ax1.set_title('All Sky Average')
                    ax1.set_ylabel('T[K]')

                    pos2=ax2.scatter(xaxis_wvr, sky_std, c='g')
                    ax2.plot(xaxis_wvr, sky_std, c='k', alpha=0.5)
                    ax2.set_xticks(xaxis_wvr[::10])
                    ax2.set_xticklabels([])
                    ax2.set_title('All Sky Std')
                    ax2.set_ylabel('T[K]')

                    pos3=ax3.imshow(wvr_atmo_Trx_tcut, aspect='auto', interpolation='nearest', extent=extent_wvr, origin='lower')
                    ax3.set_title('WVR Trj '+str(rx)+'GHz')
                    cbar3 = pl.colorbar(pos3, ax=axes)
                    cbar3.set_label('Trj[K]')
                    ax3.set_ylabel('Az[deg]')
                    ax3.set_ylim(np.min(x_diff), np.max(x_diff))
                    ax3.set_xticks(xticks_wvr_mask[::6])
                    ax3.set_xticklabels(xlabels_wvr_mask[::6])

                    pl.suptitle(wvr_scan[:-4])
                    pl.savefig(pf+'WVR_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_sky_az_average.png')
                    pl.savefig(pf+'WVR_'+wvr_scan[:8]+'_rx_'+str(rx)+'_sky_az_average.png')
                    pl.close()


                    #deg_power_wvr=[]

                    # for i in range(len(tod_wvr_Trx_list)):
                    #
                    #     az_ordered_data = tod_wvr_Trx_list[i]
                    #     az_ordered_az = x_az_atmo_matchbk_labels
                    #
                    #     pl.figure()
                    #     pl.plot(az_ordered_az, az_ordered_data)
                    #     pl.suptitle('Timestream')
                    #     pl.xlabel('Az[deg]')
                    #     pl.ylabel('T[K]')
                    #     pl.savefig(pf+'WVR_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_timestream_'+str(i)+'.png')
                    #     pl.close()
                    #
                    #     freq, ps = ps_az(np.array(az_ordered_az), np.array(az_ordered_data))
                    #
                    #
                    #     pl.figure()
                    #     pl.plot(freq, ps, c='k', alpha=0.5)
                    #     pl.scatter(freq, ps, c='r')
                    #     pl.xlabel('1/az[deg^-1]')
                    #     pl.ylabel('PS')
                    #     pl.loglog()
                    #     pl.title('Power Spectrum')
                    #     pl.savefig(pf+'WVR_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_powerspectrum_'+str(i)+'.png')
                    #     pl.close()
                    #
                    #
                    # for i in range(len(tod_bk_diff_list)):
                    #
                    #     az_ordered_data = tod_bk_diff_full_list[i][:-1]
                    #     az_ordered_az = np.arange(np.min(x_diff), np.max(x_diff))
                    #
                    #     pl.figure()
                    #     pl.plot(az_ordered_az, az_ordered_data)
                    #     pl.title('ts')
                    #
                    #     pl.figure()
                    #     pl.plot(az_ordered_az, az_ordered_data)
                    #     pl.suptitle('Timestream')
                    #     pl.xlabel('Az[deg]')
                    #     pl.ylabel('T[K]')
                    #     pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_timestream_'+str(i)+'_p3_'+str(p3_filt)+'.png')
                    #     pl.close()
                    #
                    #     freq, ps = ps_az(np.array(az_ordered_az), np.array(az_ordered_data))
                    #
                    #     pl.figure()
                    #     pl.plot(freq, ps, c='k', alpha=0.5)
                    #     pl.scatter(freq, ps, c='r')
                    #     pl.xlabel('1/az[deg^-1]')
                    #     pl.ylabel('PS')
                    #     pl.loglog()
                    #     pl.title('Power Spectrum')
                    #     pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_powerspectrum_'+str(i)+'_p3_'+str(p3_filt)+'.png')
                    #     pl.close()





            except:
                print('Scan '+wvr_scan[:-4])
                print('rx = '+str(rx))
                print('det '+str(det))
                print('Failed.')


        print('corr_list_p3True = ', corr_list_p3True)
        print('corr_list_p3False = ', corr_list_p3False)

        pl.figure()
        n, x, _ = pl.hist(corr_list_p3True, bins=np.linspace(-1, 1, 50), histtype=u'step', density=True)
        pl.suptitle('Psum-Pdiff Correlation Histogram')
        pl.xlabel('r')
        pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_psum-pdiff_correlation_p3_True.png')
        pl.savefig(pf+'ScanPairDate_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_histo_psum-pdiff_correlation_p3_True.png')
        pl.close()

        pl.figure()
        n, x, _ = pl.hist(corr_list_p3False, bins=np.linspace(-1, 1, 50), histtype=u'step', density=True)
        pl.suptitle('Psum-Pdiff Correlation Histogram')
        pl.xlabel('r')
        pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_psum-pdiff_correlation_p3_False.png')
        pl.savefig(pf+'ScanPairDate_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_histo_psum-pdiff_correlation_p3_False.png')
        pl.close()








def BK_WVR_ts_corr(bk_atmo_diff, bk_atmo_sum, wvr_atmo, xlabels_mask, xlabels_wvr_mask, save_fn_list, save_txt):

    idx_s_bk=[0]
    idx_s_wvr=[0]
    idx_e_bk=[]
    idx_e_wvr=[]

    for i in range(1, len(xlabels_wvr_mask)-1):
        #print(xlabels_wvr_mask[i][-2:])
        if xlabels_wvr_mask[i+1][-2:] != xlabels_wvr_mask[i][-2:]:
            idx_s_wvr.append(i+1)
            idx_e_wvr.append(i)

    idx_e_wvr.append(len(xlabels_wvr_mask))


    for i in range(1, len(xlabels_mask)-1):
        #print(xlabels_wvr_mask[i][-2:])
        if xlabels_mask[i+1][-2:] != xlabels_mask[i][-2:]:
            idx_s_bk.append(i+1)
            idx_e_bk.append(i)

    idx_e_bk.append(len(xlabels_mask))

    wvr_atmo_binned_tods=[]
    bk_atmo_sum_binned_tods=[]
    bk_atmo_diff_binned_tods=[]

    print('idx_s_bk=', idx_s_bk)

    for i in range(len(idx_s_bk)):
    #for i in range(1):
        s_wvr=idx_s_wvr[i]
        e_wvr=idx_e_wvr[i]
        s_bk=idx_s_bk[i]
        e_bk=idx_e_bk[i]
        wvr_atmo_tmp=wvr_atmo[:,s_wvr:e_wvr]
        bk_atmo_sum_tmp=bk_atmo_sum[:,s_bk:e_bk]
        bk_atmo_diff_tmp=bk_atmo_diff[:,s_bk:e_bk]
        wvr_tod_binned=[]
        bk_tod_sum_binned=[]
        bk_tod_diff_binned=[]
        for j in range(len(wvr_atmo[:,0])):
            wvr_tod_binned.append(np.nanmean(wvr_atmo_tmp[j,:]))
        for k in range(len(bk_atmo_sum[:,0])):
            bk_tod_sum_binned.append(np.nanmean(bk_atmo_sum_tmp[k,:]))
            bk_tod_diff_binned.append(np.nanmean(bk_atmo_diff_tmp[k,:]))

        wvr_atmo_binned_tods.append(wvr_tod_binned)
        bk_atmo_sum_binned_tods.append(bk_tod_sum_binned)
        bk_atmo_diff_binned_tods.append(bk_tod_diff_binned)
        #now they have same size in t axis

    corr_coeff_psum_wvr_list=[]
    corr_coeff_pdiff_wvr_list=[]
    corr_coeff_psum_pdiff_list=[]

    xaxis_corr=xlabels_mask[idx_s_bk]

    print('bk_atmo_sum_binned_tods=', bk_atmo_sum_binned_tods)

    for i in range(len(bk_atmo_sum_binned_tods)):
    #for i in range(1):
        ts_bk_sum=np.array(bk_atmo_sum_binned_tods[i])
        ts_bk_diff=np.array(bk_atmo_diff_binned_tods[i])
        ts_wvr=np.array(wvr_atmo_binned_tods[i])

        ts_bk_sum[np.isnan(ts_bk_sum)]=0.
        ts_bk_diff[np.isnan(ts_bk_diff)]=0.
        ts_wvr[np.isnan(ts_wvr)]=0

        ts_bk_sum_interp=np.interp(np.arange(100), np.arange(len(ts_bk_sum)), ts_bk_sum)
        ts_bk_diff_interp=np.interp(np.arange(100), np.arange(len(ts_bk_diff)), ts_bk_diff)
        ts_wvr_interp=np.interp(np.arange(100), np.arange(len(ts_wvr)), ts_wvr)

        corr_coeff_psum_wvr = scipy.stats.pearsonr(ts_bk_sum_interp, ts_wvr_interp)
        corr_coeff_psum_wvr_list.append(corr_coeff_psum_wvr[0])

        corr_coeff_pdiff_wvr = scipy.stats.pearsonr(ts_bk_diff_interp, ts_wvr_interp)
        corr_coeff_pdiff_wvr_list.append(corr_coeff_pdiff_wvr[0])

        corr_coeff_psum_pdiff = scipy.stats.pearsonr(ts_bk_sum_interp, ts_bk_diff_interp)
        corr_coeff_psum_pdiff_list.append(corr_coeff_psum_pdiff[0])

    #plot results in 2 different flavours

    xticks=np.arange(len(xaxis_corr))

    pl.figure()
    pl.plot(xaxis_corr, corr_coeff_psum_wvr_list, c='k', alpha=0.5)
    pl.scatter(xaxis_corr, corr_coeff_psum_wvr_list, c='blue', label='psum-wvr')
    pl.plot(xaxis_corr, corr_coeff_pdiff_wvr_list, c='k', alpha=0.5)
    pl.scatter(xaxis_corr, corr_coeff_pdiff_wvr_list, c='red', label='pdiff-wvr')
    pl.xticks(xticks[::3])
    pl.plot(xaxis_corr, corr_coeff_psum_pdiff_list, c='k', alpha=0.5)
    pl.scatter(xaxis_corr, corr_coeff_psum_pdiff_list, c='green', label='psum-pdiff')
    pl.legend()
    pl.ylabel('r')
    pl.suptitle('Pearson Correlation Coefficient')
    for save_fn in save_fn_list:
        pl.savefig(save_fn+'_Pearson_Corr_Coeff_v1.png')
    pl.close()



    fig,(ax1, ax2, ax3) = pl.subplots(3,1, figsize=(8,6))

    pos1=ax1.plot(xaxis_corr, corr_coeff_psum_wvr_list, c='k', alpha=0.5)
    ax1.scatter(xaxis_corr, corr_coeff_psum_wvr_list, c='blue', label='psum-wvr')
    ax1.set_xticks(xticks[::3])
    ax1.set_xticklabels([])
    ax1.set_title('Psum - WVR Trj')
    ax1.set_ylabel('r')

    pos2=ax2.plot(xaxis_corr, corr_coeff_pdiff_wvr_list, c='k', alpha=0.5)
    ax2.scatter(xaxis_corr, corr_coeff_pdiff_wvr_list, c='red', label='pdiff-wvr')
    ax2.set_xticks(xticks[::3])
    ax2.set_xticklabels([])
    ax2.set_title('Pdiff - WVR Trj')
    ax2.set_ylabel('r')

    pos3=ax3.plot(xaxis_corr, corr_coeff_psum_pdiff_list, c='k', alpha=0.5)
    ax3.scatter(xaxis_corr, corr_coeff_psum_pdiff_list, c='green', label='psum-pdiff')
    ax3.set_title('Psum - Pdiff')
    ax3.set_ylabel('r')
    ax3.set_xticks(xticks[::3])

    pl.suptitle('Pearson Correlation Coefficient')
    for save_fn in save_fn_list:
        pl.savefig(save_fn+'_Pearson_Corr_Coeff_v2.png')
    pl.close()

    Pearson_Corr_Coeff = {'psum-wvr':corr_coeff_psum_wvr_list, 'pdiff-wvr':corr_coeff_pdiff_wvr_list, 'psum-pdiff':corr_coeff_psum_pdiff_list, 'time_axis':xaxis_corr}

    f = open(save_txt,'wb')
    pk.dump(Pearson_Corr_Coeff, f)
    f.close()









def corr_map_fpu(bk_tag, wvr_scan, pf):

    def r_fpu(bk_tag, wvr_scan, rx, p3_filt, nscan, pf):

        x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')

        x_b = x_pair[np.where(det_pol=='b')]
        x_a = x_pair[np.where(det_pol=='a')]
        y_b = y_pair[np.where(det_pol=='b')]
        y_a = y_pair[np.where(det_pol=='a')]

        psum_wvr = np.full(len(det_a_list), np.nan)
        pdiff_wvr = np.full(len(det_a_list), np.nan)
        psum_pdiff = np.full(len(det_a_list), np.nan)

        for det in range(len(det_a_list)):

            a_det = det_a_list[det]
            i_det_b=np.where(x_b==x_a[det])[0]
            b_det = det_b_list[i_det_b[0]]

            fn='BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_pk.txt'

            f = open(fn,'rb')
            corr = pk.load(f)
            f.close()


            psum_wvr[det]=corr['psum-wvr'][nscan]
            pdiff_wvr[det]=corr['pdiff-wvr'][nscan]
            psum_pdiff[det]=corr['psum-pdiff'][nscan]

            scan_time = corr['time_axis'][nscan]


        return psum_wvr, pdiff_wvr, psum_pdiff, scan_time, x, y, det_a_list, det_b_list, x_a, y_a, x_b, y_b




    ets=bts.extract_ts(bk_tag, wvr_scan)

    waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan)#, az_lim=(np.min(x_az), np.max(x_az)))
    time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

    tod=ets.bk_struct
    fs_bk=tod.fs
    #Extracting Trj atmograms from WVR PWV
    D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

    for rx in [270, 210]:
        wvr_atmo_Trx = wvr_atmo_Trj[str(rx)]#-np.nanmean(wvr_atmo_Trj[str(rx)])

        for p3_filt in [True, False]:

            fn_trial='BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(0)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_pk.txt'

            f = open(fn_trial,'rb')
            corr_trial = pk.load(f)
            f.close()

            nscans = len(corr_trial['time_axis'])
            psum_wvr, pdiff_wvr, psum_pdiff, scan_time, x, y, det_a_list, det_b_list, x_a, y_a, x_b, y_b = r_fpu(bk_tag, wvr_scan, rx, p3_filt, 0, pf)

            psum_wvr_tlist = []
            pdiff_wvr_tlist = []
            psum_pdiff_tlist = []

            psum_wvr_tlist = np.full((len(det_a_list), nscans), np.nan)
            pdiff_wvr_tlist = np.full((len(det_a_list), nscans), np.nan)
            psum_pdiff_tlist = np.full((len(det_a_list), nscans), np.nan)

            for iscan in range(nscans):

                psum_wvr, pdiff_wvr, psum_pdiff, scan_time, x, y, det_a_list, det_b_list, x_a, y_a, x_b, y_b = r_fpu(bk_tag, wvr_scan, rx, p3_filt, iscan, pf)

                psum_wvr_tlist[:,iscan] = psum_wvr
                pdiff_wvr_tlist[:,iscan] = pdiff_wvr
                psum_pdiff_tlist[:,iscan] = psum_pdiff


                pdiff_wvr_abs = [np.abs(j) for j in pdiff_wvr]
                psum_pdiff_abs = [np.abs(j) for j in psum_pdiff]


                time_stamp = str(scan_time[:2])+'_'+str(scan_time[3:])
                #
                # print(bk_tag+'_psum-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')

                pl.figure(figsize=(12,8))
                pl.scatter(x,y, c='k')
                pl.scatter(x_a,y_a, s=250, marker='o', c=psum_wvr, label='Psum-WVRTrj')
                pl.xlabel('Az[deg]')
                pl.ylabel('El[deg]')
                pl.legend(title='Correlation Type')
                pl.suptitle(bk_tag+'\n'+time_stamp)
                pl.title('Rx'+str(rx))
                cbar = pl.colorbar()
                cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
                pl.savefig(pf+bk_tag+'_psum-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
                pl.savefig(pf+bk_tag[:8]+'_psum-wvr_FPU_'+str(rx)+'_nscan'+str(iscan)+'_p3_'+str(p3_filt)+'.png')
                pl.close()

                pl.figure(figsize=(12,8))
                pl.scatter(x,y, c='k')
                pl.scatter(x_a,y_a, s=250, marker='o', c=pdiff_wvr, label='Pdiff-WVRTrj')
                pl.xlabel('Az[deg]')
                pl.ylabel('El[deg]')
                pl.legend(title='Correlation Type')
                pl.suptitle(bk_tag+'\n'+time_stamp)
                pl.title('Rx'+str(rx))
                cbar = pl.colorbar()
                cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
                pl.savefig(pf+bk_tag+'_pdiff-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
                pl.savefig(pf+bk_tag[:8]+'_pdiff-wvr_FPU_'+str(rx)+'_nscan'+str(iscan)+'_p3_'+str(p3_filt)+'.png')
                pl.close()

                pl.figure(figsize=(12,8))
                pl.scatter(x,y, c='k')
                pl.scatter(x_a,y_a, s=250, marker='o', c=psum_pdiff, label='Psum-Pdiff')
                pl.xlabel('Az[deg]')
                pl.ylabel('El[deg]')
                pl.legend(title='Correlation Type')
                pl.suptitle(bk_tag+'\n'+time_stamp)
                pl.title('Rx'+str(rx))
                cbar = pl.colorbar()
                cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
                pl.savefig(pf+bk_tag+'_psum-pdiff_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
                pl.savefig(pf+bk_tag[:8]+'_psum-pdiff_FPU_'+str(rx)+'_nscan'+str(iscan)+'_p3_'+str(p3_filt)+'.png')
                pl.close()


                pl.figure(figsize=(12,8))
                pl.scatter(x,y, c='k')
                pl.scatter(x_a,y_a, s=250, marker='o', c=pdiff_wvr_abs, label='Pdiff-WVRTrj_abs')
                pl.xlabel('Az[deg]')
                pl.ylabel('El[deg]')
                pl.legend(title='Correlation Type')
                pl.suptitle(bk_tag+'\n'+time_stamp)
                pl.title('Rx'+str(rx))
                cbar = pl.colorbar()
                cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
                pl.savefig(pf+bk_tag+'_pdiff-wvr_abs_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
                pl.savefig(pf+bk_tag[:8]+'_pdiff-wvr_abs_FPU_'+str(rx)+'_nscan'+str(iscan)+'_p3_'+str(p3_filt)+'.png')
                pl.close()

                pl.figure(figsize=(12,8))
                pl.scatter(x,y, c='k')
                pl.scatter(x_a,y_a, s=250, marker='o', c=psum_pdiff_abs, label='Psum-Pdiff_abs')
                pl.xlabel('Az[deg]')
                pl.ylabel('El[deg]')
                pl.legend(title='Correlation Type')
                pl.suptitle(bk_tag+'\n'+time_stamp)
                pl.title('Rx'+str(rx))
                cbar = pl.colorbar()
                cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
                pl.savefig(pf+bk_tag+'_psum-pdiff_abs_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
                pl.savefig(pf+bk_tag[:8]+'_psum-pdiff_abs_FPU_'+str(rx)+'_nscan'+str(iscan)+'_p3_'+str(p3_filt)+'.png')
                pl.close()


            #
            # pl.figure()
            # n, x, _ = pl.hist(psum_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'stepfilled', facecolor='g', alpha=0.6)
            # n, x, _ = pl.hist(psum_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'step', color='g')
            # pl.suptitle('Psum-Pdiff Correlation Histogram')
            # pl.xlabel('r')
            # pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_psum-wvr.png')
            # pl.show()
            #
            # pl.figure()
            # n, x, _ = pl.hist(pdiff_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'stepfilled', facecolor='g', alpha=0.6)
            # n, x, _ = pl.hist(pdiff_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'step', color='g')
            # pl.suptitle('Psum-Pdiff Correlation Histogram')
            # pl.xlabel('r')
            # pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_pdiff-wvr.png')
            # pl.show()
            #
            # pl.figure()
            # n, x, _ = pl.hist(psum_pdiff_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'stepfilled', facecolor='g', alpha=0.6)
            # n, x, _ = pl.hist(psum_pdiff_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'step', color='g')
            # pl.suptitle('Psum-Pdiff Correlation Histogram')
            # pl.xlabel('r')
            # pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_psum-pdiff.png')
            # pl.show()
            #

            psum_wvr_tavg = np.full(len(det_a_list), np.nan)
            pdiff_wvr_tavg = np.full(len(det_a_list), np.nan)
            psum_pdiff_tavg = np.full(len(det_a_list), np.nan)

            psum_wvr_tstd = np.full(len(det_a_list), np.nan)
            pdiff_wvr_tstd = np.full(len(det_a_list), np.nan)
            psum_pdiff_tstd = np.full(len(det_a_list), np.nan)

            for det_i in range(len(det_a_list)):

                psum_wvr_tavg[det_i] = np.nanmean(psum_wvr_tlist[det_i,:])
                pdiff_wvr_tavg[det_i] = np.nanmean(pdiff_wvr_tlist[det_i,:])
                psum_pdiff_tavg[det_i] = np.nanmean(psum_pdiff_tlist[det_i,:])

                psum_wvr_tstd[det_i] = np.nanstd(psum_wvr_tlist[det_i,:])
                pdiff_wvr_tstd[det_i] = np.nanstd(pdiff_wvr_tlist[det_i,:])
                psum_pdiff_tstd[det_i] = np.nanstd(psum_pdiff_tlist[det_i,:])



            time_stamp = 't_avg'

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=psum_wvr_tavg, label='Psum-WVRTrj')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_psum-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_psum-wvr_FPU_'+str(rx)+'_nscanavg_p3_'+str(p3_filt)+'.png')
            pl.close()

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=pdiff_wvr_tavg, label='Pdiff-WVRTrj')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_pdiff-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_pdiff-wvr_FPU_'+str(rx)+'_nscanavg_p3_'+str(p3_filt)+'.png')
            pl.close()

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=psum_pdiff_tavg, label='Psum-Pdiff')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_psum-pdiff_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_psum-pdiff_FPU_'+str(rx)+'_nscanavg_p3_'+str(p3_filt)+'.png')
            pl.close()


            pdiff_wvr_tavg_abs = [np.abs(j) for j in pdiff_wvr_tavg]
            psum_pdiff_tavg_abs = [np.abs(j) for j in psum_pdiff_tavg]

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=pdiff_wvr_tavg_abs, label='Pdiff-WVRTrj')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_pdiff-wvr_abs_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_pdiff-wvr_abs_FPU_'+str(rx)+'_nscanavg_p3_'+str(p3_filt)+'.png')
            pl.close()

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=psum_pdiff_tavg_abs, label='Psum-Pdiff')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_psum-pdiff_abs_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_psum-pdiff_abs_FPU_'+str(rx)+'_nscanavg_p3_'+str(p3_filt)+'.png')
            pl.close()

            #
            time_stamp = 't_std'
            #
            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=psum_wvr_tstd, label='Psum-WVRTrj')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_psum-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_psum-wvr_FPU_'+str(rx)+'_nscanstd_p3_'+str(p3_filt)+'.png')
            pl.close()

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=pdiff_wvr_tstd, label='Pdiff-WVRTrj')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_pdiff-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_pdiff-wvr_FPU_'+str(rx)+'_nscanstd_p3_'+str(p3_filt)+'.png')
            pl.close()

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=psum_pdiff_tstd, label='Psum-Pdiff')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_psum-pdiff_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_psum-pdiff_FPU_'+str(rx)+'_nscanstd_p3_'+str(p3_filt)+'.png')
            pl.close()

            pdiff_wvr_tstd_abs = [np.abs(j) for j in pdiff_wvr_tstd]
            psum_pdiff_tstd_abs = [np.abs(j) for j in psum_pdiff_tstd]

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=pdiff_wvr_tstd_abs, label='Pdiff-WVRTrj')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_pdiff-wvr_abs_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_pdiff-wvr_abs_FPU_'+str(rx)+'_nscanstd_p3_'+str(p3_filt)+'.png')
            pl.close()

            pl.figure(figsize=(12,8))
            pl.scatter(x,y, c='k')
            pl.scatter(x_a,y_a, s=250, marker='o', c=psum_pdiff_tstd_abs, label='Psum-Pdiff')
            pl.xlabel('Az[deg]')
            pl.ylabel('El[deg]')
            pl.legend(title='Correlation Type')
            pl.suptitle(bk_tag+'\n'+time_stamp)
            pl.title('Rx'+str(rx))
            cbar = pl.colorbar()
            cbar.ax.set_ylabel('Pearson Correlation Coefficient [r]')
            pl.savefig(pf+bk_tag+'_psum-pdiff_abs_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')
            pl.savefig(pf+bk_tag[:8]+'_psum-pdiff_abs_FPU_'+str(rx)+'_nscanstd_p3_'+str(p3_filt)+'.png')
            pl.close()




def correlation_matrix(wvr_scan, bk_tag):

    pf='../../../Postings/WVR_postings/20220210_BK_WVR_correlations/plots/'


    ets=bts.extract_ts(bk_tag, wvr_scan)

    #x_az=tod.pointing.hor.az
    waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan)#, az_lim=(np.min(x_az), np.max(x_az)))
    time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

    #Extracting Trj atmograms from WVR PWV
    D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

    #det=21
    rx_list = [270, 210]


    wvr_atmo = ets.wvr_struct.D_pwv
    az_wvr = ets.wvr_struct.az_real #az_calib
    time_ordered_az=ets.wvr_struct.az_wvr
    fs = ets.wvr_struct.fs
    t_wvr=ets.wvr_struct.tot

    tod = ets.bk_struct
    fs_bk = tod.fs


    for rx in rx_list:

        wvr_atmo_Trx = wvr_atmo_Trj[str(rx)]#-np.nanmean(wvr_atmo_Trj[str(rx)])

        x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')
        x_b = x_pair[np.where(det_pol=='b')]
        x_a = x_pair[np.where(det_pol=='a')]
        y_b = y_pair[np.where(det_pol=='b')]
        y_a = y_pair[np.where(det_pol=='a')]

        for det in range(len(det_a_list)):

            try:

                a_det = det_a_list[det]
                i_det_b=np.where(x_b==x_a[det])[0]
                b_det = det_b_list[i_det_b[0]]

                pl.figure(figsize=(10,6))
                pl.scatter(x,y, c='k')
                pl.scatter(x_a[det],y_a[det], s=150, marker='o', c='r', label='det A')
                pl.scatter(x_b[i_det_b],y_b[i_det_b], s=150, marker='*', c='y', label='det B')
                pl.xlabel('Az[deg]')
                pl.ylabel('El[deg]')
                pl.suptitle('FPU Map')
                pl.legend(title='Selected Pair')
                pl.title('Az offs = '+str(round(x_a[det],2))+' - El Offs = '+str(round(y_a[det],2)))
                pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(det)+'.png')
                pl.savefig(pf+bk_tag[:8]+'det_fpu_location_'+str(rx)+'_det'+str(det)+'.png')
                pl.close()

                #az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det)

                for p3_filt in [False, True]:

                    print('p3_filt = ', p3_filt)

                    az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det, posting_folder='None')

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

                    # full_time=[t_bk[int(fs_s_i)] for fs_s_i in fs_bk.sf]
                    # xlabels=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time])
                    # nscans_bk = len(fs_bk.sf)
                    # xticks = np.arange(nscans_bk)
                    #
                    # full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])
                    # full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr])
                    #
                    # t_min_wvr=datetime.time(full_time_wvr[0])
                    # t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
                    # t_min_bk=datetime.time(full_time[0])
                    # t_max_bk=datetime.time(full_time[len(full_time)-1])
                    #
                    # x_min_wvr_idx = np.where(calib_az>= np.min(x_diff))[0][0]
                    #
                    # x_max_wvr_idx = np.where(calib_az>= np.max(x_diff))[0][0]
                    #
                    # if (t_min_bk <= t_min_wvr):
                    #
                    #     bk_time_mask = np.where(full_time_dt>=t_min_wvr)[0]
                    #     wvr_time_mask = np.where(full_time_wvr_dt<=t_max_bk)[0]
                    #     extent_wvr = [0, len(wvr_time_mask), np.min(calib_az), np.max(calib_az)]
                    #
                    # else:
                    #
                    #     bk_time_mask = np.where(full_time_dt<=t_max_wvr)[0]
                    #     wvr_time_mask = np.where(full_time_wvr_dt>=t_min_bk)[0]
                    #
                    #     extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]



                    full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])
                    full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr])

                    t_min_wvr=datetime.time(full_time_wvr[0])
                    t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
                    t_min_bk=datetime.time(full_time[0])
                    t_max_bk=datetime.time(full_time[len(full_time)-1])

                    x_min_wvr_idx = np.where(calib_az>= np.min(x_diff))[0][0]
                    x_max_wvr_idx = np.where(calib_az>= np.max(x_diff))[0][0]





                    print(t_min_wvr, t_max_wvr)
                    print(t_min_bk, t_max_bk)

                    bk_time_mask = np.where((full_time_dt<=t_max_wvr)&(full_time_dt>=t_min_wvr))[0]
                    wvr_time_mask = np.where((full_time_wvr_dt>=t_min_bk)&(full_time_wvr_dt<=t_max_bk))[0]

                    if (t_min_bk <= t_min_wvr):
                        extent_wvr = [0, len(wvr_time_mask), np.min(calib_az), np.max(calib_az)]
                    else:
                        extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]



                    xlabels=np.array(xlabels)
                    xlabels_wvr=np.array(xlabels_wvr)

                    xticks_mask=xticks[bk_time_mask]
                    xlabels_mask=xlabels[bk_time_mask]
                    xticks_wvr_mask=xticks_wvr[wvr_time_mask]
                    xlabels_wvr_mask=xlabels_wvr[wvr_time_mask]

                    D_sum_bk_tcut = D_sum_bk[:,bk_time_mask]
                    D_diff_bk_tcut = D_diff_bk[:,bk_time_mask]
                    wvr_atmo_tcut = wvr_atmo[:,wvr_time_mask]
                    wvr_atmo_Trx_tcut = wvr_atmo_Trx[:,wvr_time_mask]

                    nscans_bk = len(D_sum_bk[10,:])

                    az_bk_tag = az_bk[int(fs_bk.sf[0]):int(fs_bk.ef[nscans_bk-1])]
                    az_bk_shifted = [(az_bk_i + az_offs) for az_bk_i in az_bk_tag]

                    nscans=len(wvr_atmo[10,:])

                    tod_wvr = np.zeros(nscans)
                    tod_wvr_list = []
                    tod_wvr_Trx_list = []

                    x_az_atmo=np.arange(len(wvr_atmo[:,10]))

                    mask_match_bk = np.where((x_az_atmo>np.nanmin(az_bk_shifted)) & (x_az_atmo<np.nanmax(az_bk_shifted)))
                    wvr_atmo_matchbk = wvr_atmo_tcut[mask_match_bk[0],:]
                    wvr_atmo_Trx_matchbk = wvr_atmo_Trx_tcut[mask_match_bk[0],:]

                    x_az_atmo_matchbk = [az_i-az_offs for az_i in x_az_atmo[mask_match_bk]]
                    x_az_atmo_matchbk_labels = [int(az_i) for az_i in x_az_atmo_matchbk]
                    #x_az_atmo_matchbk = x_az_atmo[mask_match_bk]



                    #BK-WVR correlation
                    save_fn=pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)
                    save_fn2=pf+'ScanPairDate_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3_'+str(p3_filt)
                    pk_fn='BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_pk.txt'

                    wvr_atmo_T = wvr_atmo_Trx_matchbk - np.nanmean(wvr_atmo_Trx_matchbk)

                    BK_WVR_ts_corr(D_diff_bk_tcut, D_sum_bk_tcut, wvr_atmo_T, xlabels_mask, xlabels_wvr_mask, [save_fn, save_fn2],  pk_fn)
            #
            #
            except:
                print('Scan '+wvr_scan[:-4])
                print('rx = '+str(rx))
                print('det '+str(det))
                print('Failed.')









#main
f = open('pointing_parameters_2018_fast.txt','rb')
point_par = pk.load(f)
f.close()
az_offs = point_par['az_offs']

wvr_scan1='20200415_140134_scanAz_fast.txt'
bk_tag1='20200415D03_dk293'

wvr_scan2='20200418_190135_scanAz_fast.txt'
bk_tag2='20200418B01_dk203'

wvr_scan3='20200723_190134_scanAz_fast.txt'
bk_tag3='20200723B09_dk203'

wvr_scan4='20200921_150134_scanAz_fast.txt'
bk_tag4='20200921B09_dk293'

wvr_scan_list=[wvr_scan3, wvr_scan4]
bk_tag_list=[bk_tag3, bk_tag4]

pf='../../../Postings/WVR_postings/20220210_BK_WVR_correlations/plots/'

#start_up(wvr_scan3, bk_tag3, pf)

for i_scan in range(2):
    wvr_scan=wvr_scan_list[i_scan]
    bk_tag=bk_tag_list[i_scan]

    first_plots(wvr_scan, bk_tag, pf)
    #correlation_matrix(wvr_scan, bk_tag)
    #corr_map_fpu(bk_tag, wvr_scan, pf)
