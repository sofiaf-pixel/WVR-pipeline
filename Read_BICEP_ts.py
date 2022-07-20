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
from scipy.interpolate import interp1d
from scipy import interpolate
from dateutil import parser


raz=pA.ReadAzscan()

class struct(object):
    pass

class toStruct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


class extract_ts(object):


    def __init__(self, bk_tag, wvr_scan, unit=None, verb=True):

        '''

        '''
        #initialize.__init__(self, unit, verb=verb)

        wvr_struct=struct();
        bk_struct=struct();

        bk_struct = self.load_tod(bk_tag)

        print(bk_struct)

        wvr_struct.waz, wvr_struct.az_wvr, wvr_struct.az_real, wvr_struct.fs, wvr_struct.idx, wvr_struct.pwv_full_ts, wvr_struct.D_pwv, wvr_struct.tot = raz.read_pwvatmo(wvr_scan, show=0)

        self.wvr_struct = wvr_struct
        self.bk_struct =  bk_struct
        self.wvr_scan = wvr_scan
        self.bk_tag = bk_tag

        self.pf = '../../../Postings/WVR_postings/20220210_BK_WVR_correlations/plots/'

        #
        # f = open('bk_vs_wvr_ts.txt','wb')
        # pk.dump(data, f)
        # f.close()
        #



    def load_tod(self, tag):

        if not os.path.exists(f'tod/ba/'+tag):
            os.makedirs(f'tod/ba/'+tag)

        d_dict = mat73.loadmat('tod/ba/bicep_array_'+tag+'.mat')

        d = toStruct(**d_dict)


        class struct:
            pass

        tod=struct()

        #print(d.keys())

        tod.data_pairs=toStruct(**d.data_pairs)
        tod.data_pairs_diff=toStruct(**d.data_pairs_diff)
        tod.data_pairs_diff_noavg=toStruct(**d.data_pairs_diff_noavg)
        tod.data_pairs_diff_p0 = toStruct(**d.data_pairs_diff_p0)
        tod.data_pairs_diff_p0_p3 = toStruct(**d.data_pairs_diff_p0_p3)
        tod.data_pairs_diff_noavg_p3 = toStruct(**d.data_pairs_diff_noavg_p3)
        tod.data_pairs_diff_noavg_p3_gs = toStruct(**d.data_pairs_diff_noavg_p3_gs)
        tod.data_pairs_diff_p0_p3_gs = toStruct(**d.data_pairs_diff_p0_p3_gs)
        tod.data_pairs_diff_p0_gs = toStruct(**d.data_pairs_diff_p0_gs)
        tod.data_pairs_diff_noavg_gs = toStruct(**d.data_pairs_diff_noavg_gs)

        tod.ind = toStruct(**d.ind)
        tod.ukpervolt=d.cal
        tod.t_mjd=d.t #modified julian date
        tod.p=toStruct(**d.p)
        tod.pointing=toStruct(**d.pointing)
        tod.fs=toStruct(**d.fs)
        tod.mapind=d.mapind
        tod.det_offs=toStruct(**d.det_offs)


        standard_time=[]
        for mjd_day in tod.t_mjd:
            # MJD=0 corresponds to 1858-Nov-17:00:00:00.
            date = datetime.datetime(1858, 11, 17) + datetime.timedelta(mjd_day)
            standard_time.append(date)

        tod.std=standard_time

        # f=open(f'tod/ba/BA_'+tag+'.pk',"wb")
        # pk.dump(tod, f)
        # f.close()

        return tod


    def find_rgl_idx(self, rx, ind):
        #rx is the virtual receiver
        #not yet implemented for dual band

        if rx==30:
            a=[int(i) for i in ind.rgl030a] #sum
            b=[int(j) for j in ind.rgl030b] #diff
        elif rx==40:
            a=[int(i) for i in ind.rgl040a] #sum
            b=[int(j) for j in ind.rgl040b] #diff
        elif rx==210:
            a=[int(i) for i in ind.rgl210a] #sum
            b=[int(j) for j in ind.rgl210b] #diff
        elif rx==270:
            a=[int(i) for i in ind.rgl270a] #sum
            b=[int(j) for j in ind.rgl270b] #diff
        elif rx==0: #for all receivers
            a=[int(i) for i in ind.a] #sum
            b=[int(j) for j in ind.b] #diff

        a_shift=[x-1 for x in a]
        b_shift=[x-1 for x in b]

        return  a_shift, b_shift


    def calib_tod_rx(self, ts, rx):

        #ts=data['fb']
        #ts_p3=data_p3['fb']

        # ts=tod.data.fb
        # ts_p3=tod.data_p3.fb

        a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)

        ind_rx_a=[int(i) for i in tod.p.rx[a_shift]]
        ind_rx_b=[int(i) for i in tod.p.rx[b_shift]]

        calib_a=tod.ukpervolt[ind_rx_a]*(10**(-6)) #the content of a/b has to be shifted back by one because the rgb index array has not been imported directly to python
        calib_b=tod.ukpervolt[ind_rx_b]*(10**(-6))

        class struct:
            pass

        ts_cal=struct()
        # ts_cal_p3=struct()

        ts_cal.psum = np.multiply(ts[:, a_shift], calib_a)
        ts_cal.pdiff = np.multiply(ts[:, b_shift], calib_b)

        # ts_cal_p3.psum = np.multiply(ts_p3[:, a_shift], calib_a)
        # ts_cal_p3.pdiff = np.multiply(ts_p3[:, b_shift], calib_b)


        return ts_cal



    def calib_tod_rx_onets(self, tod, ts, rx, idx_a, idx_b):

        # ts=tod.data.fb
        # ts_p3=tod.data_p3.fb

        a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)

        ind_rx_a=int(tod.p.rx[idx_a])
        ind_rx_b=int(tod.p.rx[idx_b])

        calib_a=tod.ukpervolt[ind_rx_a]*(10**(-6)) #the content of a/b has to be shifted back by one because the rgb index array has not been imported directly to python
        calib_b=tod.ukpervolt[ind_rx_b]*(10**(-6))

        class struct:
            pass

        ts_cal=struct()
        # ts_cal_p3=struct()

        ts_cal.psum = np.multiply(ts[:, idx_a], calib_a)
        ts_cal.pdiff = np.multiply(ts[:, idx_b], calib_b)
        #
        # ts_cal_p3.psum = np.multiply(ts_p3[:, idx_a], calib_a)
        # ts_cal_p3.pdiff = np.multiply(ts_p3[:, idx_b], calib_b)


        return ts_cal



    def pl_tod_all(self, tag, rx, x_axis='t', npairs=6, posting_folder='None'):

        tod = self.load_tod(tag)
        ts_cal = self.calib_tod_rx(tod.data_pairs_diff, rx)
        ts_cal_p3 = self.calib_tod_rx(tod.data_pairs_diff_p3, rx)
        a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)
        t=tod.std
        x_az=tod.pointing.hor['az']
        fs=tod.fs


        fig = pl.figure(figsize=(14,10))
        fs=14
        ax1 = pl.subplot(211) #to plot Pair Diff from all pairs
        for i in range (npairs):
            if x_axis == 't':
                ax1.scatter(t, ts_cal.pdiff[:,i], s=3, label='GCP_idx='+str(b_shift[i]+1))
            elif x_axis == 'az':
                ax1.scatter(x_az, ts_cal.pdiff[:,i], s=3, label='GCP_idx='+str(b_shift[i]+1))

            ax1.set_title('Pair Diff', fontsize=fs)
            ax1.set_ylabel('T[K]', fontsize=fs)
            pl.legend(loc='lower right')
            #ax1.set_xlim([datetime.datetime(2020, 4, 18, 18, 40, 00), datetime.datetime(2020, 4, 18, 19, 30, 00)])
            #ax1.set_xlim([datetime.datetime(2020, 4, 9, 22, 35, 00), datetime.datetime(2020, 4, 9, 23, 20, 00)])
            ax1.tick_params(axis='both', labelsize=fs)

        ax2 = pl.subplot(212) #to plot Pair Sum from all pairs
        for i in range (npairs):
            if x_axis == 't':
                ax2.scatter(t, ts_cal.psum[:,i], s=3, label='GCP_idx='+str(a_shift[i]+1))
            elif x_axis == 'az':
                ax2.scatter(x_az, ts_cal.psum[:,i], s=3, label='GCP_idx='+str(a_shift[i]+1))
            ax2.set_title('Pair Sum', fontsize=fs)
            ax2.set_ylabel('T[K]', fontsize=fs)
            #ax2.set_xlim([datetime.datetime(2020, 4, 18, 18, 40, 00), datetime.datetime(2020, 4, 18, 19, 30, 00)])
            #ax2.set_xlim([datetime.datetime(2020, 4, 9, 22, 35, 00), datetime.datetime(2020, 4, 9, 23, 20, 00)])
            ax2.tick_params(axis='both', labelsize=fs)

        pl.suptitle(tag, fontsize=16)
        pl.savefig('tod/ba/'+tag+'/rx_'+str(rx)+'_psum-pdiff_all.png')
        if not posting_folder=='None':
            pl.savefig(posting_folder+'/'+tag+'_rx_'+str(rx)+'_psum-pdiff_all.png')
        #pl.show()
        pl.close()



    def det_fpu_location(self, rx, fn_save, el_lim=1.5, show_plots=0):

        tod=self.bk_struct

        x = tod.det_offs.x[np.where(tod.p.band == rx)]
        y = tod.det_offs.y[np.where(tod.p.band == rx)]
        idx_mat = tod.ind.e[np.where(tod.p.band == rx)]
        idx_a_mat = np.array(tod.ind.rgla)#[freq_mask])
        idx_b_mat = np.array(tod.ind.rglb)#[freq_mask])

        idx=np.array([x-1 for x in idx_mat]) #shifting from matlab to python indeces [-1]
        idx_a=np.array([x-1 for x in idx_a_mat])
        idx_b=np.array([x-1 for x in idx_b_mat])

        bk_el = tod.pointing.hor['el']
        bk_el = bk_el[np.where(bk_el>0)]

        El_fpu_center = np.nanmean(bk_el)

        y = np.array(y)
        el_lim_up = (55-El_fpu_center)+el_lim
        el_lim_dn = (55-El_fpu_center)-el_lim
        #picking det at el within 1 FWHM from the center el (55deg)
        mask_zeroel = np.where((y>=el_lim_dn) & (y<=el_lim_up))

        idx_a_masked=[]
        idx_b_masked=[]
        x_pair=[]
        y_pair=[]
        offset_from_bkel = []
        det_pol=[]

        for gcp_idx in idx[mask_zeroel]:
            if gcp_idx in idx_a:
                idx_a_masked.append(int(gcp_idx))
                x_pair.append(x[np.where(idx==gcp_idx)][0])
                y_det = y[np.where(idx==gcp_idx)][0]
                y_pair.append(y_det)
                offset_from_bkel.append(y_det-(55-El_fpu_center))
                det_pol.append('a')
            elif gcp_idx in idx_b:
                idx_b_masked.append(int(gcp_idx))
                x_pair.append(x[np.where(idx==gcp_idx)][0])
                y_det = y[np.where(idx==gcp_idx)][0]
                y_pair.append(y_det)
                offset_from_bkel.append(y_det-(55-El_fpu_center))
                det_pol.append('b')
            else:
                print('not pair idx:', gcp_idx)

        print('idx_a_masked = ', np.shape(idx_a_masked))
        print('idx_b_masked = ', np.shape(idx_b_masked))
        print(np.shape(x_pair))

        pl.figure(figsize=(10,6))
        pl.scatter(x,y)
        pl.scatter(x_pair,y_pair, s=100, marker='*', c=offset_from_bkel)
        pl.title('FPU Map')
        pl.colorbar()
        pl.savefig(fn_save)
        if show_plots==0:
            pl.close()
        else:
            pl.show()
        return x, y, np.array(idx_a_masked), np.array(idx_b_masked), np.array(x_pair), np.array(y_pair), np.array(det_pol), El_fpu_center




    def pl_tod_pair(self, tag, rx, x_axis='t', npairs=100, posting_folder='None'):

        tod = self.load_tod(tag)
        ts_cal = self.calib_tod_rx(tod.data_pairs_diff, rx)
        ts_cal_p3 = self.calib_tod_rx(tod.data_pairs_diff_p3, rx)
        a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)
        t=tod.std
        x_az=tod.pointing.hor['az']


        if npairs<=len(ts_cal.pdiff[1,:]):
            nmax=npairs
        else:
            nmax=len(ts_cal.pdiff[1,:])
            print('Plotting the full set of good detectors.')


        for i in range (nmax):
            fs=14
            try:
                pl.figure(figsize=(14,10))
                #to plot Pair Diff from all pairs
                ax1 = pl.subplot(211)
                if x_axis == 't':
                    ax1.plot(t, ts_cal.pdiff[:,i], label='Before p3 filt')
                    ax1.plot(t, ts_cal_p3.pdiff[:,i], label='After p3 filt')
                elif x_axis == 'az':
                    ax1.plot(x_az, ts_cal.pdiff[:,i], label='Before p3 filt')
                    ax1.plot(x_az, ts_cal_p3.pdiff[:,i], label='After p3 filt')
                ax1.set_title('Pair Diff\nGCP_idx='+str(b_shift[i]+1), fontsize=fs)
                ax1.set_ylabel('T[K]', fontsize=fs)
                ax1.tick_params(axis='both', labelsize=fs)
                #ax1.set_xlim([datetime.datetime(2020, 4, 18, 18, 30, 00), datetime.datetime(2020, 4, 18, 19, 30, 00)])
                pl.legend(loc='lower right')
                #to plot Pair Sum from all pairs
                ax2 = pl.subplot(212)
                # ax2.plot(t, ts_cal.psum[:,i], label='Before p3 filt')
                # ax2.plot(t, ts_cal_p3.psum[:,i], label='After p3 filt')
                ax2.plot(x_az, ts_cal.psum[:,i], label='Before p3 filt')
                ax2.plot(x_az, ts_cal_p3.psum[:,i], label='After p3 filt')
                ax2.set_title('Pair Sum\nGCP_idx='+str(a_shift[i]+1), fontsize=fs)
                ax2.set_ylabel('T[K]', fontsize=fs)
                ax2.set_ylim(np.mean(ts_cal.psum[:,i])-5*np.std(), np.mean()+5*np.std())
                #ax2.set_xlim([datetime.datetime(2020, 4, 18, 18, 30, 00), datetime.datetime(2020, 4, 18, 19, 30, 00)])
                for ax in (ax1, ax2):
                    ax.tick_params(axis='both', labelsize=fs)

                pl.legend(loc='lower right')
                pl.suptitle(tag, fontsize=16)
                pl.savefig('tod/ba/'+tag+'/rx_'+str(rx)+'_gcp_idx_'+str(b_shift[i]+1)+'-'+str(a_shift[i]+1)+'.png')
                if not posting_folder=='None':
                    pl.savefig(posting_folder+'/'+tag+'_rx_'+str(rx)+'_ab-idx_'+str(i)+'.png')
                #pl.show()
                pl.close()
            except:
                print('i='+str(i)+' failed.')




    def pl_tod_atmo(self, tag, tod, rx, i_det, az_offs_det, p3=False, gs=False, i_det_savefig=0, posting_folder='None', wvr_fn='None', showplots=0):

        if posting_folder == 'None':
            posting_folder = self.pf

        # ts_cal, ts_cal_p3 = self.calib_tod_rx(tod, rx)
        # a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)
        t=np.array(tod.std)
        #t=t[tod.mapind] #not using mapind as fs.sf;fs.se should take care of this

        #x_az=x_az[tod.mapind]
        #T_rx=ts_cal.psum[tod.mapind,det]

        if p3==False:
            if gs==False:
                ts = tod.data_pairs_diff_noavg.fb
            else:
                ts = tod.data_pairs_diff_noavg_gs.fb
        elif p3==True:
            if gs==False:
                ts = tod.data_pairs_diff_noavg_p3.fb
            else:
                ts = tod.data_pairs_diff_noavg_p3_gs.fb

        x_az=tod.pointing.hor['az']

        x_az = np.array([x_az_i + az_offs_det for x_az_i in x_az])

        i_det_list = [i_det]

        for i_det in i_det_list:
            (det_a, det_b) = i_det

            ts_cal = self.calib_tod_rx_onets(tod, ts, rx, det_a, det_b)

            T_rx_pdiff=ts_cal.pdiff
            T_rx_psum=ts_cal.psum

            fs=tod.fs
            x_sum, D_sum=raz.interpToImage_BK(x_az, T_rx_psum, fs)
            x_diff, D_diff=raz.interpToImage_BK(x_az, T_rx_pdiff, fs)

            # x_labels=[t[fs_s_i] for fs_s_i in fs.s]

            #print('Nscan=', len(fs.sf))

            #cuttting to match PWV x x_axis
            #D_sum=D_sum[:,57:]
            #D_diff=D_diff[:,57:]

            if not wvr_fn == 'None':

                waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_fn)#, az_lim=(np.min(x_az), np.max(x_az)))

                fig,(ax1, ax2, ax3) = pl.subplots(3,1)

                pos1=ax1.imshow(D_sum, aspect='auto', interpolation='nearest', extent=[0, len(D_sum[1,:]), np.min(x_sum), np.max(x_sum)], origin='lower')
                ax1.set_title('Pair Sum')
                cbar1 = pl.colorbar(pos1, ax=ax1)
                if p3 == True:
                    cbar1.set_label('T[K]\np3 filtered')#, fontsize=fs)
                else:
                    cbar1.set_label('T[K]')
                ax1.set_ylabel('Az')
                #ax1.set_xticks([-0.75,-0.25,0.25,0.75])
                #ax1.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                pos2=ax2.imshow(D_diff, aspect='auto', interpolation='nearest', extent=[0, len(D_diff[1,:]), np.min(x_diff), np.max(x_diff)], origin='lower')
                ax2.set_title('Pair Diff')
                fig.suptitle(tag+'- rx = '+str(rx)+' - gcp_idx[a-b]='+str(a_shift[i_det]+1)+'-'+str(b_shift[i_det]+1))
                cbar2 = pl.colorbar(pos2, ax=ax2)
                if p3 == True:
                    cbar2.set_label('T[K]\np3 filtered')#, fontsize=fs)
                else:
                    cbar2.set_label('T[K]')
                ax2.set_xlabel('Nscan')
                ax2.set_ylabel('Az')
                #ax2.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                pos3=ax3.imshow(D_pwv, aspect='auto', interpolation='nearest', extent=[0, len(D_pwv[1,:]), np.min(x_diff), np.max(x_diff)], origin='lower')
                ax2.set_title('PWV')
                #fig.suptitle(tag+'- rx = '+str(rx)+' - gcp_idx[a-b]='+str(a_shift[i_det]+1)+'-'+str(b_shift[i_det]+1))
                cbar3 = pl.colorbar(pos3, ax=ax3)
                cbar3.set_label('PWV[um]')

                ax3.set_xlabel('Nscan')
                ax3.set_ylabel('Az')
                ax3.set_ylim(np.min(x_az), np.max(x_az))
                #pl.yticks(fontsize=fs_ticks)



                #pl.savefig('tod/ba/'+tag+'/rx_'+str(rx)+'_gcp_idx_'+str(det_b+1)+'-'+str(det_a+1)+'_zoom_off.png')
                if not posting_folder=='None':
                    pl.savefig(posting_folder+'/'+tag+'_rx_'+str(rx)+'_ab-i_'+str(i_det_savefig)+'_p3-'+str(p3)+'_zoom_off.png')

                #pl.show()
                pl.close()

            else:

                fig,(ax1, ax2) = pl.subplots(2,1)

                nscans_bk = len(D_sum[10,:])

                pos1=ax1.imshow(D_sum, aspect='auto', interpolation='nearest', extent=[0, len(D_sum[1,:]), np.min(x_sum), np.max(x_sum)], origin='lower')
                ax1.set_title('Pair Sum')
                cbar1 = pl.colorbar(pos1, ax=ax1)
                if p3 == True:
                    cbar1.set_label('T[K]\np3 filtered')#, fontsize=fs)
                else:
                    cbar1.set_label('T[K]')
                ax1.set_ylabel('Az')
                ax1.set_xlim(int(nscans_bk/2.), nscans_bk)
                #ax1.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                pos2=ax2.imshow(D_diff, aspect='auto', interpolation='nearest', extent=[0, len(D_diff[1,:]), np.min(x_diff), np.max(x_diff)], origin='lower')
                ax2.set_title('Pair Diff')
                fig.suptitle(tag+'- rx = '+str(rx)+' - gcp_idx[a-b]='+str(det_a+1)+'-'+str(det_b+1))
                cbar2 = pl.colorbar(pos2, ax=ax2)
                if p3 == True:
                    cbar2.set_label('T[K]\np3 filtered')#, fontsize=fs)
                else:
                    cbar2.set_label('T[K]')
                ax2.set_xlabel('Nscan')
                ax2.set_ylabel('Az')
                ax2.set_xlim(int(nscans_bk/2.),nscans_bk)
                #ax2.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                #pl.savefig('tod/ba/'+tag+'/rx_'+str(rx)+'_gcp_idx_'+str(b_shift[i_det]+1)+'-'+str(a_shift[i_det]+1)+'_zoom_off.png')
                if not posting_folder=='None':
                    pl.savefig(posting_folder+tag+'_rx_'+str(rx)+'_ab-i_'+str(i_det_savefig)+'_p3-'+str(p3)+'_zoom_off.png')


                #pl.show()
                pl.close()

        if not i_det == 'all':
            return x_az, T_rx_psum, T_rx_pdiff, D_sum, D_diff

    def pl_tod_atmo_pairs(self, tag, tod, rx, i_det, az_offs_det, i_det_savefig=0, posting_folder='None', wvr_fn='None', showplots=0):

        if posting_folder == 'None':
            posting_folder = self.pf

        t=np.array(tod.std)

        ts = tod.data_pairs.fb

        x_az=tod.pointing.hor['az']

        x_az = np.array([x_az_i + az_offs_det for x_az_i in x_az])

        i_det_list = [i_det]

        for i_det in i_det_list:
            (det_a, det_b) = i_det

            ts_cal = self.calib_tod_rx_onets(tod, ts, rx, det_a, det_b)

            T_rx_pdiff=ts_cal.pdiff #bdet
            T_rx_psum=ts_cal.psum   #adet

            fs=tod.fs
            x_sum, D_sum=raz.interpToImage_BK(x_az, T_rx_psum, fs)
            x_diff, D_diff=raz.interpToImage_BK(x_az, T_rx_pdiff, fs)

            # x_labels=[t[fs_s_i] for fs_s_i in fs.s]

            #print('Nscan=', len(fs.sf))

            #cuttting to match PWV x x_axis
            #D_sum=D_sum[:,57:]
            #D_diff=D_diff[:,57:]

            if not wvr_fn == 'None':

                waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_fn)#, az_lim=(np.min(x_az), np.max(x_az)))

                fig,(ax1, ax2, ax3) = pl.subplots(3,1)

                pos1=ax1.imshow(D_sum, aspect='auto', interpolation='nearest', extent=[0, len(D_sum[1,:]), np.min(x_sum), np.max(x_sum)], origin='lower')
                ax1.set_title('detpol A')
                cbar1 = pl.colorbar(pos1, ax=ax1)
                if p3 == True:
                    cbar1.set_label('T[K]\np3 filtered')#, fontsize=fs)
                else:
                    cbar1.set_label('T[K]')
                ax1.set_ylabel('Az')
                #ax1.set_xticks([-0.75,-0.25,0.25,0.75])
                #ax1.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                pos2=ax2.imshow(D_diff, aspect='auto', interpolation='nearest', extent=[0, len(D_diff[1,:]), np.min(x_diff), np.max(x_diff)], origin='lower')
                ax2.set_title('detpol B')
                fig.suptitle(tag+'- rx = '+str(rx)+' - gcp_idx[a-b]='+str(a_shift[i_det]+1)+'-'+str(b_shift[i_det]+1))
                cbar2 = pl.colorbar(pos2, ax=ax2)
                if p3 == True:
                    cbar2.set_label('T[K]\np3 filtered')#, fontsize=fs)
                else:
                    cbar2.set_label('T[K]')
                ax2.set_xlabel('Nscan')
                ax2.set_ylabel('Az')
                #ax2.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                pos3=ax3.imshow(D_pwv, aspect='auto', interpolation='nearest', extent=[0, len(D_pwv[1,:]), np.min(x_diff), np.max(x_diff)], origin='lower')
                ax2.set_title('PWV')
                #fig.suptitle(tag+'- rx = '+str(rx)+' - gcp_idx[a-b]='+str(a_shift[i_det]+1)+'-'+str(b_shift[i_det]+1))
                cbar3 = pl.colorbar(pos3, ax=ax3)
                cbar3.set_label('PWV[um]')

                ax3.set_xlabel('Nscan')
                ax3.set_ylabel('Az')
                ax3.set_ylim(np.min(x_az), np.max(x_az))
                #pl.yticks(fontsize=fs_ticks)



                #pl.savefig('tod/ba/'+tag+'/rx_'+str(rx)+'_gcp_idx_'+str(det_b+1)+'-'+str(det_a+1)+'_zoom_off.png')
                if not posting_folder=='None':
                    pl.savefig(posting_folder+'/singledet_atmo_'+tag+'_rx_'+str(rx)+'_ab-i_'+str(i_det_savefig)+'.png')

                #pl.show()
                pl.close()

            else:

                fig,(ax1, ax2) = pl.subplots(2,1)

                nscans_bk = len(D_sum[10,:])

                pos1=ax1.imshow(D_sum, aspect='auto', interpolation='nearest', extent=[0, len(D_sum[1,:]), np.min(x_sum), np.max(x_sum)], origin='lower')
                ax1.set_title('detpol A')
                cbar1 = pl.colorbar(pos1, ax=ax1)
                cbar1.set_label('T[K]')
                ax1.set_ylabel('Az')
                ax1.set_xlim(int(nscans_bk/2.), nscans_bk)
                #ax1.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                pos2=ax2.imshow(D_diff, aspect='auto', interpolation='nearest', extent=[0, len(D_diff[1,:]), np.min(x_diff), np.max(x_diff)], origin='lower')
                ax2.set_title('detpol B')
                fig.suptitle(tag+'- rx = '+str(rx)+' - gcp_idx[a-b]='+str(det_a+1)+'-'+str(det_b+1))
                cbar2 = pl.colorbar(pos2, ax=ax2)
                cbar2.set_label('T[K]')
                ax2.set_xlabel('Nscan')
                ax2.set_ylabel('Az')
                ax2.set_xlim(int(nscans_bk/2.),nscans_bk)
                #ax2.set_xticklabels(x_labels)
                #pl.yticks(fontsize=fs_ticks)

                #pl.savefig('tod/ba/'+tag+'/rx_'+str(rx)+'_gcp_idx_'+str(b_shift[i_det]+1)+'-'+str(a_shift[i_det]+1)+'_zoom_off.png')
                if not posting_folder=='None':
                    pl.savefig(posting_folder+'/singledet_atmo_'+tag+'_rx_'+str(rx)+'_ab-i_'+str(i_det_savefig)+'.png')


                #pl.show()
                pl.close()

        if not i_det == 'all':
            return x_az, T_rx_psum, T_rx_pdiff, D_sum, D_diff




    def plot_ts(self, rx_list, det, p3=True):

        fig, ax = pl.subplots(2, 1, figsize=(12,8))

        tod=self.bk_struct
        x_az=tod.pointing.hor['az']

        colors = ['c', 'm', 'y', 'k']
        i=0

        bk_time_ordered_az=[]
        bk_time_ordered_time=[]
        bk_time_ordered_data_psum=[]
        bk_time_ordered_data_pdiff=[]

        for rx in rx_list:
            try:

                ts_cal = self.calib_tod_rx(tod.data_pairs_diff_noavg, rx)
                ts_cal_p3 = self.calib_tod_rx(tod.data_pairs_diff_noavg_p3, rx)

                a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)
                t=np.array(tod.std)

                time_ordered_az=x_az[tod.mapind]

                if p3 == True:
                    time_ordered_data_pdiff=ts_cal_p3.pdiff[tod.mapind, det]
                    time_ordered_data_psum=ts_cal_p3.psum[tod.mapind, det]
                else:
                    time_ordered_data_pdiff=ts_cal.pdiff[tod.mapind, det]
                    time_ordered_data_psum=ts_cal.pdiff[tod.mapind, det]


                time_ordered_time=t[tod.mapind]
                dt = (time_ordered_time[1]-time_ordered_time[0]).total_seconds()
                daz = time_ordered_az[1]-time_ordered_az[0]

                print(time_ordered_az)
                while (np.min(time_ordered_az)<0.):
                    time_ordered_az = [toaz_i+360. for toaz_i in time_ordered_az]
                    print('adding 360')

                az_bk_max = np.max(time_ordered_az)
                az_bk_min = np.min(time_ordered_az)

                print('bk time ordered az = ', time_ordered_az)

                bk_time_ordered_az.append(time_ordered_az)
                bk_time_ordered_data_pdiff.append(time_ordered_data_pdiff)
                bk_time_ordered_data_psum.append(time_ordered_data_psum)
                bk_time_ordered_time.append(time_ordered_time)

                #rx_res = res_150 * (150./rx)

                ax[0].plot(time_ordered_az, time_ordered_data_pdiff, label='rx = '+str(rx)+' - psum')
                ax[0].plot(time_ordered_az, time_ordered_data_psum, label='rx = '+str(rx)+' - pdiff')
                i+=1

            except:
                print('rx '+str(rx)+' failed.')

        ax[0].set_title('BK\n'+self.bk_tag)
        ax[0].set_ylabel('T[?]')
        ax[0].legend()



        #WVR Timestream

        time_ordered_time=self.wvr_struct.tot
        dt_wvr=30.
        #time_ordered_az=time_ordered_time*(360./30.)
        time_ordered_az=self.wvr_struct.az_wvr

        ps_matrix_az=[]

        date=self.wvr_scan[:8]
        time=self.wvr_scan[9:11]

        rms_i=[]
        daz=time_ordered_az[10]-time_ordered_az[9]

        wvr_time_ordered_az=[]
        wvr_time_ordered_data=[]
        wvr_time_ordered_time=[]

        pwv_full_ts = self.wvr_struct.pwv_full_ts


        for scan_i in range(110):

            pwv_ts_onescan = pwv_full_ts[self.wvr_struct.fs.s[scan_i]:self.wvr_struct.fs.e[scan_i]]
            time_ordered_az_onescan = time_ordered_az[self.wvr_struct.fs.s[scan_i]:self.wvr_struct.fs.e[scan_i]]
            time_ordered_time_onescan=time_ordered_time[self.wvr_struct.fs.s[scan_i]:self.wvr_struct.fs.e[scan_i]]-time_ordered_time[self.wvr_struct.fs.s[scan_i]]


            bk_match_mask = np.where((time_ordered_az_onescan>=az_bk_min) & (time_ordered_az_onescan<=az_bk_max))

            pwv_ts_onescan = pwv_ts_onescan[bk_match_mask]
            time_ordered_az_onescan = time_ordered_az_onescan[bk_match_mask]
            time_ordered_time_onescan = time_ordered_time_onescan[bk_match_mask]

            rms_i.append(np.nanstd(pwv_ts_onescan))

            mask = np.isfinite(pwv_ts_onescan)
            pwv_ts_onescan_masked = pwv_ts_onescan[mask]
            time_ordered_az_onescan_masked = time_ordered_az_onescan[mask]
            time_ordered_time_onescan_masked = time_ordered_time_onescan[mask]

            time_interp = np.linspace(0., 30., 3000)
            tck_time = interpolate.splrep(time_ordered_time_onescan_masked, pwv_ts_onescan_masked)
            ts_interp = interpolate.splev(time_interp, tck_time)

            tck_az = interpolate.splrep(time_ordered_time_onescan_masked, time_ordered_az_onescan_masked)
            az_interp = interpolate.splev(time_interp, tck_az)

            wvr_time_ordered_az.append(az_interp[:-100])
            wvr_time_ordered_data.append(ts_interp[:-100])
            wvr_time_ordered_time.append(time_interp[:-100])

            #ax[1].plot(az_interp[:-100], ts_interp[:-100])
            ax[1].plot(time_ordered_az_onescan_masked, pwv_ts_onescan_masked)

        ax[1].set_title('WVR\n'+self.wvr_scan)
        ax[1].set_ylabel('PWV[um]')
        pl.suptitle('Timestreams')
        pl.close()

        return wvr_time_ordered_time, wvr_time_ordered_az, wvr_time_ordered_data, bk_time_ordered_time, bk_time_ordered_az, bk_time_ordered_data_psum, bk_time_ordered_data_pdiff



    def plot_PS(self, rx_list, det):


        #
        # for rx in rx_list:
        #     ts_cal, ts_cal_p3 = self.calib_tod_rx(tod, rx)
        #     i_det_list = [j for j in range (0, len(ts_cal.pdiff[1,:]))]
        #     for j in i_det_list:
        #         D_sum, D_diff = self.pl_tod_atmo(tag, tod, rx, i_det=j, p3=False, posting_folder='None')
        #         D_sum_p3, D_diff_p3 = self.pl_tod_atmo(tag, tod, rx, i_det=j, p3=True, posting_folder='None')
        #         #raz.extract_PowerSpectrum(D_pwv, D_sum, D_diff, az_pwv, x_az)
        #


        wvr_time_ordered_time, wvr_time_ordered_az, wvr_time_ordered_data, bk_time_ordered_time, bk_time_ordered_az, bk_time_ordered_data_psum, bk_time_ordered_data_pdiff = self.plot_ts(rx_list, det)

        colors = ['c', 'm', 'y', 'k']
        i=0

        fig, ax = pl.subplots(2, 1, sharex=True, figsize=(12,8))
        for rx in rx_list:
            try:
                # ts_cal, ts_cal_p3 = self.calib_tod_rx(tod, rx)
                # a_shift, b_shift = self.find_rgl_idx(rx, tod.ind)
                # t=np.array(tod.std)

                # time_ordered_az=x_az[tod.mapind]
                # time_ordered_data=ts_cal_p3.pdiff[tod.mapind, det]
                # time_ordered_time=t[tod.mapind]

                time_ordered_time = bk_time_ordered_time[i]
                time_ordered_data_psum = bk_time_ordered_data_psum[i]
                time_ordered_data_pdiff = bk_time_ordered_data_pdiff[i]
                time_ordered_az = bk_time_ordered_az[i]

                dt = (time_ordered_time[1]-time_ordered_time[0]).total_seconds()
                daz = time_ordered_az[1]-time_ordered_az[0]
                #
                # pl.figure()
                # pl.plot(time_ordered_time, time_ordered_data)
                # pl.xlabel('t[s]')
                # pl.title('dt = '+str(dt))
                # pl.suptitle('TOD (vs time)')

                # pl.plot(time_ordered_az, time_ordered_data, alpha=0.5)
                # pl.xlabel('Az[deg]')
                # pl.title('daz = '+str(daz))
                # pl.suptitle('TOD (vs az)')


                ps_psum=raz.extract_PowerSpectrum(time_ordered_time, dt, time_ordered_az, daz, time_ordered_data_psum, data='BK')
                ps_pdiff=raz.extract_PowerSpectrum(time_ordered_time, dt, time_ordered_az, daz, time_ordered_data_pdiff, data='BK')
                # #the dt for BK is not homogeneus in the TOD --> the time domain PS is not accurate
                #
                #
                # pl.plot(ps.t.freq, ps.t.ps, c='k', alpha=0.2)
                # pl.scatter(ps.t.freq, ps.t.ps)
                # pl.loglog()
                # pl.xlabel('t[s]')
                # pl.suptitle('Time Power Spectrum')

                # rx_res = res_150 * (150./rx)

                ax[1].plot(ps_psum.az.freq, ps_psum.az.ps, alpha=0.2)
                #ax[1].axvline(x = 1./rx_res, c=colors[i], label='rx'+str(rx)+' ')
                #ax[1].axvline(x = 1./az_span, c='blue')
                ax[1].scatter(ps_psum.az.freq, ps_psum.az.ps, label='psum')#, c=colors[i])

                ax[1].plot(ps_pdiff.az.freq, ps_pdiff.az.ps, alpha=0.2)
                #ax[1].axvline(x = 1./rx_res, c=colors[i], label='rx'+str(rx)+' ')
                #ax[1].axvline(x = 1./az_span, c='blue')
                ax[1].scatter(ps_pdiff.az.freq, ps_pdiff.az.ps, label='pdiff')#, c=colors[i])


                # ax[0].set_title('BK ts')
                # pl.show()
                i+=1

            except:
                print('rx '+str(rx)+' failed.')

        #pl.loglog()
        #pl.show()

        ax[1].set_yscale('log')
        ax[1].set_xscale('log')

        ax[1].set_xlabel('Az[deg]')
        ax[1].set_title('BK Power Spectrum')
        ax[1].legend()

        # data.bk_struct.tod=bk_time_ordered_data
        # data.bk_struct.toaz=bk_time_ordered_az
        # data.bk_struct.rx=rx_list

        #data['bk_struct']=bk_struct

        # ax[1].legend(title='Rx Resolution')
        # ax[1].set_yscale('log')


        #WVR Timestream

        # fn_posting=pf+'/'+wvr_scan+'_pwv_zoom_on.png'

        #waz, az_wvr, az_real, fs, idx, pwv_full_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan, show=0)

        #Extract Time Domain PS
        #time_ordered_time=np.array([i*30. for i in D_pwv[1,:]])
        # time_ordered_time=tot
        # dt_wvr=30.
        # #time_ordered_az=time_ordered_time*(360./30.)
        # time_ordered_az=az_wvr
        #
        ps_matrix_az=[]
        #
        # date=wvr_scan[:8]
        # time=wvr_scan[9:11]
        #
        # # time_ordered_az = az_real
        #
        # rms_i=[]
        # daz=time_ordered_az[10]-time_ordered_az[9]
        #
        # wvr_time_ordered_az=[]
        # wvr_time_ordered_data=[]


        for scan_i in range(110):

            #pwv_ts_onescan= D_pwv[:,scan_i]
            # pwv_ts_onescan = pwv_full_ts[fs.s[scan_i]:fs.e[scan_i]]
            # time_ordered_az_onescan = time_ordered_az[fs.s[scan_i]:fs.e[scan_i]]
            # time_ordered_time_onescan=time_ordered_time[fs.s[scan_i]:fs.e[scan_i]]-time_ordered_time[fs.s[scan_i]]

            pwv = wvr_time_ordered_data[scan_i]
            az = wvr_time_ordered_az[scan_i]
            time = wvr_time_ordered_time[scan_i]


            # mask_az_bk = np.where((time_ordered_az_onescan > az_bk_min) & (time_ordered_az_onescan < az_bk_max))

            # rms_i.append(np.nanstd(pwv_ts_onescan))
            # #pl.plot(pwv_ts_onescan)
            # # pl.scatter(time_ordered_az_onescan, pwv_ts_onescan, label='rms='+str(np.nanstd(pwv_ts_onescan)))
            # # pl.legend()
            # # pl.show()
            #
            # mask = np.isfinite(pwv_ts_onescan)
            # # print('tot_time = ', time_ordered_time[fs.e[scan_i]]-time_ordered_time[fs.s[scan_i]])
            # pwv_ts_onescan_masked = pwv_ts_onescan[mask]
            # time_ordered_az_onescan_masked = time_ordered_az_onescan[mask]
            # time_ordered_time_onescan_masked = time_ordered_time_onescan[mask]
            #
            # time_interp = np.linspace(0., 30., 3000)
            # tck_time = interpolate.splrep(time_ordered_time_onescan_masked, pwv_ts_onescan_masked)
            # ts_interp = interpolate.splev(time_interp, tck_time)
            # #
            # # pl.plot(time_interp[:-100], ts_interp[:-100], c='r', alpha=0.5)
            # # pl.scatter(time_ordered_time_onescan_masked[:-10], pwv_ts_onescan_masked[:-10], c='k')
            #
            # tck_az = interpolate.splrep(time_ordered_time_onescan_masked, time_ordered_az_onescan_masked)
            # az_interp = interpolate.splev(time_interp, tck_az)

            #(scannum, wrapped_az) = divmod(time_ordered_az_masked_oneturn,360)
            #time_ordered_time below is wrong but I don't care because I just want the space domain

            #ps_pwv=raz.extract_PowerSpectrum(time_ordered_time_onescan_masked, dt_wvr, time_ordered_az_onescan_masked, daz, pwv_ts_onescan_masked)

            #daz=time_ordered_az[10]-time_ordered_az[9]

            ps_pwv=raz.extract_PowerSpectrum(time, time[10]-time[9], az, az[10]-az[9], pwv)

            ps_matrix_az.append(ps_pwv.az.ps)

            #f_ps_interp = interp1d(ps_pwv.az.freq, ps_pwv.az.ps, kind='cubic')
            # freq_axis_interp=np.linspace(1./360., 1./daz, 1000)
            # #ps_interp = f_ps_interp(freq_axis_interp)
            #
            # tck = interpolate.splrep(ps_pwv.az.freq, ps_pwv.az.ps) #should interrp timestream
            #
            # ps_interp = interpolate.splev(freq_axis_interp, tck, der=0)
            ax[0].plot(ps_pwv.az.freq[:-10], ps_pwv.az.ps[:-10], alpha=0.2)
            #ax[1].axvline(x = 1./rx_res, c=colors[i], label='rx'+str(rx)+' ')
            ax[0].scatter(ps_pwv.az.freq[:-10], ps_pwv.az.ps[:-10])


        #pl.show()

            # pl.scatter(ps_pwv.az.freq, ps_pwv.az.ps)#, label='az='+str(a))
            # pl.plot(freq_axis_interp, ps_interp, c='r')



        freq_avg=ps_pwv.az.freq
        # ps_matrix_az_array=np.full((len(ps_matrix_az[1]), len(ps_matrix_az)), np.nan)
        #
        # ps_avg=np.zeros(len(freq_avg))
        # #
        # for i in range (len(ps_matrix_az)-1):
        #         #try:
        #         ps_matrix_az_array[:,i]=ps_matrix_az[i]
        #         #except:
        #         #    print(str(i)+' failed.')
        #
        # for k in range (len(freq_avg)-1):
        #     ps_avg[k]=np.nanmean(ps_matrix_az_array[k,:])

        #
        # data.wvr_struct.tod=wvr_time_ordered_data
        # data.wvr_struct.toaz=wvr_time_ordered_az

        #data['wvr_struct']=wvr_struct

        #print(data)

        # ax[0].plot(freq_avg[:-1], ps_avg[:-1], c='r')#, label='az='+str(a))
        #ax[1].axvline(x = 1/2.5, c='y', label='WVR beam')
        #ax[1].legend()
        ax[0].set_yscale('log')
        ax[0].set_xscale('log')
        #ax[1].set_title('WVR power spectrum')

        pl.loglog()

        pl.show()



    def plot_bk_atmo(self, rx_list, wvr_fn='None'):
        tod = self.bk_struct
        for rx in rx_list:
            self.pl_tod_all(self.bk_tag, rx,  npairs=10, posting_folder=self.pf)
            self.pl_tod_pair(self.bk_tag, rx, npairs=10, posting_folder=self.pf)
            self.pl_tod_atmo(self.bk_tag, tod, rx, posting_folder=self.pf,  wvr_fn=wvr_fn)



    def correlate_az_tod(self, rx):

        wvr_atmo = self.wvr_struct.D_pwv
        az_wvr = self.wvr_struct.az_real #az_calib
        time_ordered_az=self.wvr_struct.az_wvr
        fs = self.wvr_struct.fs


        az_star=250
        f = open('pointing_parameters_2018_fast.txt','rb')
        point_par = pk.load(f)
        f.close()

        az_offs = point_par['az_offs']
        az_star_wvr = az_star + az_offs
        if az_star_wvr > 360:
            az_star_wvr = az_star_wvr - 360.

        nscans=len(wvr_atmo[10,:])
        tod_wvr = np.zeros(nscans)

        for iscan in range(nscans-1):

            az_wvr_onescan = az_wvr[fs.s[iscan]:fs.e[iscan]]
            time_ordered_az_onescan = time_ordered_az[fs.s[iscan]:fs.e[iscan]]
            x_az_atmo=np.arange(len(wvr_atmo[:,10]))
            mask_wvr=np.where(x_az_atmo>az_star_wvr)[0]
            idx_wvr=mask_wvr[0]

            print('mask_wvr=', mask_wvr)
            print('idx_wvr=', idx_wvr)

            scan_tod = wvr_atmo[:,iscan]
            scan_tod_az = scan_tod[idx_wvr]

            print('scan_tod_az = ', scan_tod_az)

            tod_wvr[iscan] = scan_tod_az


        # pl.plot(tod_wvr)
        # pl.xlabel('scan_n')
        # pl.show()

        tod = self.bk_struct
        fs_bk = tod.fs

        az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = self.pl_tod_atmo(self.bk_tag, tod, rx, i_det, az_offs_det)

        print(az_bk)
        while (np.min(az_bk)<0.):
            az_bk = [toaz_i+360. for toaz_i in az_bk]
            print('adding 360')
        print(az_bk)

        az_bk = np.array(az_bk)
        bk_atmo = D_sum_bk

        print('az_star = ', az_star)

        nscans_bk = len(bk_atmo[10,:])
        tod_bk = np.zeros(nscans_bk)

        for iscan in range(nscans_bk-1):

            az_bk_onescan = az_bk[fs_bk.sf[iscan]:fs_bk.ef[iscan]]
            mask_bk=np.where(az_bk>az_star)[0]
            idx_bk=mask_bk[0]

            scan_tod_bk = bk_atmo[:,iscan]
            pl.plot(scan_tod_bk)
            #pl.title('scan '+str(iscan))
            scan_tod_az_bk = scan_tod_bk[idx_bk]


            print('scan_tod_az_bk = ', scan_tod_az_bk)

            tod_bk[iscan] = scan_tod_az_bk

        pl.xlabel('Az')
        #pl.show()
        pl.close()

        # pl.plot(tod_bk)
        # pl.xlabel('scan_n')
        # pl.show()
        #
        # pl.plot(tod_bk, label='bk')
        # pl.plot(tod_wvr, label='wvr')
        # pl.legend()
        # pl.xlabel('scan_n')
        # pl.title('TOD compared')
        # pl.show()

    def BK_WVR_ts_corr(self, bk_atmo_diff, bk_atmo_sum, wvr_scan, save_fn, save_txt):


        fs_wvr = self.wvr_struct.fs
        tod = self.bk_struct
        fs_bk = tod.fs
        waz, az, calib_az, fs, idx, pwv_ts, wvr_atmo, tot = raz.read_pwvatmo(wvr_scan)
        time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)
        time_UTC=[]
        time_UTC=[parser.parse(time_str) for time_str in time_UTC_str]
        full_time_wvr=[time_UTC[int(fs_s_i)] for fs_s_i in fs_wvr.s]
        xlabels_wvr=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time_wvr])

        t_min_wvr=datetime.time(full_time_wvr[0])
        t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
        t_min_bk=datetime.time(full_time[0])
        t_max_bk=datetime.time(full_time[len(full_time)-1])
        x_diff, D_diff = raz.interpToImage_BK(az_bk, T_rx_pdiff_bk, fs_bk)
        x_min_wvr_idx = np.where(calib_az>= np.min(x_diff))[0][0]
        x_max_wvr_idx = np.where(calib_az>= np.max(x_diff))[0][0]

        t_bk=np.array(tod.std)

        full_time=[t_bk[int(fs_s_i)] for fs_s_i in fs_bk.sf]
        xlabels=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time])

        full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])
        full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr])

        if (t_min_bk <= t_min_wvr):

            bk_time_mask = np.where(full_time_dt>=t_min_wvr)[0]
            wvr_time_mask = np.where(full_time_wvr_dt<=t_max_bk)[0]

        else:

            bk_time_mask = np.where(full_time_dt<=t_max_wvr)[0]
            wvr_time_mask = np.where(full_time_wvr_dt>=t_min_bk)[0]


        xlabels_mask=xlabels[bk_time_mask]
        xlabels_wvr_mask=xlabels_wvr[wvr_time_mask]

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
        pl.savefig(save_fn+'_Pearson_Corr_Coeff_v2.png')
        pl.close()

        Pearson_Corr_Coeff = {'psum-wvr':corr_coeff_psum_wvr_list, 'pdiff-wvr':corr_coeff_pdiff_wvr_list, 'psum-pdiff':corr_coeff_psum_pdiff_list, 'time_axis':xaxis_corr}

        f = open(save_txt,'wb')
        pk.dump(Pearson_Corr_Coeff, f)
        f.close()

    def corr_map_fpu(self, wvr_atmo_Trj, pf, rx_list=[210, 270]):

        def r_fpu(bk_tag, wvr_scan, rx, p3_filt, nscan, pf):

            x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= self.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')

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

        tod=self.bk_struct
        fs_bk=tod.fs
        #Extracting Trj atmograms from WVR PWV
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
        wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

        for rx in rx_list:
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

                    print(bk_tag+'_psum-wvr_FPU_'+str(rx)+'_time'+str(time_stamp)+'_p3_'+str(p3_filt)+'.png')

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
                    pl.close()



                pl.figure()
                n, x, _ = pl.hist(psum_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'stepfilled', facecolor='g', alpha=0.6)
                n, x, _ = pl.hist(psum_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'step', color='g')
                pl.suptitle('Psum-Pdiff Correlation Histogram')
                pl.xlabel('r')
                pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_psum-wvr.png')
                pl.show()

                pl.figure()
                n, x, _ = pl.hist(pdiff_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'stepfilled', facecolor='g', alpha=0.6)
                n, x, _ = pl.hist(pdiff_wvr_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'step', color='g')
                pl.suptitle('Psum-Pdiff Correlation Histogram')
                pl.xlabel('r')
                pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_pdiff-wvr.png')
                pl.show()

                pl.figure()
                n, x, _ = pl.hist(psum_pdiff_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'stepfilled', facecolor='g', alpha=0.6)
                n, x, _ = pl.hist(psum_pdiff_tlist.flatten(), bins=np.linspace(-1, 1, 50), histtype=u'step', color='g')
                pl.suptitle('Psum-Pdiff Correlation Histogram')
                pl.xlabel('r')
                pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_histo_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_psum-pdiff.png')
                pl.show()


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
                pl.close()

                #
                time_stamp = 't_std'

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
                pl.close()






        #
        # def gauss(x, *p):
        #     A, mu, sigma, offs = p
        #     return A*np.exp(-(x-mu)**2/(2.*sigma**2))+(offs)
        #
        # def compare_one_az:
        #
        #
        #
        #
        # def correlate_tod(tod_wvr, tod_bk, az_wvr, az_bk, az_star):
        #     #for k in range(1,npts-1):
        #
        #     rowcor1=np.correlate(atmogram[k,:],atmogram[k+delta_az,:], mode='full')/(atmogram[k,:].std()*atmogram[k+delta_az,:].std())
        #     #delay[k]=rowcor1.argmax()
        #     pl.plot(lag, rowcor1, label='Az'+str(k)+'-'+str(k+20))
        #     try:
        #         fit_mask=np.where((lag>=-6) & (lag<=6))
        #         p0=[100., 0., 5., 0.]
        #         coeff, var_matrix = curve_fit(gauss, lag[fit_mask], rowcor1[fit_mask], p0=p0, bounds=((0, -100, 0, 0), (1000, 100, 10, 100)))
        #         x_model=np.linspace(np.min(lag[fit_mask]), np.max(lag[fit_mask]), 100)
        #         g=gauss(x_model, *coeff)
        #         delay[k]=coeff[1]
        #             #delay[k]=nscans-rowcor1.argmax()
        #     except:
        #         delay[k]=nscans-rowcor1.argmax()
        #         #pl.plot(x_model, g, c=c_set[k], label='gaussian model:\nmu[deg] ='+str(round(coeff[1],1)))
        #         #pl.xlim(-200,200)
        #
        #
        # mask_bk=np.where(az_bk>az_star)[0]
        # idx_bk=mask_bk[0]
        #
        #
        #
        #
        #
        # lag=np.arange(2*nscans-1)-nscans
        # npts=len(atmogram[:,1])-delta_az
        # #npts=5 #just for test
        #
        # delay=np.zeros(npts)
        # phase=np.arange(npts)
        #
        # c_set=['r', 'g', 'y', 'c', 'b']
        #
        # pl.figure(figsize=(15,8))
        # for k in range(1,npts-1):
        #     rowcor1=np.correlate(atmogram[k,:],atmogram[k+delta_az,:], mode='full')/(atmogram[k,:].std()*atmogram[k+delta_az,:].std())
        #     #delay[k]=rowcor1.argmax()
        #     pl.plot(lag, rowcor1, label='Az'+str(k)+'-'+str(k+20))
        #     try:
        #         fit_mask=np.where((lag>=-6) & (lag<=6))
        #         p0=[100., 0., 5., 0.]
        #         coeff, var_matrix = curve_fit(gauss, lag[fit_mask], rowcor1[fit_mask], p0=p0, bounds=((0, -100, 0, 0), (1000, 100, 10, 100)))
        #         x_model=np.linspace(np.min(lag[fit_mask]), np.max(lag[fit_mask]), 100)
        #         g=gauss(x_model, *coeff)
        #         delay[k]=coeff[1]
        #         #delay[k]=nscans-rowcor1.argmax()
        #     except:
        #         delay[k]=nscans-rowcor1.argmax()
        #     #pl.plot(x_model, g, c=c_set[k], label='gaussian model:\nmu[deg] ='+str(round(coeff[1],1)))
        #     #pl.xlim(-200,200)
        #
        # #pl.legend()
        # pl.xlabel('lag', fontsize=18)
        # pl.ylabel('rowcor/nscans', fontsize=18)
        # pl.suptitle('Correlation Function', fontsize=18)
        # pl.title('delta_az = '+str(delta_az), fontsize=18)
        # pl.xticks(fontsize=18)
        # pl.yticks(fontsize=18)
        # pl.suptitle(wvr_scan[:-4], fontsize=18)
        # pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_correlation_plots_deltaAz_'+str(delta_az)+'.png')
        # if self.show==1:
        #     pl.show()
        # else:
        #     pl.close()
        #
        # #mask_good=np.where(((nscans-delay)<7) & ((nscans-delay)>-7))
        # #mask_good=np.where((delay<10) & (delay>-10))
        #
        # def sin_func(x, A, phi, C):
        #     return A*np.sin((x*2.*np.pi/np.max(phase))+phi)+C
        #
        # p0=[1., 0., 0.]
        # #coeff, var_matrix = curve_fit(sin_func, phase[mask_good], (nscans-delay)[mask_good], p0=p0, bounds=((0., -2.*np.pi, 0.), (np.inf, 2.*np.pi, np.inf)))
        # coeff, var_matrix = curve_fit(sin_func, phase, delay, p0=p0, bounds=((0., -2.*np.pi, -np.inf), (np.inf, 2.*np.pi, np.inf)))
        # model=sin_func(phase, *coeff)
        #
        #
        # noise = delay - model
        # sigma = np.std(noise)
        # pl.figure()
        # n, x, _ = pl.hist(noise, bins=np.linspace(-10, 10, 40), histtype=u'step', density=True)
        # pl.suptitle('Noise 1')
        # pl.title('sigma = '+ str(round(sigma,2)))
        # #pl.show()
        # pl.close()
        #
        # pl.figure(figsize=(15,8))
        # pl.suptitle('Noise Cut', fontsize=18)
        # pl.plot(phase, noise, label='noise')
        # sigma2_cut = np.where((noise <= 2*sigma) & (noise >= (-2*sigma)))
        # pl.axhline(y=1*sigma, ls='--', alpha=0.5, label='2sigma')
        # noise2= noise[sigma2_cut]
        # sigma2=np.std(noise2)
        # pl.close()
        # pl.figure()
        # n, x, _ = pl.hist(noise2, bins=np.linspace(-10, 10, 40), histtype=u'step', density=True)
        # pl.suptitle('Noise 2')
        # pl.title('sigma = '+ str(round(sigma2,2)))
        # pl.close()
        #
        # pl.plot(phase[sigma2_cut], noise[sigma2_cut], label='noise_cut1')
        # pl.axhline(y=1*sigma2, ls='--', alpha=0.5, label='2sigma2')
        # stage2_sigma2_cut = np.where((noise <= 2*sigma2) & (noise >= (-2*sigma2)))
        # pl.plot(phase[stage2_sigma2_cut], noise[stage2_sigma2_cut], label='noise_cut2')
        # pl.xticks(fontsize=18)
        # pl.yticks(fontsize=18)
        # # noise3=noise[stage2_sigma2_cut]
        # # sigma3=np.std(noise3)
        # # pl.axhline(y=2*sigma3, ls='--', alpha=0.5, label='2sigma2')
        # # new_sigma2_cut = np.where((noise <= 2*sigma2) & (noise >= (-2*sigma2)))
        # # pl.plot(phase[new_sigma2_cut], noise[new_sigma2_cut], label='noise_cut2')
        #
        # pl.legend(fontsize=18)
        # #pl.show()
        # pl.close()
        #
        # data_fit_mask = stage2_sigma2_cut
        # p0=[1., 0., 0.]
        # coeff, var_matrix = curve_fit(sin_func, phase[data_fit_mask], delay[data_fit_mask], p0=p0, bounds=((0., -2.*np.pi, -np.inf), (np.inf, 2.*np.pi, np.inf)))
        # print('offset=', coeff[2])
        # print('avg=', np.mean(delay[data_fit_mask]))
        # model=sin_func(phase, *coeff)
        # tau = (coeff[0]*30) #1scan=30s
        # h_clouds=1900#m
        # delta_x_sky = h_clouds*np.sin(np.radians(delta_az))*(1./np.tan(math.radians(self.el)))
        # ws = delta_x_sky/tau
        # print('wind speed =', ws)
        # wd = np.degrees(coeff[1])+(delta_az/2.) #midpoint between 2 correlated timestreams
        # if wd<0:
        #     wd=360.+wd
        # print('wind dir =', wd)
        #
        #
        # pl.figure(figsize=(15,8))
        # #pl.plot(phase[mask_good],(nscans-delay)[mask_good],'-bs', markersize=5, linewidth=0.7, c='blue')
        # pl.plot(phase[data_fit_mask], delay[data_fit_mask],'-bs', markersize=5, linewidth=0.7, c='blue')
        # pl.ylim(-10,10)
        # pl.xlabel('Azimuthal Phase [deg]', fontsize=18)
        # pl.ylabel('Delay [s/30]', fontsize=18)
        # pl.xticks(fontsize=18)
        # pl.yticks(fontsize=18)
        # pl.suptitle(wvr_scan[:-4], fontsize=18)
        # pl.title('ws[m/s] ='+str(round(ws,2))+'\nwd[deg] = '+str(round(wd,2)), fontsize=18)
        #
        # label_real='Amp[s/30]='+str(round(coeff[0],2))+'\n phase[deg]='+str(round(np.degrees(coeff[1]),2))
        # pl.plot(phase, model, c='r', label=label_real)
        # pl.legend(fontsize=18)
        # pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_delay_deltaAz_'+str(delta_az)+'.png')
        # if self.show==1:
        #     pl.show()
        # else:
        #     pl.close()
        #
        # return ws, wd
