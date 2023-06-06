
import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import numpy as np
import pickle as pk
import pylab as pl
import scipy.optimize as sp
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec
from itertools import groupby
from operator import itemgetter
from math import atan, atan2
import math

#from initialize import initialize
#from wvr_pipeline import wvrScience

#wvrS=wvrScience.wvrScience()

class ReadAzscan(object):
    '''

    '''

    def __init__(self, unit=None, verb=True):

        '''


        '''
        #initialize.__init__(self, unit, verb=verb)
        self.T_hot=363.15
        self.T_cold=283.15
        self.LO=91.65
        self.IF=2*self.LO

        self.ch=[1.05,1.9,3.2,6.1]
        self.freqs=[self.IF-self.ch[3],self.IF-self.ch[2],self.IF-self.ch[1],self.IF-self.ch[0],self.IF+self.ch[0],self.IF+self.ch[1],self.IF+self.ch[2],self.IF+self.ch[3]]
        self.path_to_all = ''
        self.path_to_wvr1data = ''

    def findScans(self, az, nwrap=360):
        '''
        Temporary copied from wvrScience
        '''
        class stru:
            def __init__(self):
                self.num = []
                self.s = []
                self.e = []

        naz = size(az)
        (scannum, az) = divmod(az,nwrap)
        print('scannum=', scannum)
        s = [next(group)[0] for key, group in groupby(enumerate(scannum), key=itemgetter(1))]  # indices of start of scan
        e = s[1:] ; e.append(naz) # indices of end of scan

        fs = stru()
        fs.num = scannum
        fs.s = s
        fs.e = e

        #TODO: remove first and last scan ?
        return az, fs


    def interpToImage(self, waz, tsrc, fs):

        nscans = len(fs.s)

        D = zeros([361,nscans])
        y = copy(tsrc)
        # for each 360-scan
        for i in range(nscans):
            s = fs.s[i];e=fs.e[i]
            idx=range(s,e)
            yp = np.interp(arange(0,361), waz[idx], y[idx], right=nan, left=nan)
            D[:,i]=yp

        return D


    def interpToImage_BK(self, x_az, T_rx, fs):
        #nscans = len(fs.s)
        az_range = np.abs(x_az[int(fs.ef[1])] -x_az[int(fs.sf[1])])
        D = zeros([int(az_range)+1, len(fs.sf)])
        y = copy(T_rx)
        for i in range (len(fs.sf)):
            s = int(fs.sf[i]);e=int(fs.ef[i])
            idx=range(s,e)
            x_interp=arange(np.min(x_az[idx]), np.max(x_az[idx]))
            if (i % 2) == 0:
                yp = np.interp(x_interp, x_az[idx], y[idx], right=nan, left=nan)
            else:
                yp = np.interp(x_interp, np.flip(x_az[idx]), np.flip(y[idx]), right=nan, left=nan)
            D[:,i]=yp
        return x_interp, D

    def return_waz(self, x_az, fs):
        waz=np.zeros(len(x_az))
        for i in range (len(fs.sf)):
            s = int(fs.sf[i]);e=int(fs.ef[i])
            idx=range(s,e)
            x_interp=arange(np.min(x_az[idx]), np.max(x_az[idx]))
            if (i % 2) == 0:
                waz[idx] = x_az[idx]
            else:
                waz[idx] = np.flip(x_az[idx])

        return np.array(x_interp), waz



    def read_Az_slow(self, filename, pathtofn='', show_plots=0):
        with open('wvr1_data/'+pathtofn+filename) as f:
             lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
             FH = np.loadtxt(lines, dtype='str', delimiter=' ', skiprows=1, usecols = (0,1,19,20,21,22,23,24))
             time=FH[:,0]
             print('time_slow=', time[:6])

             data=np.zeros((len(time), 7))
             for i in range(7):
                 FH_i=FH[:,i+1]
                 list_of_floats=[]
                 for item in FH_i:
                     list_of_floats.append(float(item))
                 data[:,i]=list_of_floats

             #FH = np.loadtxt(lines, delimiter=' ', skiprows=1, usecols = (1,19,20,21,22,23,24))
             #time = np.loadtxt(lines, dtype='str', delimiter=' ', skiprows=1, usecols = (1))
             #timeUTC timeWVR T0 T1 T2 T3 El Az



        waz_1=data[:,6]
        #print(waz_1)
        waz,fs =  self.findScans(waz_1)


        D = {0:None, 1:None, 2:None, 3:None}
        for ch in range(4):
            D[ch] = self.interpToImage(waz, data[:,1+ch], fs)

        sd = shape(D[0])

        fig, axs = plt.subplots(4, 1, figsize=(12,8))
        for i in range(0,4):
            im=axs[i].imshow(D[i],aspect='auto',interpolation='nearest', origin='lower')
            axs[i].set_xticks(range(0,sd[1],10))
            #axs[i].set_xticklabels('')
            axs[i].set_yticks(range(0,sd[0],60))
            cbar = fig.colorbar(im, extend='both', ax=axs[i])
            cbar.set_label('T_sky[K]\nCh%s'%i)

            if i==3:
                subplots_adjust(hspace=0.01)
                title = filename[-31:].replace('.txt','_atmogram')
                suptitle(title,y=0.97, fontsize=10)
                print("Saving %s.png"%title)
                if not os.path.exists('output_plots/T_atmograms/'):
                    os.makedirs('output_plots/T_atmograms/')
                pl.savefig('output_plots/T_atmograms/'+title+'_BeforeModSub.png')

            if show_plots==1:
                plt.show()
            else:
                pl.close()

        return time, data, D, waz, fs

    def interp_fastdata(self, time, data_w_zeros):
        good_ind=np.where(data_w_zeros > 0.)
        data_wout_zeros=data_w_zeros[good_ind]
        x_axis=time[good_ind]
        f=interp.interp1d(x_axis,data_wout_zeros, kind='cubic', fill_value='extrapolate')
        data_interp=f(time)
        return data_interp

    def calibrate_fast_data(self, interp_data):
        def calibrate_one_ch(data_C, data_A, data_H, data_B):
            calib_ch_data=np.zeros(len(data_A))
            T_hot_arr=np.full(len(data_A), self.T_hot)
            T_cold_arr=np.full(len(data_A), self.T_cold)

            T_ref=(T_hot_arr+T_cold_arr)/2.
            V_ref=(data_H+data_C)/2.
            G=(data_H-data_C)/(T_hot_arr-T_cold_arr)
            V_sky=(data_A+data_B)/2.
            calib_ch_data=T_ref+((V_sky-V_ref)/G)
            return calib_ch_data, G

        def calibrate_one_ch_Gavg(data_C, data_A, data_H, data_B, G_avg_ch):
            calib_ch_data=np.zeros(len(data_A))
            T_hot_arr=np.full(len(data_A), self.T_hot)
            T_cold_arr=np.full(len(data_A), self.T_cold)
            G_avg_ch_arr=np.full(len(data_A), G_avg_ch)

            T_ref=(T_hot_arr+T_cold_arr)/2.
            V_ref=(data_H+data_C)/2.
            V_sky=(data_A+data_B)/2.
            calib_ch_data=T_ref+((V_sky-V_ref)/G_avg_ch_arr)
            return calib_ch_data

        data_calibrated=np.zeros((len(interp_data[:,0]),7))
        data_calibrated_Gavg=np.zeros((len(interp_data[:,0]),7))
        G_array=np.zeros((len(interp_data[:,0]),4))
        data_calibrated[:,0]=interp_data[:,0]
        data_calibrated_Gavg[:,0]=interp_data[:,0]
        data_calibrated[:,5]=interp_data[:,17]
        data_calibrated_Gavg[:,5]=interp_data[:,17]
        data_calibrated[:,6]=interp_data[:,18]
        data_calibrated_Gavg[:,6]=interp_data[:,18]
        for i in range (4):
            data_calibrated[:,1+i],G_array[:,i]=calibrate_one_ch(interp_data[:,4*i+1], interp_data[:,4*i+2], interp_data[:,4*i+3], interp_data[:,4*i+4])
            data_calibrated_Gavg[:,1+i]=calibrate_one_ch_Gavg(interp_data[:,4*i+1], interp_data[:,4*i+2], interp_data[:,4*i+3], interp_data[:,4*i+4], np.mean(G_array[:,i]))

        return data_calibrated, data_calibrated_Gavg, G_array



    def remove_mod(self, calib_data, fn, pathtofn, dTdEl, tilt_par, clean_mod=3, showplots=1, method='fit'): #method=fit/fourier_proj/import_model
    #clean_mod=1--> removes single mod
    #clean_mod=2--> removes double mod
    #clean_mod=3--> removes single and double mod

        if method=='import_model': #removes single and double

            year=fn[:4]
            print('year=', year)
            month=fn[4:6]
            print('month=', month)
            #day=fn[6:8]
            day='01'
            print('day=', day)
            time=fn[9:15]
            print('time=', time)
            hh=fn[9:11]
            print('hh=', hh)
            mm=fn[11:13]
            print('mm=', mm)
            ss=fn[13:15]
            print('ss=', ss)

            fn_date=datetime.datetime(int(year), int(month), int(day), int(hh))

            f = open('mod_param/modulation_parameters_'+year+month+'.txt','rb') #this file has to be created in ModAmptoTilt.py
            p = pk.load(f)
            f.close()

            f = open('mod_param/tilt_angle_fullmonth'+year+month+'.txt','rb') #this file has to be created in ModAmptoTilt.py
            t = pk.load(f)
            f.close()

            for i in range (len(t['x_axis_onemonth'])):
                t_i=t['x_axis_onemonth'][i]
                print('date[i]=', t_i)
                if (t_i.date()==fn_date.date()):
                    if (t_i.hour==fn_date.hour):
                        print('index found.')
                        index=i
                        print('index=', index)
                        print('check_date=',t_i)


            dT0_fn=t['dT0'][index]
            dT1_fn=t['dT1'][index]
            dT2_fn=t['dT2'][index]
            dT3_fn=t['dT3'][index]

            tilt0=p['single'][0]
            tilt1=p['single'][1]
            tilt2=p['single'][2]
            tilt3=p['single'][3]

            tilt=np.mean(np.array(p['single'][:4]))

            print('tilt=', tilt)

            phi_single=p['single'][4]

            A0_fn=tilt*dT0_fn
            A1_fn=tilt*dT1_fn
            A2_fn=tilt*dT2_fn
            A3_fn=tilt*dT3_fn

            p_single=[A0_fn, A1_fn, A2_fn, A3_fn, radians(phi_single)]
            print('p_single=', p_single)
            p_double=p['double']
            p_double[4]=radians(p_double[4])

            print('p_double=', p_double)

            double_mod_removed_data=np.zeros((len(calib_data[:,1]),4))
            single_mod_removed_data=np.zeros((len(calib_data[:,1]),4))
            model_double=np.zeros((len(calib_data[:,1]),4))
            model_single=np.zeros((len(calib_data[:,1]),4))

            T_matrix=calib_data[:,1:5]


            def remove_one_mod(theta_az, T_matrix, mod, p_mod): #mod=1 for Single Mod - mod=2 for Double Mod
                model=np.zeros((len(calib_data[:,1]),4))
                mod_removed_data=np.zeros((len(calib_data[:,1]),4))
                model_old=np.zeros((len(calib_data[:,1]),4))
                mod_removed_data_old=np.zeros((len(calib_data[:,1]),4))

                p=np.zeros(9)
                p[:5]=p_mod

                def modulation(x, a, phi, C):
                    return C + (a * np.sin(np.radians(mod*x) + phi))

                p0=p_mod[0], p_mod[4] #Amp, phase
                p1=p_mod[1], p_mod[4]
                p2=p_mod[2], p_mod[4]
                p3=p_mod[3], p_mod[4]

                model[:,0] = modulation(theta_az, *p0, C=np.mean(T_matrix[:,0]))
                mod_removed_data[:,0] = T_matrix[:,0] - modulation(theta_az, p0[0], p0[1], 0)

                model[:,1] = modulation(theta_az, *p1, C=np.mean(T_matrix[:,1]))
                mod_removed_data[:,1] = T_matrix[:,1] - modulation(theta_az, p1[0], p1[1], 0)

                model[:,2] = modulation(theta_az, *p2, C=np.mean(T_matrix[:,2]))
                mod_removed_data[:,2] = T_matrix[:,2] - modulation(theta_az, p2[0], p2[1], 0)

                model[:,3] = modulation(theta_az, *p3, C=np.mean(T_matrix[:,3]))
                mod_removed_data[:,3] = T_matrix[:,3] - modulation(theta_az, p3[0], p3[1], 0)

                p_err=np.nan

                return p, p_err, mod_removed_data, model

            if clean_mod==3:

                p2, p2_err, double_mod_removed_data, model_double = remove_one_mod(calib_data[:,6], T_matrix, mod=2, p_mod=p_double)
                p1, p1_err, single_mod_removed_data, model_single = remove_one_mod(calib_data[:,6], double_mod_removed_data, mod=1, p_mod=p_single)

            if clean_mod==2:
                p2, p2_err, double_mod_removed_data, model_double = remove_one_mod(calib_data[:,6], T_matrix, mod=2, p_mod=p_double)
                single_mod_removed_data=double_mod_removed_data
                model_single=np.full(np.shape(model_double), 0.)

            if clean_mod==1:
                p1, p1_err, single_mod_removed_data, model_single = remove_one_mod(calib_data[:,6], T_matrix, mod=1, p_mod=p_single)
                double_mod_removed_data=single_mod_removed_data
                model_double=np.full(np.shape(model_single), 0.)



        if method=='fit':

            p1=np.zeros(9)
            p1_err=np.zeros(9)
            p2=np.zeros(9)
            p2_err=np.zeros(9)
            double_mod_removed_data=np.zeros((len(calib_data[:,1]),4))
            single_mod_removed_data=np.zeros((len(calib_data[:,1]),4))
            model_double=np.zeros((len(calib_data[:,1]),4))
            model_single=np.zeros((len(calib_data[:,1]),4))

            T_matrix=calib_data[:,1:5]

            def remove_one_mod(theta_az, T_matrix, mod): #mod=1 for Single Mod - mod=2 for Double Mod
                model=np.zeros((len(calib_data[:,1]),4))
                mod_removed_data=np.zeros((len(calib_data[:,1]),4))
                model_old=np.zeros((len(calib_data[:,1]),4))
                mod_removed_data_old=np.zeros((len(calib_data[:,1]),4))


                #def modulation(x, a, phi, C, mod):
                def modulation(x, a, phi, C, mod_frac):
                    #return C + (a * np.sin(np.radians(mod*x) + phi))
                    return C + (a * np.sin(np.radians(mod_frac*x) + phi))

                def err_modulation(p, x, y):
                    #return modulation(x, p[0], p[1], p[2]) - y
                    return modulation(x, p[0], p[1], p[2], p[3]) - y
                def err_global_modulation(p, x, y1, y2, y3, y4):
                    p1=p[0], p[4], p[5], p[9]
                    p2=p[1], p[4], p[6], p[9]
                    p3=p[2], p[4], p[7], p[9]
                    p4=p[3], p[4], p[8], p[9]

                    err1 = err_modulation(p1, x, y1)
                    err2 = err_modulation(p2, x, y2)
                    err3 = err_modulation(p3, x, y3)
                    err4 = err_modulation(p4, x, y4)

                    return np.concatenate((err1, err2, err3, err4))

                #p=[a1, a2, a3, a4, phi, C1, C2, C3, C4, mod_frac]
                p=[2., 1.5, 1., 0.5, 0., np.mean(T_matrix[:,0]), np.mean(T_matrix[:,1]), np.mean(T_matrix[:,2]), np.mean(T_matrix[:,3]), mod]
                #p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global_modulation, p, args=(theta_az, T_matrix[:,0],  T_matrix[:,1],  T_matrix[:,2],  T_matrix[:,3]), full_output=1)
                res = sp.least_squares(err_global_modulation, p, bounds=((0,0,0,0,-np.inf,0,0,0,0,mod-0.4), (np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,mod+0.4)), args=(theta_az, T_matrix[:,0],  T_matrix[:,1],  T_matrix[:,2],  T_matrix[:,3]))
                print('res=', res)
                p_best=res.x
                print('p_best=', p_best)
                cov_x=res.jac
                p_err=np.sqrt(np.diag(cov_x))

                phi_output=p_best[4]

                p0=p_best[0], p_best[4], p_best[5], p_best[9]
                p1=p_best[1], p_best[4], p_best[6], p_best[9]
                p2=p_best[2], p_best[4], p_best[7], p_best[9]
                p3=p_best[3], p_best[4], p_best[8], p_best[9]


                model[:,0] = modulation(theta_az, *p0)
                mod_removed_data[:,0] = T_matrix[:,0] - modulation(theta_az, p0[0], p0[1], 0, p0[3])

                model[:,1] = modulation(theta_az, *p1)
                mod_removed_data[:,1] = T_matrix[:,1] - modulation(theta_az, p1[0], p1[1], 0, p1[3])

                model[:,2] = modulation(theta_az, *p2)
                mod_removed_data[:,2] = T_matrix[:,2] - modulation(theta_az, p2[0], p2[1], 0, p2[3])

                model[:,3] = modulation(theta_az, *p3)
                mod_removed_data[:,3] = T_matrix[:,3] - modulation(theta_az, p3[0], p3[1], 0, p3[3])


                return p_best, p_err, mod_removed_data, model

            if clean_mod==3:
                p2, p2_err, double_mod_removed_data, model_double = remove_one_mod(calib_data[:,6], T_matrix, mod=2)
                p1, p1_err, single_mod_removed_data, model_single = remove_one_mod(calib_data[:,6], double_mod_removed_data, mod=1)

            if clean_mod==2:
                p2, p2_err, double_mod_removed_data, model_double = remove_one_mod(calib_data[:,6], T_matrix, mod=2)
                p1=np.full(len(p2), np.nan)
                p1_err=np.full(len(p2_err), np.nan)
                single_mod_removed_data=double_mod_removed_data
                model_single=np.full(np.shape(model_double), 0.)

            if clean_mod==1:
                p1, p1_err, single_mod_removed_data, model_single = remove_one_mod(calib_data[:,6], T_matrix, mod=1)
                p2=np.full(len(p1), np.nan)
                p2_err=np.full(len(p1_err), np.nan)
                double_mod_removed_data=single_mod_removed_data
                model_double=np.full(np.shape(model_single), 0.)



        if method=='fit_both': #just makes sense if clean_mod==3
        #removes single and double mod simoultaneously allowing for an extra free parameter l_scale that scales the az axis

            p=np.zeros(15)
            p_err=np.zeros(15)
            mod_removed_data=np.zeros((len(calib_data[:,1]),4))
            model=np.zeros((len(calib_data[:,1]),4))

            T_matrix=calib_data[:,1:5]

            def remove_both_mod(theta_az, T_matrix):
                model=np.zeros((len(calib_data[:,1]),4))
                mod_removed_data=np.zeros((len(calib_data[:,1]),4))

                #def modulation(x, a, phi, C, mod):
                def modulation(x, a, b, phi_a, phi_b, C, l_scale):
                    #return C + (a * np.sin(np.radians(mod*x) + phi))
                    return C + (a * np.sin(np.radians(l_scale*x) + phi_a))+(b * np.sin(np.radians(l_scale*2.*x) + phi_b))

                def err_modulation(p, x, y):
                    #return modulation(x, p[0], p[1], p[2]) - y
                    return modulation(x, p[0], p[1], p[2], p[3], p[4], p[5]) - y
                def err_global_modulation(p, x, y1, y2, y3, y4):
                    p1=p[0], p[4], p[8], p[9], p[10], p[14]
                    p2=p[1], p[5], p[8], p[9], p[11], p[14]
                    p3=p[2], p[6], p[8], p[9], p[12], p[14]
                    p4=p[3], p[7], p[8], p[9], p[13], p[14]

                    err1 = err_modulation(p1, x, y1)
                    err2 = err_modulation(p2, x, y2)
                    err3 = err_modulation(p3, x, y3)
                    err4 = err_modulation(p4, x, y4)

                    return np.concatenate((err1, err2, err3, err4))

                #p=[a1, a2, a3, a4, b1, b2, b3, b4, phi_a, phi_b, C1, C2, C3, C4, l_scale]
                p=[2., 1.5, 1., 0.5, 2., 1.5, 1., 0.5, 0., 0., np.mean(T_matrix[:,0]), np.mean(T_matrix[:,1]), np.mean(T_matrix[:,2]), np.mean(T_matrix[:,3]), 1]
                #p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global_modulation, p, args=(theta_az, T_matrix[:,0],  T_matrix[:,1],  T_matrix[:,2],  T_matrix[:,3]), full_output=1)
                res = sp.least_squares(err_global_modulation, p,bounds=((0,0,0,0,0,0,0,0,-np.inf,-np.inf,0,0,0,0,0), (np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf, np.inf,np.inf,np.inf,np.inf,np.inf,np.inf)), args=(theta_az, T_matrix[:,0],  T_matrix[:,1],  T_matrix[:,2],  T_matrix[:,3]))
                print('res=', res)
                p_best=res.x
                print('p_best=', p_best)
                cov_x=res.jac
                p_err=np.sqrt(np.diag(cov_x))


                p0=p_best[0], p_best[4], p_best[8], p_best[9], p_best[10], p_best[14]
                p1=p_best[1], p_best[5], p_best[8], p_best[9], p_best[11], p_best[14]
                p2=p_best[2], p_best[6], p_best[8], p_best[9], p_best[12], p_best[14]
                p3=p_best[3], p_best[7], p_best[8], p_best[9], p_best[13], p_best[14]


                p0_zerooffs=p_best[0], p_best[4], p_best[8], p_best[9], 0, p_best[14]
                p1_zerooffs=p_best[1], p_best[5], p_best[8], p_best[9], 0, p_best[14]
                p2_zerooffs=p_best[2], p_best[6], p_best[8], p_best[9], 0, p_best[14]
                p3_zerooffs=p_best[3], p_best[7], p_best[8], p_best[9], 0, p_best[14]



                model[:,0] = modulation(theta_az, *p0)
                mod_removed_data[:,0] = T_matrix[:,0] - modulation(theta_az, *p0_zerooffs)#[0], p0[1], 0, p0[3])

                model[:,1] = modulation(theta_az, *p1)
                mod_removed_data[:,1] = T_matrix[:,1] - modulation(theta_az, *p1_zerooffs)#[0], p1[1], 0, p1[3])

                model[:,2] = modulation(theta_az, *p2)
                mod_removed_data[:,2] = T_matrix[:,2] - modulation(theta_az, *p2_zerooffs)#[0], p2[1], 0, p2[3])

                model[:,3] = modulation(theta_az, *p3)
                mod_removed_data[:,3] = T_matrix[:,3] - modulation(theta_az, *p3_zerooffs)#[0], p3[1], 0, p3[3])


                if clean_mod==3:
                    p2, p2_err, double_mod_removed_data, model_double = np.full(np.shape(p_best),np.nan), np.full(np.shape(p_err),np.nan), np.full(np.shape(mod_removed_data),np.nan), np.full(np.shape(model),np.nan)
                    p1, p1_err, single_mod_removed_data, model_single = p_best, p_err, mod_removed_data, model #single and double mod removed data
                else:
                    print('Error: Method fit_both just available for clean_mod=3.')
                    sys.exit()


        if method == 'fit_both':
            #Both Mods plotted at once
            fig=plt.figure(figsize=(12,8))
            gs = GridSpec(4, 2, width_ratios=[1, 1], height_ratios=[3, 1, 3, 1], hspace=0.5)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = fig.add_subplot(gs[3])
            ax5 = fig.add_subplot(gs[4])
            ax6 = fig.add_subplot(gs[5])
            ax7 = fig.add_subplot(gs[6])
            ax8 = fig.add_subplot(gs[7])

            ax1.scatter(calib_data[:,6], T_matrix[:,0], s=2, c='c', label="Calibrated Data") #T[ch]=calib_data[:, ch+1]
            ax1.scatter(calib_data[:,6], model[:,0], s=1, c='k', label="Mod Fit\nl_fit="+str(round(p[14],4)))
            #ax1.scatter(calib_data[:,6], model_double_old[:,0], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[0],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[5],2)))
            ax1.legend(loc='upper right')
            ax1.set_title('Channel 0')

            ax3.scatter(calib_data[:,6], mod_removed_data[:,0], s=2, c='r', label="Residuals")
            ax3.legend(loc='upper right')

            ax2.scatter(calib_data[:,6], T_matrix[:,1], s=2, c='c', label="Calibrated Data")
            ax2.scatter(calib_data[:,6], model[:,1], s=1, c='k', label="Mod Fit\nl_fit="+str(round(p[14],4)))
            #ax2.scatter(calib_data[:,6], model_double_old[:,1], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[1],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[6],2)))
            ax2.legend(loc='upper right')
            ax2.set_title('Channel 1')

            ax4.scatter(calib_data[:,6], mod_removed_data[:,1], s=2, c='r', label="Residuals")
            ax4.legend(loc='upper right')

            ax5.scatter(calib_data[:,6], calib_data[:,3], s=2, c='c', label="Calibrated Data")
            ax5.scatter(calib_data[:,6], model[:,2], s=1, c='k', label="Mod Fit\nl_fit"+str(round(p[14],4)))
            #ax5.scatter(calib_data[:,6], model_double_old[:,2], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[2],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[7],2)))
            ax5.legend(loc='upper right')
            ax5.set_title('Channel 2')

            ax7.scatter(calib_data[:,6], mod_removed_data[:,2], s=2, c='r', label="Residuals")
            ax7.legend(loc='upper right')

            ax6.scatter(calib_data[:,6], calib_data[:,4], s=2, c='c', label="Calibrated Data")
            ax6.scatter(calib_data[:,6], model[:,3], s=1, c='k', label="Mod Fit\nl_fit="+str(round(p[14],4)))
            #ax6.scatter(calib_data[:,6], model_double_old[:,3], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[3],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[8],2)))
            ax6.legend(loc='upper right')
            ax6.set_title('Channel 3')

            ax8.scatter(calib_data[:,6], mod_removed_data[:,3], s=2, c='r', label="Residuals")
            ax8.legend(loc='upper right')

            pl.suptitle('Single and Double Modulation Fit and Subtraction\nFitted Function: C + A⋅sin(2⋅l_scale⋅theta_az + phiA)+ B⋅sin(l_scale⋅theta_az + phiB)\n'+fn[:-4])
            if not os.path.exists('output_plots/T_atmograms/'):
                os.makedirs('output_plots/T_atmograms/')
            pl.savefig('output_plots/T_atmograms/'+fn[-31:-4]+'_ModFit.png')

            if showplots==1:
                pl.show()
            else:
                pl.close()


            #Double Mod
            fig=plt.figure(figsize=(12,8))
            gs = GridSpec(4, 2, width_ratios=[1, 1], height_ratios=[3, 1, 3, 1], hspace=0.5)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = fig.add_subplot(gs[3])
            ax5 = fig.add_subplot(gs[4])
            ax6 = fig.add_subplot(gs[5])
            ax7 = fig.add_subplot(gs[6])
            ax8 = fig.add_subplot(gs[7])

            ax1.scatter(calib_data[:,6], T_matrix[:,0], s=2, c='c', label="Calibrated Data") #T[ch]=calib_data[:, ch+1]
            ax1.scatter(calib_data[:,6], model_double[:,0], s=1, c='k', label="Double Mod Fit\nA[K]="+str(round(p2[0],2))+'\nPhi[rad]='+str(round(p2[4],2))+'\nC[K]='+str(round(p2[5],2)))
            #ax1.scatter(calib_data[:,6], model_double_old[:,0], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[0],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[5],2)))
            ax1.legend(loc='upper right')
            ax1.set_title('Channel 0')

            ax3.scatter(calib_data[:,6], double_mod_removed_data[:,0], s=2, c='r', label="Residuals")
            ax3.legend(loc='upper right')

            ax2.scatter(calib_data[:,6], T_matrix[:,1], s=2, c='c', label="Calibrated Data")
            ax2.scatter(calib_data[:,6], model_double[:,1], s=1, c='k', label="Double Mod Fit\nA[K]="+str(round(p2[1],2))+'\nPhi[rad]='+str(round(p2[4],2))+'\nC[K]='+str(round(p2[6],2)))
            #ax2.scatter(calib_data[:,6], model_double_old[:,1], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[1],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[6],2)))
            ax2.legend(loc='upper right')
            ax2.set_title('Channel 1')

            ax4.scatter(calib_data[:,6], double_mod_removed_data[:,1], s=2, c='r', label="Residuals")
            ax4.legend(loc='upper right')

            ax5.scatter(calib_data[:,6], calib_data[:,3], s=2, c='c', label="Calibrated Data")
            ax5.scatter(calib_data[:,6], model_double[:,2], s=1, c='k', label="Double Mod Fit\nA[K]="+str(round(p2[2],2))+'\nPhi[rad]='+str(round(p2[4],2))+'\nC[K]='+str(round(p2[7],2)))
            #ax5.scatter(calib_data[:,6], model_double_old[:,2], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[2],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[7],2)))
            ax5.legend(loc='upper right')
            ax5.set_title('Channel 2')

            ax7.scatter(calib_data[:,6], double_mod_removed_data[:,2], s=2, c='r', label="Residuals")
            ax7.legend(loc='upper right')

            ax6.scatter(calib_data[:,6], calib_data[:,4], s=2, c='c', label="Calibrated Data")
            ax6.scatter(calib_data[:,6], model_double[:,3], s=1, c='k', label="Double Mod Fit\nA[K]="+str(round(p2[3],2))+'\nPhi[rad]='+str(round(p2[4],2))+'\nC[K]='+str(round(p2[8],2)))
            #ax6.scatter(calib_data[:,6], model_double_old[:,3], s=1, c='r', label="Double Mod Fit\nA[K]="+str(round(p2_old[3],2))+'\nPhi[rad]='+str(round(p2_old[4],2))+'\nC[K]='+str(round(p2_old[8],2)))
            ax6.legend(loc='upper right')
            ax6.set_title('Channel 3')

            ax8.scatter(calib_data[:,6], double_mod_removed_data[:,3], s=2, c='r', label="Residuals")
            ax8.legend(loc='upper right')

            pl.suptitle('Double Modulation Fit and Subtraction\nFitted Function: C + A⋅sin(2⋅theta_az + phi)\n'+fn[:-4])
            if not os.path.exists('output_plots/T_atmograms/'):
                os.makedirs('output_plots/T_atmograms/')
            pl.savefig('output_plots/T_atmograms/'+fn[-31:-4]+'_DoubleModFit.png')

            if showplots==1:
                pl.show()
            else:
                pl.close()

            #Single Mod

            #from El el_correction
            #(tilt, phi)=tilt_par
            #dT0=dTdEl[0]*tilt*(np.sin(np.radians(np.asarray(calib_data[:,6])+phi)))
            #dT1=dTdEl[1]*tilt*(np.sin(np.radians(np.asarray(calib_data[:,6])+phi)))
            #dT2=dTdEl[2]*tilt*(np.sin(np.radians(np.asarray(calib_data[:,6])+phi)))
            #dT3=dTdEl[3]*tilt*(np.sin(np.radians(np.asarray(calib_data[:,6])+phi)))

            fig=plt.figure(figsize=(12,8))
            gs = GridSpec(6, 2, width_ratios=[1, 1], height_ratios=[3, 1, 1, 3, 1, 1], hspace=0.5)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = fig.add_subplot(gs[3])
            #ax5 = fig.add_subplot(gs[4])
            #ax6 = fig.add_subplot(gs[5])
            ax7 = fig.add_subplot(gs[6])
            ax8 = fig.add_subplot(gs[7])
            ax9 = fig.add_subplot(gs[8])
            ax10 = fig.add_subplot(gs[9])
            #ax11 = fig.add_subplot(gs[10])
            #ax12 = fig.add_subplot(gs[11])

            phi_fit=np.degrees(p1[4])-360.

            ax1.scatter(calib_data[:,6], double_mod_removed_data[:,0], s=2, c='c', label="Double mod removed Data") #T[ch]=calib_data[:, ch+1]
            ax1.scatter(calib_data[:,6], model_single[:,0], s=1, c='k', label="Single Mod Fit\nA="+str(round(p1[0],2))+'\nPhi='+str(round(phi_fit,2))+'\nC='+str(round(p1[5],2))+'\nl_fit='+str(round(p1[9],2)))
            #ax1.scatter(calib_data[:,6], p1[5]+dT0, s=1, c='r', label="Single Mod from El corr\nA="+str(round(dTdEl[0]*tilt, 2))+'\nPhi='+str(round(phi,2)))
            ax1.legend(loc='upper right')
            ax1.set_title('Channel 0')

            ax3.scatter(calib_data[:,6], single_mod_removed_data[:,0], s=1, c='r', label="Residuals from fit")
            ax3.legend(loc='upper right')

            #ax5.scatter(calib_data[:,6], double_mod_removed_data[:,0]-dT0, s=1, c='r', label="Residuals from El correction")
            #ax5.legend(loc='upper right')

            ax2.scatter(calib_data[:,6], double_mod_removed_data[:,1], s=2, c='c', label="Double mod removed Data")
            ax2.scatter(calib_data[:,6], model_single[:,1], s=1, c='k', label="Single Mod Fit\nA="+str(round(p1[1],2))+'\nPhi='+str(round(phi_fit,2))+'\nC='+str(round(p1[6],2))+'\nl_fit='+str(round(p1[9],2)))
            #ax2.scatter(calib_data[:,6], p1[6]+dT1, s=1, c='r', label="Single Mod from El corr\nA="+str(round(dTdEl[1]*tilt,2))+'\nPhi='+str(round(phi,2)))
            ax2.legend(loc='upper right')
            ax2.set_title('Channel 1')

            ax4.scatter(calib_data[:,6], single_mod_removed_data[:,1], s=1, c='r', label="Residuals from fit")
            ax4.legend(loc='upper right')

            #ax6.scatter(calib_data[:,6], double_mod_removed_data[:,1]-dT1, s=1, c='r', label="Residuals from El correction")
            #ax6.legend(loc='upper right')

            ax7.scatter(calib_data[:,6], double_mod_removed_data[:,2], s=2, c='c', label="Double mod removed Data")
            ax7.scatter(calib_data[:,6], model_single[:,2], s=1, c='k', label="Single Mod Fit\nA="+str(round(p1[2],2))+'\nPhi='+str(round(phi_fit,2))+'\nC='+str(round(p1[7],2))+'\nl_fit='+str(round(p1[9],2)))
            #ax7.scatter(calib_data[:,6], p1[7]+dT2, s=1, c='r', label="Single Mod from El corr\nA="+str(round(dTdEl[2]*tilt,2))+'\nPhi='+str(round(phi,2)))
            ax7.legend(loc='upper right')
            ax7.set_title('Channel 2')

            ax9.scatter(calib_data[:,6], single_mod_removed_data[:,2], s=1, c='r', label="Residuals from fit")
            ax9.legend(loc='upper right')

            #ax11.scatter(calib_data[:,6], double_mod_removed_data[:,2]-dT2, s=1, c='r', label="Residuals from El correction")
            #ax11.legend(loc='upper right')

            ax8.scatter(calib_data[:,6], double_mod_removed_data[:,3], s=2, c='c', label="Double mod removed Data")
            ax8.scatter(calib_data[:,6], model_single[:,3], s=1, c='k', label="Single Mod Fit\nA="+str(round(p1[3],2))+'\nPhi='+str(round(phi_fit,2))+'\nC='+str(round(p1[8],2))+'\nl_fit='+str(round(p1[9],2)))
            #ax8.scatter(calib_data[:,6], p1[8]+dT3, s=1, c='r', label="Single Mod from El corr\nA="+str(round(dTdEl[3]*tilt,2))+'\nPhi='+str(round(phi,2)))
            ax8.legend(loc='upper right')
            ax8.set_title('Channel 3')

            ax10.scatter(calib_data[:,6], single_mod_removed_data[:,3], s=1, c='r', label="Residuals from fit")
            ax10.legend(loc='upper right')

            #ax12.scatter(calib_data[:,6], double_mod_removed_data[:,3]-dT3, s=1, c='r', label="Residuals from El correction")
            #ax12.legend(loc='upper right')

            pl.suptitle('Single Modulation Fit and Subtraction\nFitted Function: C + A⋅sin(theta_az + phi)\n'+fn[:-4])
            if not os.path.exists('output_plots/T_atmograms/'):
                os.makedirs('output_plots/T_atmograms/')
            pl.savefig('output_plots/T_atmograms/'+fn[-31:-4]+'_SingleModFit'+method+'.png')
            if showplots==1:
                pl.show()
            else:
                pl.close()


        return p2, p2_err, p1, p1_err, double_mod_removed_data, single_mod_removed_data, model_double, model_single



    def read_Az_fast(self, filename, pathtofn='', dTdEl=0, tilt_par=(0.,0.), clean_mod=3, test_plot_int=1, test_plot_cal=1, test_plot_single_mod=1, clean_method='fit', show_mod_plots=0):
    #this function read fast files and returns as an array the
    #interpolated values for [t, T0_A, T1_A, T2_A, T3_A, El, Az]
    #plus it writes them into a file *_interp.txt

    #clean_mod=1--> removes single mod
    #clean_mod=2--> removes double mod
    #clean_mod=3--> removes single and double mod
        # with open('wvr1_data/'+pathtofn+filename) as f:
        #      lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
        #      FH = np.loadtxt(lines, dtype='str', delimiter=' ', skiprows=1, usecols = (0,1,19,20,21,22,23,24))
        #      time=FH[:,0]

        with open(pathtofn+filename[:-9]+'/'+filename) as f:
            lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
            time = np.loadtxt(lines, dtype='str', delimiter=' ', skiprows=1, usecols = (0))
            #The data order is: TIME [TIMEWVR CH0C CH0A CH0H CH0B CH1C CH1A CH1H CH1B CH2C CH2A CH2H CH2B CH3C CH3A CH3H CH3B EL AZ]
            #time=FH_all[:,0]

        with open(pathtofn+filename[:-9]+'/'+filename) as f:
            lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
            FH = np.loadtxt(lines, delimiter=' ', skiprows=1, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
            #The data order is: TIMEWVR CH0C CH0A CH0H CH0B CH1C CH1A CH1H CH1B CH2C CH2A CH2H CH2B CH3C CH3A CH3H CH3B EL AZ

        interp_data=np.zeros((len(FH[:,0]),19))
        interp_data[:,0] = FH[:,0]
        interp_data[:,17] =FH[:,17]
        interp_data[:,18] = FH[:,18]%(360.)
        waz_1=FH[:,18] #wrapped Azimuth

        plt.plot(FH[:,0], interp_data[:,18])
        plt.xlabel('time[s]')
        plt.ylabel('Az[deg]')
        plt.title('Az Scan')
        #plt.show()
        plt.close()

        for i in range(len(interp_data[0,:])-3):
            interp_data[:,i+1]=self.interp_fastdata(FH[:,0], FH[:,1+i])

        if test_plot_int==1:
            #fig, axes = plt.subplots(2, 2)
            fig=plt.figure()
            gs = GridSpec(4, 2, width_ratios=[1, 1], height_ratios=[3, 1, 3, 1], hspace=0.5)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = fig.add_subplot(gs[3])
            ax5 = fig.add_subplot(gs[4])
            ax6 = fig.add_subplot(gs[5])
            ax7 = fig.add_subplot(gs[6])
            ax8 = fig.add_subplot(gs[7])

            T0=np.concatenate((interp_data[:,2], interp_data[:,4]), axis=None)
            el_fit=np.concatenate((interp_data[:,18], interp_data[:,18]), axis=None)
            T0_coef=np.polyfit(el_fit, T0, 10)
            T0_fit=np.poly1d(T0_coef)
            T0_plot=T0_fit(interp_data[:,18])

            ax1.scatter(interp_data[:,18], interp_data[:,2], s=2, c='r', label="T_A")
            ax1.scatter(interp_data[:,18], interp_data[:,4], s=2, c='c', label="T_B")
            ax1.plot(interp_data[:,18], T0_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax1.legend(loc='upper right')
            ax1.set_title('Channel 0')

            ax3.scatter(interp_data[:,18], interp_data[:,2]-T0_plot, s=2, c='r', label="T_A - PolyFit")
            ax3.scatter(interp_data[:,18], interp_data[:,4]-T0_plot, s=2, c='c', label="T_B - PolyFit")
            ax3.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax3.legend(loc='upper right')

            T1=np.concatenate((interp_data[:,6], interp_data[:,8]), axis=None)
            T1_coef=np.polyfit(el_fit, T1, 10)
            T1_fit=np.poly1d(T1_coef)
            T1_plot=T1_fit(interp_data[:,18])

            ax2.scatter(interp_data[:,18], interp_data[:,6], s=2, c='r', label="T_A")
            ax2.scatter(interp_data[:,18], interp_data[:,8], c='c', s=2,label="T_B")
            ax2.plot(interp_data[:,18], T1_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax2.legend(loc='upper right')
            ax2.set_title('Channel 1')

            ax4.scatter(interp_data[:,18], interp_data[:,6]-T1_plot, s=2, c='r', label="T_A - PolyFit")
            ax4.scatter(interp_data[:,18], interp_data[:,8]-T1_plot, c='c', s=2,label="T_B - PolyFit")
            ax4.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax4.legend(loc='upper right')

            T2=np.concatenate((interp_data[:,10], interp_data[:,12]), axis=None)
            T2_coef=np.polyfit(el_fit, T2, 10)
            T2_fit=np.poly1d(T2_coef)
            T2_plot=T2_fit(interp_data[:,18])

            ax5.scatter(interp_data[:,18], interp_data[:,10], c='r',s=2,label="T_A")
            ax5.scatter(interp_data[:,18], interp_data[:,12], c='c',s=2, label="T_B")
            ax5.plot(interp_data[:,18], T2_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax5.legend(loc='upper right')
            ax5.set_title('Channel 2')

            ax7.scatter(interp_data[:,18], interp_data[:,10]-T2_plot, c='r', s=2, label="T_A - PolyFit")
            ax7.scatter(interp_data[:,18], interp_data[:,12]-T2_plot, c='c', s=2, label="T_B - PolyFit")
            ax7.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax7.legend(loc='upper right')

            T3=np.concatenate((interp_data[:,14], interp_data[:,16]), axis=None)
            T3_coef=np.polyfit(el_fit, T3, 10)
            T3_fit=np.poly1d(T3_coef)
            T3_plot=T3_fit(interp_data[:,18])

            ax6.scatter(interp_data[:,18], interp_data[:,14], c='r', s=2, label="T_A")
            ax6.scatter(interp_data[:,18], interp_data[:,16], c='c', s=2, label="T_B")
            ax6.plot(interp_data[:,18], T3_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax6.legend(loc='upper right')
            ax6.set_title('Channel 3')

            ax8.scatter(interp_data[:,18], interp_data[:,14]-T3_plot, c='r', s=2, label="T_A - PolyFit")
            ax8.scatter(interp_data[:,18], interp_data[:,16]-T3_plot, c='c', s=2, label="T_B - PolyFit")
            ax8.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax8.legend(loc='upper right')

            pl.suptitle('Raw Temperature readings at mirror position A and B\n'+filename[:-4])
            if not os.path.exists(f'../Output/'+filename[:-4]):
                os.makedirs(f'../Output/'+filename[:-4]+'/')
            pl.savefig('../Output/'+filename[:-4]+'/Backlash-corrected_data.png')
            #pl.show()
            pl.close()


        calib_data, calib_data_Gavg, G = self.calibrate_fast_data(interp_data)

        if test_plot_cal==1:
            fig = pl.figure()
            ax1 = pl.subplot2grid((3,1), (0,0), rowspan=2)
            ax2 = pl.subplot2grid((3,1), (2,0))
            ax1.scatter(calib_data[:,6], calib_data[:,1], s=2, c='k')
            ax1.set_title('Calibrated Temperature')
            ax2.scatter(calib_data[:,6], G[:,0], s=2, c='r')
            ax2.set_title('Gain')
            # pl.plot(calib_data[:,0], calib_data[:,2], label="T_src-Ch1")
            # pl.plot(calib_data[:,0], calib_data[:,3], label="T_src-Ch2")
            # pl.plot(calib_data[:,0], calib_data[:,4], label="T_src-Ch3")
            pl.suptitle('Calibrated T0 data\n'+filename[:-4])
            #pl.show()
            pl.close()

            G_avg_ch=np.full(len(G[:,0]), np.mean(G[:,0]))
            fig = pl.figure()
            ax1 = pl.subplot2grid((3,1), (0,0), rowspan=2)
            ax2 = pl.subplot2grid((3,1), (2,0))
            ax1.scatter(calib_data_Gavg[:,6], calib_data_Gavg[:,1], s=2, c='k')
            ax1.set_title('Calibrated Temperature G_avg')
            ax2.scatter(calib_data_Gavg[:,6], G_avg_ch, s=2, c='r')
            ax2.set_title('Gain')
            # pl.plot(calib_data[:,0], calib_data[:,2], label="T_src-Ch1")
            # pl.plot(calib_data[:,0], calib_data[:,3], label="T_src-Ch2")
            # pl.plot(calib_data[:,0], calib_data[:,4], label="T_src-Ch3")
            pl.suptitle('Calibrated T0 data Averaging Gain\n'+filename[:-4])
            #pl.show()
            pl.close()


        mat = np.matrix(interp_data)
        with open(filename[:-4]+'_interp.txt','wb') as f:
            for line in mat:
                np.savetxt(f, line, fmt='%.5f')
        f.close()

        mat_c = np.matrix(calib_data)
        with open(filename[:-4]+'_calib.txt','wb') as f:
            for line in mat_c:
                np.savetxt(f, line, fmt='%.5f')
        f.close()

        #Create atmogram before removing mod
        waz,fs =  self.findScans(waz_1)

        D = {0:None, 1:None, 2:None, 3:None}
        for ch in range(4):
            D[ch] = self.interpToImage(waz, calib_data_Gavg[:,1+ch], fs)

        sd = shape(D[0])

        fig, axs = plt.subplots(4, 1, figsize=(12,8))
        for i in range(0,4):
            im=axs[i].imshow(D[i],aspect='auto',interpolation='nearest', origin='lower')
            axs[i].set_xticks(range(0,sd[1],10))
            #axs[i].set_xticklabels('')
            axs[i].set_yticks(range(0,sd[0],60))
            cbar = fig.colorbar(im, extend='both', ax=axs[i])
            cbar.set_label('T_sky[K]\nCh%s'%i)

            if i==3:
                subplots_adjust(hspace=0.01)
                title = filename[-31:].replace('.txt','_atmogram')
                print(title)
                suptitle(title+'\nBefore mod subtraction',y=0.97, fontsize=10)
                print("Saving %s.png"%title)
                if not os.path.exists('output_plots/T_atmograms/'):
                    os.makedirs('output_plots/T_atmograms/')
                pl.savefig(self.path_to_all+'output_plots/T_atmograms/'+title+'_BeforeModSub.png')
                if show_mod_plots==1:
                    plt.show()
                else:
                    plt.close()

        if clean_mod>=1:

            p2, p2_err, p1, p1_err, double_mod_removed_data, single_mod_removed_data, model_double, model_single = self.remove_mod(calib_data_Gavg, filename, pathtofn, dTdEl, tilt_par, clean_mod, showplots=show_mod_plots, method=clean_method)
        #
            #return waz_1, p, p_err, mod_removed_data, model, calib_data_Gavg
            D_clean = {0:None, 1:None, 2:None, 3:None}
            D_double = {0:None, 1:None, 2:None, 3:None}
            D_single = {0:None, 1:None, 2:None, 3:None}
            for ch in range(4):
                D_clean[ch] = self.interpToImage(waz, single_mod_removed_data[:,ch], fs)
                D_double[ch] = self.interpToImage(waz, model_double[:,ch], fs)
                D_single[ch] = self.interpToImage(waz, model_single[:,ch], fs)
            sd = shape(D[0])

            fig, axs = plt.subplots(4, 1, figsize=(8,10))
            az_offs = 250.04 #deg
            im=axs[0].imshow(D[0], aspect='auto',interpolation='nearest', origin='lower')
            axs[0].set_xticks(range(0,sd[1],10))
            axs[0].set_ylabel('Az[deg]', fontsize=12)
            axs[0].set_yticklabels([str(int(i-az_offs)) for i in range(0,sd[1],10)], fontsize=12)
            axs[0].set_xticklabels('')
            axs[0].set_yticks(range(0,sd[0],60))
            axs[0].set_title('T measured')
            cbar = fig.colorbar(im, extend='both', ax=axs[0])
            cbar.set_label('T[K]', fontsize=12)

            im=axs[1].imshow(D_single[0], aspect='auto',interpolation='nearest', origin='lower')
            axs[1].set_xticks(range(0,sd[1],10))
            axs[1].set_ylabel('Az[deg]', fontsize=12)
            axs[1].set_yticklabels([str(int(i-az_offs)) for i in range(0,sd[1],10)], fontsize=12)
            axs[1].set_yticks(range(0,sd[0],60))
            axs[1].set_xticklabels('')
            axs[1].set_title('Dipole model')
            cbar = fig.colorbar(im, extend='both', ax=axs[1])
            cbar.set_label('T[K]', fontsize=12)

            im=axs[2].imshow(D_double[0], aspect='auto',interpolation='nearest', origin='lower')
            axs[2].set_xticks(range(0,sd[1],10))
            axs[2].set_ylabel('Az[deg]', fontsize=12)
            axs[2].set_yticklabels([str(int(i-az_offs)) for i in range(0,sd[1],10)], fontsize=12)
            axs[2].set_yticks(range(0,sd[0],60))
            axs[2].set_xticklabels('')
            axs[2].set_title('Quadrupole Model')
            cbar = fig.colorbar(im, extend='both', ax=axs[2])
            cbar.set_label('T[K]', fontsize=12)

            im=axs[3].imshow(D_clean[0], aspect='auto',interpolation='nearest', origin='lower')
            axs[3].set_xticks(range(0,sd[1],10))
            axs[3].set_yticklabels([str(int(i-az_offs)) for i in range(0,sd[1],10)])
            axs[3].set_yticks(range(0,sd[0],60))
            axs[3].set_ylabel('Az[deg]', fontsize=12)
            axs[3].set_xlabel('Scan Number', fontsize=12)
            axs[3].set_title('T - (Dipole+Quadrupole)')
            cbar = fig.colorbar(im, extend='both', ax=axs[3])
            cbar.set_label('T[K]', fontsize=12)

            subplots_adjust(hspace=0.2)
            title = filename[-31:].replace('_scanAz_fast.txt','')
            #suptitle(title+'\nChannel 0',y=0.97, fontsize=18)
            if not os.path.exists('output_plots/T_atmograms/'):
                os.makedirs('output_plots/T_atmograms/')
            pl.savefig(self.path_to_all+'output_plots/T_atmograms/'+title+'_AfterModSub_'+clean_method+'_allsteps_Tch0.png')
            plt.close()


            fig, axs = plt.subplots(4, 1, figsize=(12,8))
            for i in range(0,4):
                im=axs[i].imshow(D_clean[i], aspect='auto',interpolation='nearest', origin='lower')
                axs[i].set_xticks(range(0,sd[1],10))
                #axs[i].set_xticklabels('')
                axs[i].set_yticks(range(0,sd[0],60))
                #axs[i].set_title('TSRC%s'%i)
                cbar = fig.colorbar(im, extend='both', ax=axs[i])
                cbar.set_label('T_sky[K]\nCh%s'%i)
                if i==3:
                    subplots_adjust(hspace=0.01)
                    title = filename[-31:].replace('.txt','_atmogram_modclean')
                    suptitle(title+'\nAfter mod subtraction',y=0.97, fontsize=10)
                    print("Saving %s.png"%title)
                    if not os.path.exists('output_plots/T_atmograms/'):
                        os.makedirs('output_plots/T_atmograms/')
                    pl.savefig(self.path_to_all+'output_plots/T_atmograms/'+title+'_AfterModSub_'+clean_method+'.png')
                    if show_mod_plots==1:
                        plt.show()
                    else:
                        plt.close()

            fig, axs = plt.subplots(4, 1, figsize=(12,8))
            for i in range(0,4):
                im=axs[i].imshow(D_double[i], aspect='auto',interpolation='nearest', origin='lower')
                axs[i].set_xticks(range(0,sd[1],10))
                #axs[i].set_xticklabels('')
                axs[i].set_yticks(range(0,sd[0],60))
                #axs[i].set_title('TSRC%s'%i)
                cbar = fig.colorbar(im, extend='both', ax=axs[i])
                cbar.set_label('T_sky[K]\nCh%s'%i)
                if i==3:
                    subplots_adjust(hspace=0.01)
                    title = filename[-31:].replace('.txt','_atmogram_doublemod')
                    suptitle(title+'\nDouble modulation model',y=0.97, fontsize=10)
                    print("Saving %s.png"%title)
                    if not os.path.exists('output_plots/T_atmograms/'):
                        os.makedirs('output_plots/T_atmograms/')
                    pl.savefig(self.path_to_all+'output_plots/T_atmograms/'+title+'_'+clean_method+'.png')
                    if show_mod_plots==1:
                        plt.show()
                    else:
                        plt.close()

            fig, axs = plt.subplots(4, 1, figsize=(12,8))
            for i in range(0,4):
                im=axs[i].imshow(D_single[i], aspect='auto',interpolation='nearest', origin='lower')
                axs[i].set_xticks(range(0,sd[1],10))
                #axs[i].set_xticklabels('')
                axs[i].set_yticks(range(0,sd[0],60))
                #axs[i].set_title('TSRC%s'%i)
                cbar = fig.colorbar(im, extend='both', ax=axs[i])
                cbar.set_label('T_sky[K]\nCh%s'%i)
                if i==3:
                    subplots_adjust(hspace=0.01)
                    title = filename[-31:].replace('.txt','_atmogram_singlemod')
                    suptitle(title+'\nSingle modulation model',y=0.97, fontsize=10)
                    if not os.path.exists('output_plots/T_atmograms/'):
                        os.makedirs('output_plots/T_atmograms/')
                    pl.savefig(self.path_to_all+'output_plots/T_atmograms/'+title+'_'+clean_method+'.png')
                    #savefig(title+'.png')
                    if show_mod_plots==1:
                        plt.show()
                    else:
                        plt.close()

            print('Returning single/double modulation clear data.')
            return D_clean, waz_1, single_mod_removed_data, p2, p2_err, p1, p1_err, calib_data_Gavg, FH

        elif clean_mod==0:
            print('Returning raw calibrated data.')
            return time, D, waz, calib_data_Gavg, FH

        else:
            print('remove_mod parameter value is not acceptable. \nChoose 0 if you want raw calibrated data.\nChoose 1 if you want single and double modulation removed data. ')

    def read_Az(self, filename, pathtofn=''):
        if filename[-8:-4]=='slow':
            FH=self.read_Az_slow(filename)
        elif filename[-8:-4]=='fast':
            FH=self.read_Az_fast(filename, pathtofn)
        else:
            print("ERROR in the given filename.")
            FH=0
        return FH



    def read_pwvatmo(self, fn, az_lim='None', template='SPole_annual_50', savefig_fn='None', savefig_fn_posting='None', show=0):

            def plot_atmo(waz, pwv_list, title=''):
                az, fs =  self.findScans(waz)
                idx=np.where(waz<=360)
                #pl.scatter(az, pwv_list)
                #pl.show()
                D_pwv = self.interpToImage(az, pwv_list, fs)
                print('pwv_list_shape=', np.shape(pwv_list))
                print('az_shape=', np.shape(az))
                fig=pl.figure(figsize=(10,6))
                #to match BK atmo
                D_pwv=D_pwv#[:,:66]
                im=pl.imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
                if not az_lim == 'None':
                    (az_min, az_max) = az_lim
                    pl.ylim(az_min, az_max)
                pl.suptitle('PWV Atmogram\n'+fn[:-4])
                pl.title('PWV avg ='+str(round(np.nanmean(D_pwv),2)))
                pl.xlabel('Scan Number')
                pl.ylabel('Az[deg]')
                cbar = fig.colorbar(im, extend='both')
                cbar.set_label('PWV[um]')
                if not savefig_fn == 'None':
                    pl.savefig(savefig_fn)
                if not savefig_fn_posting=='None':
                    pl.savefig(savefig_fn_posting)
                if show==0:
                    pl.close()
                else:
                    pl.show()
                return az, fs, D_pwv

            if os.path.exists(self.path_to_all+'am_datfiles_Az/'+template+'/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_import_model_pickle_temps.txt'):
                pickle_fn_temps=self.path_to_all+'am_datfiles_Az/'+template+'/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_import_model_pickle_temps.txt'
            else:
                pickle_fn_temps=self.path_to_all+'am_datfiles_Az/'+template+'/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_fit_pickle_temps.txt'

            if os.path.exists(self.path_to_all+'am_datfiles_Az/'+template+'/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'):
                pickle_fn=self.path_to_all+'am_datfiles_Az/'+template+'/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'
            else:
                pickle_fn=self.path_to_all+'am_datfiles_Az/'+template+'/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_fit_fitoutput.txt'

            print(fn)
            D_clean, waz_1, single_mod_removed_data, p2, p2_err, p1, p1_err, calib_data_Gavg, FH = self.read_Az_fast(fn)
            print('FH=', FH)
            print(np.shape(FH))
            tot=FH[:,0] #time ordered time

            #using fit or import_model should be the same here as I am just using it to extract waz
            f = open(pickle_fn_temps,'rb')
            pickle_Temps = pk.load(f)
            f.close()
            waz=pickle_Temps['wAz']

            print('pickle temps = ', pickle_Temps.keys())

            f = open(pickle_fn,'rb')
            fit_output = pk.load(f)
            f.close()
            pwv_ts=fit_output['pwv']

            print('fit output = ', fit_output.keys())

            az, fs, D_pwv = plot_atmo(waz, pwv_ts)

            calib_az=self.calib_pointing(az)
            #calib_az=az #not calibrating for now

            idx=np.where(waz<=360)

            return waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot


    def calib_pointing(self, wvr_az, point_par_fn='pointing_parameters_2018_fast.txt'):

        path = self.path_to_all + '../txtfiles/'

        f = open(path+point_par_fn,'rb')
        point_par = pk.load(f)
        f.close()

        real_az = wvr_az - point_par['az_offs']

        return real_az


    def extract_PowerSpectrum(self, time_ordered_time, dt, time_ordered_az, daz, time_ordered_data, data='WVR'):

        class struct(object):
            pass

        ps=struct();
        ps.t=struct();
        ps.az=struct();


        def ps_t(time_ordered_time, time_ordered_data):

            if data=='WVR':
                dtime=dt
            elif data=='BK':
                dtime=(time_ordered_time[1]-time_ordered_time[0]).total_seconds()

            #time_ordered_time=time_ordered_time[np.where(time_ordered_time != 0)]
            #time_ordered_data=time_ordered_data[np.where(time_ordered_time != 0)]


            #Power Spectrum in time domain
            time_ = time_ordered_time
            #T, dt = len(time_), np.median(np.gradient(time_))
            time = len(time_)
            df=1./dt
            # fmids = np.geomspace(1/(T*dt),.5/dt,256)
            # df    = np.exp(np.gradient(np.log(fmids)).mean())
            # fbins = np.append(fmids[0]/np.sqrt(df),fmids*np.sqrt(df))

            ps    = np.square(np.abs(np.fft.fft(time_ordered_data)/time))#*np.hanning(len(time_ordered_data)))))
            freq_ = np.fft.fftfreq(time,dtime)

            idx_pos=np.where(freq_>=0)

            return freq_[idx_pos], ps[idx_pos]



        def ps_az(az_ordered_az, time_ordered_az, az_ordered_data):

            #daz = time_ordered_az[1]-time_ordered_az[0]

            #az_ordered_az=az_ordered_az[np.where(az_ordered_az != 0)]
            #az_ordered_data=az_ordered_data[np.where(az_ordered_az != 0)]

            #Power Spectrum in time domain
            az_ = az_ordered_az

            #T, dt = len(time_), np.median(np.gradient(time_))
            az = len(az_)
            df=1./daz
            # fmids = np.geomspace(1/(T*dt),.5/dt,256)
            # df    = np.exp(np.gradient(np.log(fmids)).mean())
            # fbins = np.append(fmids[0]/np.sqrt(df),fmids*np.sqrt(df))
            ps    = np.square(np.abs(np.fft.fft(az_ordered_data)/az))#*np.hanning(len(az_ordered_data)))))
            freq_ = np.fft.fftfreq(az,daz)

            idx_pos=np.where(freq_>=0)

            return freq_[idx_pos], ps[idx_pos]



        def ps_one_channel(time_ordered_time, time_ordered_az, time_ordered_data):

            #time domain ps
            ps.t.freq, ps.t.ps = ps_t(time_ordered_time, time_ordered_data)

            az_ordered_az=np.sort(time_ordered_az)
            idx=np.argsort(time_ordered_az)
            az_ordered_data=time_ordered_data[idx]

            #az_ordered_az=az_ordered_az[np.where(az_ordered_az != 0)]
            #az_ordered_data=az_ordered_data[np.where(az_ordered_az != 0)]

            #space domain ps
            ps.az.freq, ps.az.ps= ps_az(az_ordered_az, time_ordered_az, az_ordered_data)

            return ps


        ps_one_channel = ps_one_channel(time_ordered_time, time_ordered_az, time_ordered_data)


        return ps_one_channel


    #
    #
    #
    # def extract_PowerSpectrum(self, D_pwv, D_bk_sum, D_bk_diff, az_pwv, x_az, ):
    #
    #     class struct(object):
    #         pass
    #
    #     #fig, ax = pl.subplots(7,1)
    #     var_i=0
    #
    #     D_dict={'pwv': [], 'T_bk_sum':[], 'T_bk_diff':[]}
    #     ps_norm_dict={'pwv': [], 'T_bk_sum':[], 'T_bk_diff':[]}
    #
    #     D_dict['pwv'] = D_pwv
    #     D_dict['T_bk_sum'] = D_bk_sum
    #     D_dict['T_bk_diff'] = D_bk_diff
    #
    #     if x_az[int(len(x_az)/2)] <= 0:
    #         az_bk = [(360.+i) for i in x_az]
    #     else:
    #         az_bk = x_az
    #
    #     for key in D_dict:
    #
    #         D_var=D_dict[key]
    #         ps_list=[]
    #         ang_freqs_list=[]
    #
    #         if key == 'pwv':
    #             az = az_pwv
    #         else:
    #             az = az_bk
    #
    #         ps_az=np.zeros(len(D_var[:,0]))
    #         for i in range(len(D_var[:,0])):
    #             var_az=D_var[i,:]
    #             ps_az[i]=np.nanmean(var_az)
    #
    #         ps_az=ps_az-np.nanmean(ps_az)
    #
    #         nan_array = np.isnan(ps_az)
    #         not_nan_array = ~ nan_array
    #         ps_az = ps_az[not_nan_array]
    #
    #         pl.plot(ps_az)
    #         pl.title('ps_az -'+str(key))
    #         pl.show()
    #
    #         FT=np.fft.fft(ps_az)
    #         ps = np.abs(FT)**2
    #         #ang_step = 1/2.
    #         ang_step=np.abs(az[1]-az[0])
    #         #Norm=np.sum(ps)*ang_step
    #         Norm=1 #not normalized
    #         ps_norm=ps/Norm
    #         ang_freqs = np.fft.fftfreq(ps_az.size, ang_step)
    #         print('key=', key)
    #         print('ang_step=', ang_step)
    #         print('ps_az.size=', ps_az.size)
    #         print('az=', az[0:10])
    #         idx = np.argsort(ang_freqs[np.where(ang_freqs>0)])
    #
    #         ang_freqs=ang_freqs[idx]
    #         ps_norm=ps_norm[idx]
    #         #ps_norm_smooth = butter_lowpass_filtfilt(ps_norm, normal_cutoff)
    #         #print('var_i=', var_i)
    #         #print('ps_norm=', ps_norm)
    #         data_struct = struct()
    #         data_struct.freqs = ang_freqs
    #         data_struct.ps=ps_norm
    #         ps_norm_dict[key]=data_struct
    #
    #         var_i+=1
    #
    #     fig, (ax, ax2) = plt.subplots(2,1)
    #     for key in ps_norm_dict:
    #         if not (key == 'pwv'):
    #             data = ps_norm_dict[key]
    #             ax.scatter(data.freqs, data.ps, s=3, label=key)
    #             ax.plot(data.freqs, data.ps, c='k', alpha=0.5)
    #             ax.set_xscale('log')
    #             ax.set_yscale('log')
    #             ax.set_ylabel('T_bk\n[K^2]')
    #
    #     pwv_data=ps_norm_dict['pwv']
    #     #ax2=ax.twinx()
    #     ax2.scatter(pwv_data.freqs, pwv_data.ps, s=3, c='r', label='pwv')
    #     ax2.plot(pwv_data.freqs, pwv_data.ps, c='k', alpha=0.5)
    #     ax2.set_ylabel('PWV\n[um^2]')
    #     #ax[var_i].set_ylim(0,0.0012)
    #     #ax[var_i].set_xlim(0.05, np.max(ang_freqs))
    #     ax2.set_xscale('log')
    #     ax2.set_yscale('log')
    #
    #     pl.legend(loc='lower left')
    #     pl.suptitle('Power Spectrum')
    #     pl.xlabel('theta^(-1)[deg^(-1)]')
    #     pl.show()
