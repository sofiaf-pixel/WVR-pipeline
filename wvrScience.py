import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import matplotlib.cm as cmap
from matplotlib.dates import DateFormatter, DateLocator, AutoDateLocator, num2date, date2num
from operator import itemgetter
from itertools import groupby
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import PercentFormatter
import wvrReadData

import analysisUtils as au
import wvrReadData as wrd
from initialize import initialize
import reduc_wvr_pager as rw
import numpy as np

import pickle

#import scipy.io
# data = {
#    'bigdata' : {
#        'a' : array([1, 2, 3]),
#        'b' : array([1, 2, 3]),
#        'c' : array([1, 2, 3]),
#     }
#}
#scipy.io.savemat('test.mat', data)

def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
        
class wvrScience(initialize):
    '''

    '''
    def __init__(self, unit=None, verb=True):

        '''

        
        '''
        initialize.__init__(self, unit, verb=verb)
        self.wvrR = wrd.wvrReadData(self.unit, verb=verb)



        
    def plotAtmogram(self, fileList, td0, tw0, inter=False, verb=True, fitphase = 0):
        '''
        Created by NL 20161010
        Re-written by dB 20161110
        TODO: Store data in pickle as intermediate product

        '''
        
        if inter:
            ion()
        else:
            ioff()


        spath='/n/home02/sfatigoni/work/plots/InputLossCorrection/'
            
        fileList = au.aggregate(fileList)
        nfiles = size(fileList)                 # nfiles to analyze
        pcoef= {0:None, 1:None, 2:None, 3:None} # init fit coef for 4 chans
        D = {0:None, 1:None, 2:None, 3:None}    # init maps for 4 chans
        Dres = {0:None, 1:None, 2:None, 3:None} # init map residuals for 4 chans
        Dbl = {0:None, 1:None, 2:None, 3:None}  # init maps baseline for 4 chans

        pcoef_double= {0:None, 1:None, 2:None, 3:None} # init fit coef for 4 chans                           
        D_double = {0:None, 1:None, 2:None, 3:None}    # init maps for 4 chans                              
        Dres_double = {0:None, 1:None, 2:None, 3:None} # init map residuals for 4 chans
        Dbl_double = {0:None, 1:None, 2:None, 3:None}  # init maps baseline for 4 chans    

        Dtau = {0:None, 1:None, 2:None, 3:None}
        Dw = {0:None, 1:None, 2:None, 3:None} 
        Dw_c = {0:None, 1:None, 2:None, 3:None}
        
        c = ['b','r','g','m']

        if verb: print "Loading %d slow files:"%nfiles
        utslow,tslow,d,azslow,elslow,tsrc = self.wvrR.readSlowFile(fileList)
        nchan = shape(tsrc)[1]
        dres = zeros(shape(tsrc))  # init residuals on the 4 tsrc
        dres_double = zeros(shape(tsrc))

        if size(tslow) == 1: return
        waz, fs = self.findScans(azslow)
        fname = fileList[0].split('_')
        
        if nfiles > 1:
            fileslow = '%s_2400.txt'%fname[0]
            figsize=(36,12)
            trange=[utslow[0].replace(hour=0,minute=0,second=0), utslow[-1].replace(hour=23,minute=59,second=59)]
        else: 
            fileslow = '%s_%s.txt'%(fname[0],fname[1][0:4])
            figsize=(12,10)
            trange=[utslow[0].replace(minute=0,second=0),utslow[-1].replace(minute=59,second=59)]

        # majorloc = AutoDateLocator(minticks=5, maxticks=12, interval_multiples=True) 
        # df = DateFormatter('%H:%M')
        phase_ch=[109.18115635, 107.62194176, 103.16825431,  94.22004318]
        phase_double_ch=[188.32 , 197.49 , 200.075, 202.66]

        #Loop through channels
        nscans = len(fs.s)
        amp=np.zeros(nscans)
        phase_data = zeros(4)
        phase_data_double = zeros(4)
        new_data=np.zeros(shape(tsrc))


        #td0_list=[0.021,0.025,0.027,0.029,0.033]                                                                                              
        #tw0_list=[1.04,1.14,1.19,1.25,1.34] 
        
        td=[0.027, 0.020, 0.010, 0.010]#0.027
        tw=[1.19, 0.40, 0.18, 0.12]#1.19

        #td[0]=td0
        #tw[0]=tw0

        print('td=', td)
        print('tw=', tw)
        
        for i in range(nchan):
                
            phase=phase_ch[i]
            phase_double=phase_double_ch[i]

            #Double mod cleaning
            print("FitDouble_1")
            res_double, pcoef0_double, baseline_double, T_double = self.filterScans(waz, tsrc[:,i], fs, 'sin', phase_double_ch[i], amp, fitphase=0, fitamp=0, k=2)
            dres_double[:,i]=res_double
            offset=pcoef0_double[:,0,0]
            amp_double,p=self.FindAmpDoubleMod(season='2018',offs=offset)
            print("FitDouble_2")
            res_double, pcoef0_double, baseline_double, T_double = self.filterScans(waz, res_double, fs, 'sin', phase_double_ch[i], amp_double, fitphase=0, fitamp=0, k=2)
            
            #new_data[:,i]=res-ffourier_double #just Sky + Single Mod
            pcoef_double[i] = pcoef0_double
            D_double[i] = self.interpToImage(waz, tsrc[:,i], fs)
            Dres_double[i] = self.interpToImage(waz, res_double, fs)
            Dbl_double[i] = self.interpToImage(waz, baseline_double, fs)
            sd = shape(D[0])
                        
            #Single mod cleaning

            amp_single,p=self.FindAmpSingleMod(season='2018',offs=offset)
            print("FitSingle_1")
            res, pcoef0, baseline, Tb = self.filterScans(waz, res_double, fs, 'sin', phase_ch[i], amp_single, fitphase=0, fitamp=0, k=1)

            tau, w, tau_c, w_c = self.TemperaturetoPWV_inputtw(waz, Tb, 258.7, 0.02, 290, td[i], tw[i])
            #tau_c,w_c=self.TemperaturetoPWV(i, waz, Tb, Tatm=258.7, L=0.02, Tloss=290)
            #tau,w=self.TemperaturetoPWV(i, waz, Tb, Tatm=258.7, L=0.0, Tloss=290)
            
            
            dres[:,i]=res
            pcoef[i] = pcoef0
            D[i] = self.interpToImage(waz, T_double, fs)
            Dres[i] = self.interpToImage(waz, Tb, fs)
            Dbl[i] = self.interpToImage(waz, baseline, fs)
            sd = shape(D[0])

            Dtau[i]=self.interpToImage(waz, tau, fs)
            Dw[i]=self.interpToImage(waz, w, fs)
            Dw_c[i]=self.interpToImage(waz, w_c, fs)

            print('Dw_c', Dw_c[i])


            
        #plt scatter trial plot
        plt.figure()
        plt.scatter(Dw_c[0], Dw_c[1], s=0.5, label='Channel 1')
        plt.scatter(Dw_c[0], Dw_c[2], s=0.5, label='Channel 2')
        plt.scatter(Dw_c[0], Dw_c[3], s=0.5, label='Channel 3')
        plt.legend()
        plt.xlabel('Channel 0')
        plt.suptitle('PWV Scatter Plot')
        plt.title('Corrected L=0.02 tw/td nominal')
        plt.savefig(spath+'L=0p02_correction_'+str(fname[0])+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
        plt.close()
        

        #plot scatter plot for all combinations of tw0 and td0
        #self.tdtw_correction(fs,waz,Tb,Tatm=258.7,L=0.02,Tloss=290)

                    

        
        figsize=(12,10)
        figure(10, figsize=figsize);clf()
        for i in range(0,4):
                sp1 = subplot2grid((7,2), (4*(i/2)+0, mod(i,2)), colspan=1)
                imshow(D_double[i],aspect='auto',interpolation='nearest', origin='lower')
                sp1.set_xticks(range(0,sd[1],10))
                sp1.set_xticklabels('')
                sp1.set_yticks(range(0,sd[0],60))
                sp1.set_title('TSRC%s'%i)
                plt.colorbar()

                sp2 = subplot2grid((7,2), (4*(i/2)+1, mod(i,2)), colspan=1)
                imshow(Dres_double[i],aspect='auto',interpolation='nearest', origin='lower')
                sp2.set_xticks(range(0,sd[1],10))
                sp2.set_xticklabels('')
                sp2.set_yticks(range(0,sd[0],60))
                plt.colorbar()
                
                sp3 = subplot2grid((7,2), (4*(i/2)+2, mod(i,2)), colspan=1)
                imshow(Dbl_double[i],aspect='auto',interpolation='nearest', origin='lower')
                sp3.set_xticks(range(0,sd[1],10))
                sp3.set_yticks(range(0,sd[0],60))
                sp3.set_ylabel('Az( [deg]')
                sp3.set_xlabel('scan number')
                sp3.grid(False)
                plt.colorbar()

                if i==3:
                        subplots_adjust(hspace=0.01)
                        title = fileslow.replace('.txt','_atmogram')
                        suptitle(title,y=0.97, fontsize=20)
                        if verb: print "Saving %s.png"%title
                        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
                        plt.show()

            
        figure(11, figsize=figsize);clf()
        for i in range(0,4):
                
                sp4 = subplot2grid((7,2), (4*(i/2)+0, mod(i,2)), colspan=1)
                imshow(D[i],aspect='auto',interpolation='nearest', origin='lower')
                sp4.set_xticks(range(0,sd[1],10))
                sp4.set_xticklabels('')
                sp4.set_yticks(range(0,sd[0],60))
                sp4.set_title('TSRC%s'%i)

                sp5 = subplot2grid((7,2), (4*(i/2)+1, mod(i,2)), colspan=1)
                imshow(Dres[i],aspect='auto',interpolation='nearest', origin='lower')
                sp5.set_xticks(range(0,sd[1],10))
                sp5.set_xticklabels('')
                sp5.set_yticks(range(0,sd[0],60))

                sp6 = subplot2grid((7,2), (4*(i/2)+2, mod(i,2)), colspan=1)
                imshow(Dbl[i],aspect='auto',interpolation='nearest', origin='lower')
                sp5.set_xticks(range(0,sd[1],10))
                sp5.set_yticks(range(0,sd[0],60))
                sp5.set_ylabel('Az( [deg]')
                sp5.set_xlabel('scan number')
                sp5.grid(False)

                if i==3:
                        subplots_adjust(hspace=0.01)
                        title = fileslow.replace('.txt','_atmogram')
                        suptitle(title,y=0.97, fontsize=20)
                        if verb: print "Saving %s.png"%title
                        savefig(title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
                        plt.show()


                figure(11, figsize=figsize);clf()
        for i in range(0,4):

                sp4 = subplot2grid((7,2), (4*(i/2)+0, mod(i,2)), colspan=1)
                imshow(D[i],aspect='auto',interpolation='nearest', origin='lower')
                sp4.set_xticks(range(0,sd[1],10))
                sp4.set_xticklabels('')
                sp4.set_yticks(range(0,sd[0],60))
                sp4.set_title('TSRC%s'%i)
                plt.colorbar()

                sp5 = subplot2grid((7,2), (4*(i/2)+1, mod(i,2)), colspan=1)
                imshow(Dres[i],aspect='auto',interpolation='nearest', origin='lower')
                sp5.set_xticks(range(0,sd[1],10))
                sp5.set_xticklabels('')
                sp5.set_yticks(range(0,sd[0],60))
                plt.colorbar()
                
                sp6 = subplot2grid((7,2), (4*(i/2)+2, mod(i,2)), colspan=1)
                imshow(Dbl[i],aspect='auto',interpolation='nearest', origin='lower')
                sp5.set_xticks(range(0,sd[1],10))
                sp5.set_yticks(range(0,sd[0],60))
                sp5.set_ylabel('Az( [deg]')
                sp5.set_xlabel('scan number')
                sp5.grid(False)
                plt.colorbar()

                if i==3:
                        subplots_adjust(hspace=0.01)
                        title = fileslow.replace('.txt','_atmogram')
                        suptitle(title,y=0.97, fontsize=20)
                        if verb: print "Saving %s.png"%title
                        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
                        plt.show()


                    
        #plot the tau/w

        figure(14, figsize=figsize);clf()
        for i in range(0,4):

                sp4 = subplot2grid((7,2), (4*(i/2)+0, mod(i,2)), colspan=1)
		imshow(Dres[i],aspect='auto',interpolation='nearest', origin='lower')
		sp4.set_xticks(range(0,sd[1],10))
                sp4.set_xticklabels('')
                sp4.set_yticks(range(0,sd[0],60))
                sp4.set_title('TSRC%s'%i)
                plt.colorbar()
                
                sp5 = subplot2grid((7,2), (4*(i/2)+1, mod(i,2)), colspan=1)
                imshow(Dtau[i],aspect='auto',interpolation='nearest', origin='lower')
                sp5.set_xticks(range(0,sd[1],10))
                sp5.set_xticklabels('')
                sp5.set_yticks(range(0,sd[0],60))
                plt.colorbar()
                
                sp6 = subplot2grid((7,2), (4*(i/2)+2, mod(i,2)), colspan=1)
                imshow(Dw[i],aspect='auto',interpolation='nearest', origin='lower')
                sp5.set_xticks(range(0,sd[1],10))
                sp5.set_yticks(range(0,sd[0],60))
                sp5.set_ylabel('Az( [deg]')
                sp5.set_xlabel('scan number')
                sp5.grid(False)
                plt.colorbar()

                if i==3:
                        subplots_adjust(hspace=0.01)
                        title = fileslow.replace('.txt','_PWVatmogram')
                        suptitle(title,y=0.97, fontsize=20)
                        if verb: print "Saving %s.png"%title
                        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
                        plt.show()

                        

        figure(13, figsize=figsize);clf()
        ax1 = plt.subplot(411)
        plt.imshow(Dres[0],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        plt.ylabel('Az [deg]')
        plt.title('TSRC0')
        plt.colorbar()

        ax2 = plt.subplot(412, sharex=ax1)
        plt.imshow(Dres[1],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel('Az [deg]')
        plt.title('TSRC1')
        plt.colorbar()
        
        ax3 = plt.subplot(413, sharex=ax1)
        plt.imshow(Dres[2],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.title('TSRC2')
        plt.colorbar()
        
        ax4 = plt.subplot(414, sharex=ax1)
	plt.imshow(Dres[3],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.xlabel('Scan_n')
        plt.title('TSRC3')
        plt.colorbar()


        title = fileslow.replace('.txt','_TSRC')
        suptitle(title,y=0.97, fontsize=20)
        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
        
        plt.show()


        figure(14, figsize=figsize);clf()
        ax1 = plt.subplot(411)
        plt.imshow(Dw_c[0],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        plt.ylabel('Az [deg]')
        plt.title('PWV0')
        plt.suptitle('PWV after Correction for Input Losses')
        plt.colorbar()

        ax2 = plt.subplot(412, sharex=ax1)
        plt.imshow(Dw_c[1],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel('Az [deg]')
        plt.title('PWV1')
	plt.colorbar()

        ax3 = plt.subplot(413, sharex=ax1)
        plt.imshow(Dw_c[2],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.title('PWV2')
	plt.colorbar()

        ax4 = plt.subplot(414, sharex=ax1)
        plt.imshow(Dw_c[3],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.xlabel('Scan_n')
        plt.title('PWV3')
        plt.colorbar()

        title = fileslow.replace('.txt','_PWV_corrected')
	suptitle(title,y=0.97, fontsize=20)
        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
        
        
        plt.show()



        figure(15, figsize=figsize);clf()
        ax1 = plt.subplot(411)
        plt.imshow(Dw[0],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        plt.ylabel('Az [deg]')
        plt.title('PWV0')
        plt.suptitle('PWV No Correction for Input Losses')
        plt.colorbar()

        ax2 = plt.subplot(412, sharex=ax1)
        plt.imshow(Dw[1],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel('Az [deg]')
        plt.title('PWV1')
        plt.colorbar()

        ax3 = plt.subplot(413, sharex=ax1)
        plt.imshow(Dw[2],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.title('PWV2')
        plt.colorbar()

        ax4 = plt.subplot(414, sharex=ax1)
        plt.imshow(Dw[3],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.xlabel('Scan_n')
        plt.title('PWV3')
        plt.colorbar()

        title = fileslow.replace('.txt','_PWV')
        suptitle(title,y=0.97, fontsize=20)
        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
        
        plt.show()


        figure(16, figsize=figsize);clf()
        ax1 = plt.subplot(411)
        plt.imshow(Dw_c[0],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax1.get_xticklabels(), fontsize=6)
        plt.ylabel('Az [deg]')
        plt.title('PWV0')
        plt.colorbar()

        ax2 = plt.subplot(412, sharex=ax1)
        plt.imshow(Dw_c[1],aspect='auto',interpolation='nearest', origin='lower')
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.ylabel('Az [deg]')
        plt.title('PWV1')
        plt.colorbar()

        ax3 = plt.subplot(413, sharex=ax1)
        plt.imshow(Dw_c[2],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
        plt.title('PWV2')
        plt.colorbar()

        ax4 = plt.subplot(414, sharex=ax1)
        plt.imshow(Dw_c[3],aspect='auto',interpolation='nearest', origin='lower')
        plt.ylabel('Az [deg]')
	plt.xlabel('Scan_n')
        plt.title('PWV3')
        plt.colorbar()

        title = fileslow.replace('.txt','_PWV_corrected')
	suptitle(title,y=0.97, fontsize=20)
        savefig(spath+title+'_td0_'+str(td[0])+'_tw0_'+str(tw[0])+'.png')
        
        plt.show()
        

        
        if not inter:
            close('all')
        au.movePlotsToReducDir(self.reducDir)

       # for n in range (0, len(fileList)-1):
        amp=np.zeros((4,len(pcoef[0][:,0,0])))
        offs=np.zeros((4,len(pcoef[0][:,1,0])))

        for i in range (0,4):
                amp[i]=pcoef[i][:,0,0]
                offs[i]=pcoef[i][:,1,0]

        if fitphase==0:
                a = open('Amplitude_fixed.txt','w')
                p = open('Phase_fixed.txt','w')
                o = open('Offs_fixed.txt','w')

                pickle.dump(amp, a)
                pickle.dump(phase_ch, p)
                pickle.dump(offs, o)

        if fitphase==1:
                a = open('Amplitude.txt','w')
                p = open('Phase.txt','w')
                o = open('Offs.txt','w')

	        pickle.dump(amp, a)
                pickle.dump(phase_data, p)
                pickle.dump(offs, o)
         
        a.close()
        p.close()
        o.close()

        return waz, D, Dres, Dbl, pcoef, D_double, Dres_double, Dbl_double, Dw, Dw_c, pcoef_double, phase_data


    def TemperaturetoPWV(self, i, waz, res, Tatm, L, Tloss):

            td=[0.027, 0.020, 0.010, 0.010]
            tw=[1.19, 0.40, 0.18, 0.12]
            
            tau=np.log(Tatm/(Tatm-res))
            w=(tau-td[i])/tw[i]

            #corrected
            #Tloss=290
            #L=0.03
            Tc=res-(L*Tloss)/(1-L)
            tau_c=np.log(Tatm/(Tatm-Tc))
            w_c=(tau_c-td[i])/tw[i]
            
            return tau_c, w_c

    def TemperaturetoPWV_inputtw(self, waz, res, Tatm, L, Tloss, td, tw):

            tau=np.log(Tatm/(Tatm-res))
            w=(tau-td)/tw

            #corrected                                                                                                               
            Tc=res-(L*Tloss)/(1-L)
            tau_c=np.log(Tatm/(Tatm-Tc))
            w_c=(tau_c-td)/tw

            return tau, w, tau_c, w_c

    def tdtw_correction(self, fs, waz, res, Tatm, L, Tloss):

        spath='/n/home02/sfatigoni/work/plots/InputLossCorrection/'
            
        Dw_v = {0:None, 1:None, 2:None, 3:None}
        Dw_c_v = {0:None, 1:None, 2:None, 3:None}
            
        nchan = 4
	
        td0_list=[0.021,0.025,0.027,0.029,0.033]
        tw0_list=[1.04,1.14,1.19,1.25,1.34]

        for td0 in td0_list:
                for tw0 in tw0_list:
                        td=[td0, 0.020, 0.010, 0.010]
                        tw=[tw0, 0.40, 0.18, 0.12]

                        for i in range(nchan):
                                tau_v, w_v, tau_c_v, w_c_v=self.TemperaturetoPWV_inputtw(waz, res, Tatm, L, Tloss, td[i], tw[i])
                                Dw_c_v[i]=self.interpToImage(waz, w_c_v, fs)

                        plt.figure()
                        plt.scatter(Dw_c_v[0], Dw_c_v[1], s=0.5, label='Channel 1')
                        plt.scatter(Dw_c_v[0], Dw_c_v[2], s=0.5, label='Channel 2')
                        plt.scatter(Dw_c_v[0], Dw_c_v[3], s=0.5, label='Channel 3')
                        plt.legend()
                        plt.xlabel('Channel 0')
                        plt.suptitle('PWV Scatter Plot')
                        plt.title('Corrected L=0.02 - td0='+str(td0)+'-tw0='+str(tw0))
                        plt.savefig(spath+'L=0p02_correction_td0='+str(td0)+'-tw0='+str(tw0)+'.png')
                        plt.close()

    

    def pltAtm(self, D, Dres, Dbl, nchan, fileslow, verb):
            sd = shape(D[0])
            figsize=(10,12)
            c = ['b','r','g','m']

            for i in range (0,3):
            
                    sp4 = subplot2grid((7,2), (4*(i/2)+0, mod(i,2)), colspan=1)
                    imshow(D[i],aspect='auto',interpolation='nearest', origin='lower')
                    sp4.set_xticks(range(0,sd[1],10))
                    sp4.set_xticklabels('')
                    sp4.set_yticks(range(0,sd[0],60))
                    sp4.set_title('TSRC%s'%i)

                    sp5 = subplot2grid((7,2), (4*(i/2)+1, mod(i,2)), colspan=1)
                    imshow(Dres[i],aspect='auto',interpolation='nearest', origin='lower')
                    sp5.set_xticks(range(0,sd[1],10))
                    sp5.set_xticklabels('')
                    sp5.set_yticks(range(0,sd[0],60))

                    sp6 = subplot2grid((7,2), (4*(i/2)+2, mod(i,2)), colspan=1)
                    imshow(Dbl[i],aspect='auto',interpolation='nearest', origin='lower')
                    sp6.set_xticks(range(0,sd[1],10))
                    sp6.set_yticks(range(0,sd[0],60))
                    sp6.set_ylabel('Az( [deg]')
                    sp6.set_xlabel('scan number')
                    sp6.grid(False)

                    if i==3:     
                            subplots_adjust(hspace=0.01)
                            title = fileslow.replace('.txt','_atmogram')
                            suptitle(title,y=0.97, fontsize=20)
                            if verb: print "Saving %s.png"%title
                            savefig(title+'.png')
                            plt.show()

                            
                            
    def pltParameters(self, i, pcoef, nchan, sd, fitphase, phase_ch, fileslow):
            figsize=(10,12)
            c = ['b','r','g','m']
            
            figure(12, figsize=figsize);clf()
            sp1 = subplot(3,2,1)
            plot(pcoef[i][:,0,0],'o-',color=c[i])
            ylabel('Sin Amp [K]')
            sp1.set_xticklabels('')
            xl=sp1.set_xlim([-2,sd[1]])
            if i == 0: sp1.grid(which='both')
            if i == 0:
                    yl= sp1.set_ylim([-0.1,1.])
                    text(xl[0],yl[1],'Fit Coeffs')

            sp2 = subplot(3,2,3)
            if fitphase==1:
                    phase_plot=np.full(len(pcoef[i][:,0,0]), phase_data[i])
                    plot(phase_plot,'o-',color=c[i])
            else:
                    phase_plot=np.full(len(pcoef[i][:,0,0]), phase_ch[i])
                    plot(phase_plot,'o-',color=c[i])
            ylabel('Sin Phase [deg]')
            sp2.set_xlim([-2,sd[1]])
            sp2.set_xticklabels('')
            yl= sp2.set_ylim([180,220])
            if i == 0: sp2.grid(which='both')

            sp3 = subplot(3,2,5)
            plot(pcoef[i][:,1,0],'o-',color=c[i])
            ylabel('Sin Offset [K]')
            yl= sp3.set_ylim([0,250])
            sp3.set_xlim([-2,sd[1]])
            sp3.set_xlabel('scan number')
            if i == 0: sp3.grid(which='both')

            sp4 = subplot(3,2,2)
            plot(pcoef[i][:,0,1],'o',color=c[i])
            ylabel('Sin Amplitude err [K]')
            sp4.set_xticklabels('')
            xl=sp4.set_xlim([-2,sd[1]])
            if i == 0: sp4.grid(which='both')
            if i== 0:
                yl= sp4.set_ylim([-0.1,.2])
                text(xl[0],yl[1],'Fit Coeffs Errors')

            sp5 = subplot(3,2,4)
            plot(np.zeros(len(phase_plot)),'o-',color=c[i])
            ylabel('Sin Phase err [deg]')
            sp5.set_xticklabels('')
            sp5.set_xlim([-2,sd[1]])
            yl= sp5.set_ylim([-0.05,0.05])
            if i == 0: sp5 .grid(which='both')

            sp6 = subplot(3,2,6)
            plot(pcoef[i][:,1,1],'o',color=c[i])
            if i == 0: sp6.grid(which='both')
            ylabel('Sin Offset err[K]')
            sp6.set_xlabel('scan number')
            sp6.set_xlim([-2,sd[1]])
            yl= sp6.set_ylim([0,0.15])

            if i==3:
                sp1.legend(['Tsrc0','Tsrc1','Tsrc2','Tsrc3'],loc=1,prop={'size':8})
                subplots_adjust(hspace=0.01)
                title = fileslow.replace('.txt','_sinFits')
                suptitle(title,y=0.97, fontsize=20)
                if verb: print "Saving %s.png"%title
                savefig(title+'.png')
                plt.show()

    def pltRes(self, i, Dres, nchan, sd):
            figsize=(10,12)
            c = ['b','r','g','m']
            
            az_temp={}
            figure(13, figsize=figsize);clf()
            for i in range(nchan):
                sp = subplot(5,1,i+1)
            for j in range(sd[1]):
                plot(Dres[i][:,j],'.', color=c[i])
            az_temp[i] = nanmean(Dres[i],axis=1)
            plot(az_temp[i],'k-')
            sp.grid(which='both')
            ylabel('Tsrc%s [K]'%i)
            sp.set_xticklabels('')
            ylim([-1,1])
            sp.set_xlim([-2,sd[0]])

            sp2 = subplot(5,1,5)
            plot(az_temp[i], '-',color=c[i])
            sp2.set_xlim([-2,sd[0]])
            ylim([-1,1])
            xlabel('Az [deg]')
            ylabel('Tsrc [K]')

            if i==3:
                subplots_adjust(hspace=0.01)
                title = fileslow.replace('.txt','_residuals')
                suptitle(title,y=0.97, fontsize=20)
                if verb: print "Saving %s.png"%title
                savefig(title+'.png')
                plt.show()




    def PlotAtmogrambyDate(self, unit, start, end, inter=False, verb=True, fitphase = True):

        rwp = rw.reduc_wvr_pager(unit)
        fl = rwp.makeFileListFromData(start=start, end=end, typ='scanAz')

        waz, D, Dres, Dbl, pcoef = self.plotAtmogram(fl,inter=inter,verb=verb,fitphase=fitphase)

        return waz, D, Dres, Dbl, pcoef



        ##############################
    #### Helper functions and fitting functions
    
    def makeTiltModel(self, az):
        '''
        function to write out a tilt model which saves for each tag (scanAz tags mostly)
          tag, datetime, theth_az, phi_az 
        '''

        # get a file list of azscan
        # find the  amp ,ang and offset for each of those 1 hour tags. using testAngFit
        # find the relationship between the building tilts and the amp and ang found from the sin fits. Do this using  findAngles which does a minimization
        # then write a new function ( this one) which saves in a structure the  tiltModel along with time. The goal is to have for each 1hr tag a phase/amp to be used for the  actual data  fit.

        class stru:
            def __init__(self):
                self.num = []
                self.s = []
                self.e = []

    def findAngles(self,tilts, p, channel = 0):

        """
        utslow, tilts = wvrR.readTiltFile(fl_tilt)
        p = wvrS.testAngFit(fl)
        
        Fit theta_az = arctan(ytilt+C1/xtilt+C1) +C2
            phi_az = sqrt(xtilt^2+ytilt^2) +C3
            C3 related C1/C2
            Loop over C1 and C2
        """
        ch = channel
        el0 = 55
        
        # define the  theta_az function of tilts
        def theta_az(tilts,x0):
            C1,C2,C3 = x0
            return rad2deg(arctan2(tilts[:,2]+C1,tilts[:,0]+C2))+C3

        #define residual
        def theta_az_res(x0,tilts,ph):
            theta = theta_az(tilts,x0)
            theta_wvr = ph
            residual = theta - theta_wvr
            return residual

        # take median of 4 channels phase.
        ph = median(p[:,:,1,0],0)  
        sol = least_squares(theta_az_res, [1.,1,-80],args=(tilts,ph))
        theta_az_sol = theta_az(tilts,sol.x)
        print sol.success, sol.x

        bins = arange(-100,100,1)
        figure(2);clf()
        subplot(2,1,1)
        plot(ph,'.b')
        plot(theta_az_sol,'g.')
        xlabel('time since start of season [tags]')
        ylabel('theta az [deg]')
        legend(['per tag fitted phase of WVR data','measured building pitch/roll adjusted to fit phase'])
        subplot(2,1,2)
        hist(theta_az_sol-ph,bins=bins)
        xlabel('histogram of diff')
        #savefig('pitch_roll_fit_to_per_tag_phase.png')
        # 20161129: best fit is xf = arctan2(roll+0.945,pitch+0.945)-128

        # solve for Tatm, tau0, phi 
        # get factors D1,D2,D3
        def phi_az(x0,tilts,dc,mod):
            D1, D2, D3, D4 = x0
            phi = sqrt((tilts[:,0]+D1)**2+(tilts[:,2]+D2)**2)
            phi_wvr = rad2deg(D3*mod/(dc-D4))
            return phi, phi_wvr
        
        def phi_az_res(x0,tilts,dc,mod):
            phi,phi_wvr = phi_az(x0,tilts,dc,mod)
            return phi - phi_wvr
            
        # loop over channels
        for i in range(4):
            dc = p[i,:,2,0]
            mod = p[i,:,0,0]
            sol = least_squares(phi_az_res, [1.,1,1.0,-270],args=(tilts,dc,mod))
            phi_az_sol = phi_az(tilts,sol.x)
            print sol.success, sol.x

        print "Tatm=%3.2f, tau0=%3.2f"%(D2, tau0)
        phi_az = sqrt((tilts[:,0]+D1)**2+(tilts[:,2]+D1)**2)
        phi_az_wvr = rad2deg(D3*p[ch,:,0,0]/(p[ch,:,2,0]-D2))
        
        figure(4);clf()
        subplot(2,1,1)
        plot(phi_az,'g.')
        plot(phi_az_wvr,'b.')
        xlabel('time since start of season [tags]')
        ylabel('phi az [deg]')
        bins = arange(-.5,.5,0.01)
        subplot(2,1,2)
        hist(phi_az-phi_az_wvr,bins=bins)
        xlabel('histogram of diff')
        #savefig('pitch_roll_fit_to_per_tag_phase.png')
        return C1, C2, D1, D2, D3

    def testAngFit(self, fl):
        """
        For each 1hr of data, read the azScan dataset.
        Use the slow data.
        For each of the 4 channels
        Remove p0 (DC level) at each 360-scan
        Then fit a single sine wave to the whole 1-hr observation to obtain the phase. 
        We also get the amplitude and the DC offset but those are irrelevant.
        """
        #bounds=((0, -180, 0), (inf, 180, inf))
        
        # define the function to fit to:
        def sinusoidal(x,amplitude,phase,offset):
            return amplitude * sin(deg2rad(x+phase)) + offset
        
        nobs = size(fl)
        pcoef = zeros([4, nobs, 3, 2]) # 4chans x nobs x 3=amp,phase,offset x 2=val,error
        fit_err=zeros([4,len(fl),3])

        check_plot=0

        
        for i,f in enumerate(fl):
            utslow,tslow,d,azslow,elslow,tsrc = self.wvrR.readSlowFile(f)
            if size(d) == 1: continue
            waz,fs =  self.findScans(azslow)
            data = zeros([4,len(tslow)])

            
            for j in range(4):
                def sinusoidal_amp(x, amplitude, offset):
                        phase= 0.#load from file
                        return amplitude * sin(deg2rad(x+phase)) + offset

                res, pcoef0, baseline = self.filterScans(waz, tsrc[:,j], fs, 'p0')
                #fit, pcov= curve_fit(sinusoidal, azslow, smooth(res,5), p0=[1.0,-80.0,0.0], bounds=bounds)
                fit, pcov= curve_fit(sinusoidal_amp, azslow, smooth(res,5), p0=[1.0,0.0,0.0])
                fit_err[j,i,:] = sqrt(diag(pcov))

                pcoef[j,i,0,:]=[fit[0],fit_err[j,i,0]]  
                pcoef[j,i,1,:]=[fit[1],fit_err[j,i,1]]
                
                if pcoef[j,i,1,0]<= -180.:
                        pcoef[j,i,1,0]=180-(-pcoef[j,i,1,0]-180.)
                
                #pcoef[j,i,2,:]=[mean(pcoef0),1]
                pcoef[j,i,2,:]=[fit[2],fit_err[j,i,2]]
                data[j,:]=res

         
            if i==40 and check_plot==1:

                fig=plt.figure()

                ax1 = fig.add_subplot(411)
                ax1.scatter(tslow, data[0,:], s=2)
                ax1.plot(tslow, sinusoidal(waz,pcoef[0,i,0,0], pcoef[0,i,1,0], pcoef[0,i,2,0]), label='Ch0, phase = %i' %int(pcoef[0,i,1,0]), color='red')
                ax1.set_ylabel('Tsrc')
                ax1.set_xlabel('tslow')
                ax1.set_xlim([100, 300])
                ax1.legend(loc="upper right")

                ax2 = fig.add_subplot(412)
                ax2.scatter(tslow, data[1,:], s=2)
                ax2.plot(tslow, sinusoidal(waz,pcoef[1,i,0,0], pcoef[1,i,1,0], pcoef[1,i,2,0]), label = 'Ch1, phase= %i' %int(pcoef[1,i,1,0]), color='red')
                ax2.set_ylabel('Tsrc')
                ax2.set_xlabel('tslow')
                ax2.set_xlim([100, 300])
                ax2.legend(loc="upper right")

                ax3 = fig.add_subplot(413)
                ax3.scatter(tslow, data[2,:], s=2)
                ax3.plot(tslow, sinusoidal(waz,pcoef[2,i,0,0], pcoef[2,i,1,0], pcoef[2,i,2,0]),label='Ch3, phase=%i' %int(pcoef[2,i,1,0]), color='red')
                ax3.set_ylabel('Tsrc')
                ax3.set_xlabel('tslow')
                ax3.set_xlim([100, 300])
                ax3.legend(loc="upper right")
                         
                ax4 = fig.add_subplot(414)
                ax4.scatter(tslow, data[3,:], s=2)
                ax4.plot(tslow, sinusoidal(waz,pcoef[3,i,0,0], pcoef[3,i,1,0], pcoef[3,i,2,0]),label='Ch3, phase=%i' %int(pcoef[3,i,1,0]), color='red')  
                ax4.set_ylabel('Tsrc')
                ax4.set_xlabel('tslow')
                ax4.set_xlim([100, 300])
                ax4.legend(loc="upper right")

                fig.suptitle("Fit "+str(fl[i]))
                         
                plt.show()

                         
                debug=0
                if pcoef[j,i,1,0]==-180:
                        debug=1
                if debug:
                    clf()
                    subplot(2,1,1)
                    #plot(azslow,res)
                    plot(azslow,smooth(res,5),'g')
                    plot(azslow,sinusoidal(azslow,*fit),'r')
                    print fit
                    print fit2
                    title('%s chan:%s'%(f,j))
                    subplot(2,1,2)
                    loglog(abs(fft(sinusoidal(azslow,*fit))),'r')
                    loglog(abs(fft(res)))
                    xlim([100,200])
                    draw()
                    raw_input()
                    
        return pcoef, fit_err

                

    def TestAngFourier(self, azslow, res, k):            
            #T=27.5
            #Omega=2.*(np.pi)*k/T

            b=2.*(res*np.sin(np.deg2rad(k*azslow))).mean()
            a=2.*(res*np.cos(np.deg2rad(k*azslow))).mean()
            
            ffourier=a*np.cos(np.deg2rad(k*azslow))+b*np.sin(np.deg2rad(k*azslow))
            
            phase_rad=np.arctan2(-b,a)
            phase= math.degrees(phase_rad)

            amp=np.sqrt((a**2.)+(b**2.))

            fcos=amp*np.cos(np.deg2rad(azslow)+phase_rad)
            
            return phase_rad, amp, ffourier





    
    def findScans(self, az):
        '''
        The az reading is in degrees traveled from home position since beginning of observation.
        waz: Wrapped azimuth (using divmod( , 360)) from 0-360
        Also also enumerate each 360-degree scan and return an fs structure
        - scannum,
        - start index,
        - end index of each scan

        '''
        class stru:
            def __init__(self):
                self.num = []
                self.s = []
                self.e = []

        naz = size(az)
        (scannum, az) = divmod(az,360)
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
            yp = interp(arange(0,361), waz[idx], y[idx], right=nan, left=nan)
            D[:,i]=yp

        return D

    def filterScans(self,waz, d, fs,filttype, phase_fixed, amp_fixed, fitphase, fitamp, k):
        """
        filttype can be p0, p1, p2 p3 for poly subtraction
        filt type can be sin
        """
        if filttype[0] == 'n':
            # do nothing
            print "do nothing"
        elif filttype[0]=='p':
            # subtract poly of required order
            [d, pcoef, baseline]=self.polysub_scans(d,fs,int(filttype[1:]))
        elif filttype == 'sin':
            # do cos fit
            [d,pcoef, baseline,T] = self.sinsub_scans(waz, d, fs, fitphase, fitamp, phase_fixed, amp_fixed, k)
        elif filttype == 'skydip':
            [d,pcoef, baseline] = self.expsub_scans(waz, d, fs)   

        return  d, pcoef, baseline,T

    
    def polysub_scans(self,d,fs,porder):
        """
        """
        nscans = len(fs.s)
        pcoef=zeros([nscans,porder+1]);
        # for each 360-scan 
        y = copy(d)
        b = zeros(shape(d))
        T = zeros(shape(d))

        for i in range(nscans):
            s = fs.s[i];e=fs.e[i]
            x=arange(e-s)
            fit = polyfit(x,y[s:e],porder)
            baseline = polyval(fit,x)
            
            y[s:e] = y[s:e] - baseline
            pcoef[i,:]=fit
            b[s:e] = baseline      
            
        return y, pcoef, b
            
    def sinsub_scans(self, waz, d, fs, fitphase, fitamp, phase_fixed, amp_fixed, k):
        '''
        Slice the data in each scan and remove a sin fit to each scan
        '''
        #bounds3=((0, -360, 0), (inf, 360, inf))                                                                                   #bounds2=((0, 0), (inf,inf))
        
        nscans = len(fs.s)
        y = copy(d)  # to store the residuals
        b = zeros(shape(d)) # to store the cosine fit
        T = zeros(shape(d))
        # Assume frequency = 1 rotation in 360deg.                                                         
        # initial guesses:                                                                                                                                                  
        ph0 = -60
        
	#def cosinusoidal(x,off):
         #       return amp_fixed_new*cos(deg2rad(k*x+phase_fixed))+off

        # initialize some lists to hold the fit parameters
        if fitphase==1:
                pcoef = zeros([nscans, 3, 2])
                params0 = [1.0, ph0, mean(d)]
        else:
                pcoef = zeros([nscans, 2, 2]) # nobs x 3=Amplitudes,phase, offset x 2=values and error
                params0 = [1.0, mean(d)]
        # loop over each 360 scan
        print(nscans)
        for i in range(nscans):
                def cosinusoidal(x,off):                                                                                                                                                         
                        return amp_fixed[i]*cos(deg2rad(k*x+phase_fixed))+off
                
                #if fitamp==1:
                       # def cosinusoidal(x,off,amp):
                        #        return amp*cos(deg2rad(k*x+phase_fixed))+off
                #if fitamp==0:
                 #       def cosinusoidal(x,off):
                  #              return amp_fixed[i]*cos(deg2rad(k*x+phase_fixed))+off

                print(i)
                s = fs.s[i];e=fs.e[i]
                idx = range(s,e)
                
                fit, pcov= curve_fit(cosinusoidal, waz[idx], y[idx])#, p0=params0)
                #print("fitamp="+str(fitamp))
                if len(fit)==2:
                        fit = [fit[0],fit[1]]
                        fiterr = sqrt(diag(pcov))
                        fiterr = [fiterr[0],fiterr[1]]
                        #print("amp="+str(fit[1]))
                        #print("offs="+str(fit[0]))
                if len(fit)==1:
                        fit = [fit[0]]
                        fiterr = sqrt(diag(pcov))
                        fiterr = [fiterr[0]]
                        #print("offs="+str(fit[0]))
                        #print("amp="+str(amp_fixed[i]))
                baseline = cosinusoidal(waz[idx],*fit)
                baseline2=cosinusoidal(waz[idx],0)#baseline without removing the offset
                
                
                print("baseline")
                print(baseline)
                print("baseline2")
                print(baseline2)
                print("fit")
                print(fit)

                
                y[idx] = y[idx] - baseline2
                b[idx] = baseline
                pcoef[i,:,0]=fit
                pcoef[i,:,1]=fiterr
                T[idx]=y[idx]-baseline2
        
                
                debug=0
                if debug:
                        print '%d of %d'%(i,nscans)
                        print fit
                        clf()
                        plot(waz[idx],y[idx],'.-')
                        plot(waz[idx],baseline,'r')
                        title('%d of %d'%(i,nscans))
                        xlim([0,360])
                        raw_input()
                
        return y, pcoef, b, T





    def expsub_scans(self, waz, d, fs):
        '''
        Slice the data in each scan and remove an exp fit to each scan
        '''
        bounds = ((240, 0, -2), (280, 4, 4))
        # define the function to fit to:
        def skydip(x,Tatm,tau0,phi):
            kappa = -60
            el0 = deg2rad(55)
            c = cos(el0)
            s = sin(el0)
            Tcmb = 2.73
            eps =  phi*sin(deg2rad(x + kappa))
            DC = Tatm + ((Tcmb - Tatm))*exp(-tau0/s)
            fac = -c * eps +eps**2*(s/2+c**2/s+c**2/2.)
            MOD = DC*fac

            return DC+MOD
            
        nscans = len(fs.s)
        y = copy(d)  # to store the residuals
        b = zeros(shape(d)) # to store the cosine fit
        
        # initialize some lists to hold the fit parameters
        pcoef = zeros([nscans, 3, 2]) 

        # Assume frequency = 1 rotation in 360deg.
        # initial guesses:
        params0 = [mean(d), 2.0, 1.0]

        # loop over each 360 scan
        for i in range(nscans):
            s = fs.s[i];e=fs.e[i]
            idx = range(s,e)
            try:
                fit, pcov= curve_fit(skydip, waz[idx], y[idx], p0=params0,bounds=bounds)
            except:
                fit = [nan,nan,nan]
                fiterr = [nan,nan,nan]
            fiterr = sqrt(diag(pcov))
            baseline = skydip(waz[idx],*fit)
            print '%d of %d'%(i,nscans)
            print fit

            debug=1
            if debug:
                clf()
                plot(waz[idx],y[idx],'.-')
                plot(waz[idx],baseline,'r')
                title('%d of %d'%(i,nscans))
                xlim([0,360])
                raw_input()

            y[s:e] = y[s:e] - baseline
            b[s:e] = baseline
            pcoef[i,:,0]=fit
            pcoef[i,:,1]=fiterr
            
            # define new starting point params results of last fit.
            #params0 = fit
            
        return y, pcoef, b



    def ExtractPhaseTrend(self, phase_data): #fits the phase data with a 1st degree polynomial function and returns the coefficients (p[0]=tilt, p[1]=offset)  
        scan_n=range(len(phase_data))		
        
	p=np.zeros(2)
	p=np.polyfit(scan_n, phase_data, 1)

        return p



    def FindPhase(self,unit,season): #returns the phase for a given day of observation

        rwp = rw.reduc_wvr_pager(unit)

        start=season+"0101"
        end=season+"0601"
        
        fl = rwp.makeFileListFromData(typ='scanAz', start=start, end=end)
        pcoef, fit_err=self.testAngFit(fl)
        phase_data=pcoef[:,:,1,0]
        dt_list=rwp.dt_list
        phase=np.zeros(4)
        scan_n=range(len(phase_data[0,:]))
        fig = plt.figure()
        fig2, axs = plt.subplots(4, 1, sharex=True, tight_layout=True)
        
        #fig3 = plt.figure()
        #fig4 = plt.figure()
        
        for i in range(0,4):
                newphase = np.zeros((len(scan_n),2))
                newphase[:,0]=scan_n
                p=self.ExtractPhaseTrend(phase_data[i,:])
                newphase[:,1]=np.asarray(scan_n)*p[0]+p[1]
                phase_pickle = open("phase_season"+str(season)+"_ch"+str(i)+".pickle","wb")
                pickle.dump(newphase, phase_pickle)
                phase_pickle.close()

                axs[i].hist(phase_data[i,:], bins=200, normed=True, label='Ch %i' %i)
                axs[i].set_xlim(left=-150, right=150)
                axs[i].legend()
                axs[i].set_xlabel('phase')
                
                self.PlotPhaseFit(i, dt_list, phase_data[i,:], newphase[:,1], fig)
                #self.PlotPhaseErr(i, dt_list, fit_err, fig3)
                

        fig.suptitle("Phase fit "+str(season))
        fig2.suptitle("Phase Histogram "+str(season))
        #fig3.suptitle("Phase Err "+str(season))
        #fig4.suptitle("Phase Err Histogram "+str(season))
        
        
        plt.show()

        return fit_err


    def FindPhaseFourier(self,unit,season):

	rwp = rw.reduc_wvr_pager(unit)
        start=season+"0101"
        end=season+"0601"
        fl = rwp.makeFileListFromData(typ='scanAz', start=start, end=end)
        nobs = size(fl)
        phase_data = zeros([4, nobs])
        phase_data_double = zeros([4, nobs])
        offs_data = zeros([4, nobs])
        amp_data = zeros([4, nobs])                                     
        amp_data_double = zeros([4, nobs])
        scan_n=range(nobs)
        check_plot=0

        if check_plot==1:
                fig,axs = plt.subplots(4, 1, sharex=True, tight_layout=True)
                fig2,axs2 = plt.subplots(4, 1, sharex=True, tight_layout=True)
                fig3,axs3 = plt.subplots(4, 1, sharex=True, tight_layout=True)

        for j,f in enumerate(fl):

                dt_list=rwp.dt_list
                utslow,tslow,d,azslow,elslow,tsrc = self.wvrR.readSlowFile(f)
                waz,fs =  self.findScans(azslow)
                for i in range(4):
                        
                        res, pcoef0, baseline = self.filterScans(waz, tsrc[:,i], fs, 'p0', phase_fixed=0., amp_fixed=0, fitphase=0, fitamp=1, k=2)
                        offs_data[i,j]=baseline.mean()
                        phase_data_rad_double, ampl_double, ffourier_double = self.TestAngFourier(azslow, res, k=2)
                        amp_data_double[i,j]=ampl_double
                        phase_data_double[i,j]=np.rad2deg(phase_data_rad_double)

                        phase_data_rad, ampl, ffourier = self.TestAngFourier(azslow, res-ffourier_double, k=1)
                        amp_data[i,j]=ampl
                        phase_data[i,j]=np.rad2deg(phase_data_rad)
                        
                        #res, pcoef0, baseline = self.filterScans(waz, tsrc[:,i], fs, 'p0', phase_fixed=0., fitphase=0, k=1)
                        #print(baseline.mean())
                        #offs_data[i,j]=baseline.mean()
                        #phase_rad, ampl, ffourier = self.TestAngFourier(azslow, res, k=1)
                        #amp_data[i,j]=ampl
                        #if phase_rad <0:
                         #       phase_rad=(np.pi/2.)+((np.pi/2.)+phase_rad)  
                        #phase_data[i,j]=phase_rad
                        #phase_data[i,j]=np.rad2deg(phase_rad)
                        #double_mod=res-ffourier
                        #phase_data_double_rad, ampl_double, ffourier_double = self.TestAngFourier(azslow, double_mod, k=2)
                        #amp_data_double[i,j]=ampl_double
                        #phase_data_double[i,j]=np.rad2deg(phase_data_double_rad)

                        if phase_data_double[i,j]<0:
                                phase_data_double[i,j]=360.+phase_data_double[i,j]
                        
                           
                        if j==20 and check_plot==1:
                                f=(np.pi)/30.
                                axs[i].plot(azslow, res,label='phase=%i'%np.rad2deg(phase_rad))
                                axs[i].plot(azslow, ffourier, label='Ch %i' %i, color='red')
                                axs[i].set_xlim([0, 2000])
                                axs[i].set_ylabel('Tsrc-T0')
                                axs[i].set_xlabel('azslow')
                                axs[i].legend(loc="upper right")
                                fig.suptitle("Data")

                                axs2[i].plot(azslow, double_mod, label='phase=%i'%np.rad2deg(phase_data_double_rad))
                                axs2[i].plot(azslow, ffourier_double, label='Ch %i' %i, color='red')
                                axs2[i].set_xlim([0, 2000])
                                axs2[i].set_ylabel('Tsrc - T0')
                                axs2[i].set_xlabel('azslow')
                                axs2[i].legend(loc="upper right")
                                fig2.suptitle("Data - Single Mod ")

                                axs3[i].plot(azslow, double_mod-ffourier_double, label='Ch %i' %i) 
                                axs3[i].set_xlim([0, 2000])
                                axs3[i].set_ylabel('Tsrc - T0')
                                axs3[i].set_xlabel('azslow')
                                axs3[i].legend(loc="upper right")
                                fig3.suptitle("Data - Single and Double Mod")
                                
                                plt.show()

        return(phase_data, phase_data_double,amp_data, amp_data_double, offs_data, dt_list)


    def FindAmpDoubleMod(self, season, offs):

            f_amp = open('amp_data_double_'+season+'.txt','r')
            amp_data_double = pickle.load(f_amp)

            f_offs = open('offs_data_'+season+'.txt','r')
            offs_data  = pickle.load(f_offs)

            mask=amp_data_double[3,:]>0.1 #just for 2018

            l=len(amp_data_double[0,mask])
            
            amp_data_all=np.zeros(4*l)
            offs_data_all=np.zeros(4*l)

            for i in range (0,4):
                    amp_data_all[i*l:(i+1)*l]=amp_data_double[i,mask]
                    offs_data_all[i*l:(i+1)*l]=offs_data[i,mask]

            p=np.polyfit(offs_data_all, amp_data_all, deg=1)

            amp_double= p[0]*offs+p[1]

            return amp_double,p
            
    def FindAmpSingleMod(self, season, offs):

            f_amp = open('amp_data_'+season+'.txt','r')
            amp_data = pickle.load(f_amp)

            f_offs = open('offs_data_'+season+'.txt','r')
            offs_data  = pickle.load(f_offs)

            mask=amp_data[0,:]>0.28 #just for 2018                                                                                                                                         

            l=len(amp_data[0,mask])

            amp_data_all=np.zeros(4*l)
            offs_data_all=np.zeros(4*l)

            for i in range (0,4):
                    amp_data_all[i*l:(i+1)*l]=amp_data[i,mask]
                    offs_data_all[i*l:(i+1)*l]=offs_data[i,mask]

            def toy_model(x,T0,a):
                    tau=-np.log(1.-(x/T0))
                    T=T0*(1.-np.exp(-tau))
                    T2=T0*(1-np.exp(-tau*(1.+a)))
                    return T2-T 
            
            fit, pcov= curve_fit(toy_model, offs_data_all, amp_data_all, p0=[250.,0.01])

            #computes the value of the amp given the input offset
            tau=-np.log(1.-(offs/fit[0]))
            T=fit[0]*(1.-np.exp(-tau))
            T2=fit[0]*(1-np.exp(-tau*(1.+fit[1])))

            amp_single=T2-T

            return amp_single,fit
            


    
    def FitandPlot(self, phase_data, amp_data, offs_data, season, dt_list):

        scan_n=range(size(phase_data[0,:]))
        p=zeros([4,2])
        phase_season=np.zeros(4)
        amp_season=np.zeros(4)
        offs_season=np.zeros(4)
        fig=plt.figure()
        fig2, axs = plt.subplots(4, 1, sharex=True, tight_layout=True)

        for i in range (0,4):
                p[i,:]=self.ExtractPhaseTrend(phase_data[i,:])
                newphase=np.asarray(scan_n)*p[i,0]+p[i,1]

                phase_pickle = open("phase_double_season"+str(season)+"_ch"+str(i)+"_Fourier.pickle","wb")
                pickle.dump(newphase, phase_pickle)
                phase_pickle.close()

                self.PlotPhaseFit(i, dt_list, phase_data[i,:], newphase, fig, p)
                counts, bin_edges=np.histogram(phase_data[i,:],bins=200, density=True)
                axs[i].hist(phase_data[i,:],bins=200,density=True, label='Ch %i' %i)
                axs[i].set_xlim(left=100, right=300)
                axs[i].legend()

                phase_index=(np.where(counts==np.max(counts))[0])[0]
                phase_season[i]=(bin_edges[phase_index]+bin_edges[phase_index+1])/2.

        fig.suptitle("Phase fit Fourier"+str(season))
        fig2.suptitle("Phase fit Fourier Histogram"+str(season))
        plt.show()


        fig=plt.figure()
        fig2, axs = plt.subplots(4, 1, sharex=True, tight_layout=True)
        
        for i in range (0,4):
                self.PlotAmpFit(i, dt_list, amp_data[i,:], fig)
                counts_amp, bin_edges_amp=np.histogram(amp_data[i,:], bins=50, density=True)
                n, bins, patches = axs[i].hist(amp_data[i,:], bins=50, density=True, label='Ch %i' %i)
                axs[i].set_xlim(left=0, right=0.75)
                axs[i].legend()

                amp_index=(np.where(counts_amp==np.max(counts_amp))[0])[0]
                amp_season[i]=(bin_edges_amp[amp_index]+bin_edges_amp[amp_index+1])/2.

        fig.suptitle("Amplitude fit Fourier"+str(season))
        fig2.suptitle("Amplitude fit Fourier Histogram"+str(season))
        plt.show()


        fig=plt.figure()
        fig2, axs = plt.subplots(4, 1, sharex=True, tight_layout=True)
        
        for i in range (0,4):
                self.PlotOffsFit(i, dt_list, offs_data[i,:], fig)
                counts_offs, bin_edges_offs=np.histogram(offs_data[i,:], bins=50, density=True)
                n, bins, patches = axs[i].hist(offs_data[i,:], bins=50, density=True, label='Ch %i' %i)
                axs[i].set_xlim(left=0, right=0.75)
                axs[i].legend()
                offs_index=(np.where(counts_offs==np.max(counts_offs))[0])[0]
                offs_season[i]=(bin_edges_offs[offs_index]+bin_edges_offs[offs_index+1])/2.

        fig.suptitle("Offset "+str(season))
        fig2.suptitle("Offset Histogram"+str(season))
        plt.show()

        
        return phase_season, amp_season, offs_season
        



    def PlotPhaseErr(self, i, dt_list, fit_err, fig):
        years = mdates.YearLocator()
        months = mdates.MonthLocator()
        days = mdates.DayLocator()
            
        ax="ax"+str(i+1)
        st="41"+str(i+1)
        ax = fig.add_subplot(int(st))
        ax.scatter(dt_list, fit_err[i,:,1], s=2, label='Ch %i' %i, color='red')
        ax.set_ylabel('phase_err')
        ax.set_ylim(top=+10)
        ax.set_ylim(bottom=-10)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        ax.tick_params(axis='x', rotation=45)


        
    def PlotPhaseFit(self,i, dt_list, phase_data, newphase, fig, p):

        years = mdates.YearLocator()
        months = mdates.MonthLocator()
        days = mdates.DayLocator()

        ax="ax"+str(i+1)
        st="41"+str(i+1)
        ax = fig.add_subplot(int(st))
        np.set_printoptions(precision=3)
        ax.scatter(dt_list, phase_data,label='Ch %i' %i, s=2)
        ax.plot(dt_list, newphase, label='Offset= %i\n' %int(p[i,1])+'tilt=%f' %p[i,0], color='red')
        ax.set_ylabel('phase')
        ax.set_ylim(top=+300)
        ax.set_ylim(bottom=100)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        ax.tick_params(axis='x', rotation=45)
        if i==3:
                myFmt = DateFormatter('%m/%Y') 
                ax.xaxis.set_major_formatter(myFmt)
        else:
                myFmt = DateFormatter('')
                ax.xaxis.set_major_formatter(myFmt)


    def PlotAmpFit(self,i, dt_list, amp_data, fig):

        years = mdates.YearLocator()
        months = mdates.MonthLocator()
        days = mdates.DayLocator()

        ax="ax"+str(i+1)
        st="41"+str(i+1)
        ax = fig.add_subplot(int(st))
        np.set_printoptions(precision=3)
        ax.scatter(dt_list, amp_data,label='Ch %i' %i, s=2)
        ax.set_ylabel('amplitude')
        #ax.set_ylim(top=+300)
        #ax.set_ylim(bottom=100)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        ax.tick_params(axis='x', rotation=45)
        if i==3:
                myFmt = DateFormatter('%m/%Y')
                ax.xaxis.set_major_formatter(myFmt)
        else:
                myFmt = DateFormatter('')
		ax.xaxis.set_major_formatter(myFmt)

                
    def PlotOffsFit(self,i, dt_list, offs_data, fig):

        years = mdates.YearLocator()
        months = mdates.MonthLocator()
        days = mdates.DayLocator()

        ax="ax"+str(i+1)
        st="41"+str(i+1)
        ax = fig.add_subplot(int(st))
        np.set_printoptions(precision=3)
        ax.scatter(dt_list, offs_data,label='Ch %i' %i, s=2)
        ax.set_ylabel('offset')
        #ax.set_ylim(top=+300)                                                                             
        #ax.set_ylim(bottom=100)                                                                           
        ax.legend(loc="upper right")
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        ax.tick_params(axis='x', rotation=45)
        if i==3:
                myFmt = DateFormatter('%m/%Y')
                ax.xaxis.set_major_formatter(myFmt)
        else:
                myFmt = DateFormatter('')
                ax.xaxis.set_major_formatter(myFmt)
                

    def PlotHist(self, i, phase_data, fig):
            
        #ax="ax"+str(i+1)
        #st="41"+str(i+1)
        ax1 = fig.add_subplot(int(411))
        ax1.hist(phase_data[0,:], bins=200, normed=True, label='Ch 1')
        ax1.legend()




    def PlotPhaseRes(self,i, dt_list, res, fig):

        years = mdates.YearLocator()
        months = mdates.MonthLocator()
        days = mdates.DayLocator()

        ax="ax"+str(i+1)
        st="41"+str(i+1)
        ax = fig.add_subplot(int(st))
        ax.scatter(dt_list, phase_data, s=2)
        ax.plot(dt_list, newphase, label='Ch %i' %i, color='red')
        ax.set_ylabel('phase')
        ax.set_ylim(top=-50)
        ax.set_ylim(bottom=-200)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        ax.tick_params(axis='x', rotation=45)
        if i==3:
                myFmt = DateFormatter('%m/%Y')
                ax.xaxis.set_major_formatter(myFmt)
        else:
                myFmt = DateFormatter('')
                ax.xaxis.set_major_formatter(myFmt)









        
#from scipy.optimize import fsolve
#import math
#
#def equations(p):
#    x, y = p
#    return (x+y**2-4, math.exp(x) + x*y - 3)

#x, y =  fsolve(equations, (1, 1))

#print equations((x, y))

#nsolve([x+y**2-4, exp(x)+x*y-3], [x, y], [1, 1])




