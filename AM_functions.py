import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import scipy.interpolate as sint
import numpy as np
import pickle
import pylab as pl
from matplotlib import dates
import datetime
import scipy.optimize as sp
from pathlib import Path
import plotElnod as pE
import plotAzscan as pA
import BK_analysis as BK_an
from time import perf_counter
import glob
import re as r

re=pE.ReadElnod()
raz=pA.ReadAzscan()
bk=BK_an.extract_ts()

class AM_functions(object):

    def __init__(self, unit=None, verb=True):

        '''


        '''
        self.LO=91.65
        self.IF=2*self.LO

        self.ch=[1.05,1.9,3.2,6.1]
        self.freqs=[self.IF-self.ch[3],self.IF-self.ch[2],self.IF-self.ch[1],self.IF-self.ch[0],self.IF+self.ch[0],self.IF+self.ch[1],self.IF+self.ch[2],self.IF+self.ch[3]]
        self.path_to_all = ''

    def plot_out_file(self, filename):
        with open('am_software/cfg_files/'+filename) as f:
            lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
            data = np.loadtxt(lines)
        pl.plot(data[:,0],-np.log(data[:,1]), label='Tau[neper]')
        pl.plot(data[:,0],data[:,2], label='T')
        pl.xlim(170,200)
        pl.close()

    def create_am_datfile(self, filename, path_to_data='wvr1_data/', template='SPole_annual_50.amc', spline=0, showplots=2, write_dat=1):
        if not os.path.exists(f'am_datfiles/'+template[:-4]):
            os.makedirs(f'am_datfiles/'+template[:-4])
        if not os.path.exists(f'am_datfiles/'+template[:-4]+'/'+filename[:-4]):
            os.makedirs(f'am_datfiles/'+template[:-4]+'/'+filename[:-4])

        data=re.read_elnod_fast(path_to_data+filename[:-9]+'/'+filename)

        path='am_datfiles/'+template[:-4]+'/'+filename[:-4]
        if spline==0:
            pickle_fn=path+'/'+filename[:-4]+'_pickle_temps.txt'
        elif spline==1:
            pickle_fn=path+'/'+filename[:-4]+'_pickle_temps_red.txt'
            pickle_betas=path+'/'+filename[:-4]+'_pickle_betas.txt'
        elif spline==2:
            pickle_fn=path+'/'+filename[:-4]+'_pickle_temps_corr.txt'
            pickle_betas=path+'/'+filename[:-4]+'_pickle_betas.txt'
            pickle_deltas=path+'/'+filename[:-4]+'_pickle_deltas.txt'
            pickle_z=path+'/'+filename[:-4]+'_pickle_z_offs.txt'

        pickle_Temps= {'am_dat_filename': [], 'El':[], 'dir':[], 'T0':[], 'T1':[], 'T2':[], 'T3':[]}
        z_offs_elnod={'z_offs':[],'z_offs_err':[]}
        betas=np.zeros(4)

        Npoints=20

        if spline==0:
            open_endstring='.dat'
            with open(path_to_data+filename[:-4]+open_endstring,'w+') as f:
                time=data[:,0]
                T0=data[:,1]
                T1=data[:,2]
                T2=data[:,3]
                T3=data[:,4]
                el=data[:,5]

                time=time[el>0]
                T0=T0[el>0]
                T1=T1[el>0]
                T2=T2[el>0]
                T3=T3[el>0]
                el=el[el>0]

                previous_el=0
                dir=np.zeros(len(el))
                dir_col=['Empty']*len(el)
                for j in range (len(el)):
                    print('t=', time[j])
                    if (el[j]>=previous_el):
                        dir[j]=1
                        dir_col[j]='r'
                    elif (el[j]<previous_el):
                        dir[j]=-1
                        dir_col[j]='b'
                    previous_el=el[j]

                T0_up=T0[np.where(dir==1)]
                T0_down=T0[np.where(dir==-1)]
                T1_up=T1[np.where(dir==1)]
                T1_down=T1[np.where(dir==-1)]
                T2_up=T2[np.where(dir==1)]
                T2_down=T2[np.where(dir==-1)]
                T3_up=T3[np.where(dir==1)]
                T3_down=T3[np.where(dir==-1)]

                el_up=el[np.where(dir==1)]
                el_down=el[np.where(dir==-1)]

                if write_dat==1:
                    for file in os.scandir('am_datfiles/'+template[:-4]+'/'+filename[:-4]):
                        if file.name.endswith("_red.dat"):
                            os.unlink(file.path)

                    dir_string=['None']*len(time)
                    for i in range(len(time)):
                        if dir[i]==1:
                            dir_string[i]='up'
                        elif dir[i]==-1:
                            dir_string[i]='dn'
                        else:
                            print('Direction not identified.\ndir='+str(dir))

                        dat_fn='am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el[i],1))+'_'+dir_string[i]+'_am'+open_endstring

                        with open('am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el[i],1))+'_'+dir_string[i]+'_am'+open_endstring,'w+') as f:
                            print('writing on'+ 'am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el[i],1))+'_'+dir_string[i]+'_am'+open_endstring)
                            f.write('1.25 1.50 {0}\n'.format(T0[i]))
                            f.write('3.25 2.50 {0}\n'.format(T1[i]))
                            f.write('5.5 2.00 {0}\n'.format(T2[i]))
                            f.write('7.25 1.50 {0}'.format(T3[i]))
                            print('writing done.')
                        f.close()

                    pickle_Temps['am_dat_filename'].append(dat_fn)
                    pickle_Temps['El'].append(round(el[i],1))
                    pickle_Temps['dir'].append(dir[i])
                    pickle_Temps['T0'].append(T0[i])
                    pickle_Temps['T1'].append(T1[i])
                    pickle_Temps['T2'].append(T2[i])
                    pickle_Temps['T3'].append(T3[i])

                f = open(pickle_fn,'wb')
                pickle.dump(pickle_Temps, f)
                f.close()


        elif spline==1:

            open_endstring='_red.dat'
            with open(path_to_data+filename[:-4]+open_endstring,'w+') as f:
                time=data[:,0]
                T0=data[:,1]
                T1=data[:,2]
                T2=data[:,3]
                T3=data[:,4]
                el=data[:,5]

                time=time[el>0]
                T0=T0[el>0]
                T1=T1[el>0]
                T2=T2[el>0]
                T3=T3[el>0]
                el=el[el>0]

                T0[np.where(T0 == np.nan)]=0.
                T1[np.where(T0 == np.nan)]=0.
                T2[np.where(T0 == np.nan)]=0.
                T3[np.where(T0 == np.nan)]=0.

                previous_el=0
                dir=np.zeros(len(el))
                for j in range (len(el)):
                    print('t=', time[j])
                    if (el[j]>=previous_el):
                        dir[j]=1
                    elif (el[j]<previous_el):
                        dir[j]=-1
                    previous_el=el[j]

                T0_up=T0[np.where(dir==1)]
                T0_down=T0[np.where(dir==-1)]
                T1_up=T1[np.where(dir==1)]
                T1_down=T1[np.where(dir==-1)]
                T2_up=T2[np.where(dir==1)]
                T2_down=T2[np.where(dir==-1)]
                T3_up=T3[np.where(dir==1)]
                T3_down=T3[np.where(dir==-1)]

                el_up=el[np.where(dir==1)]
                el_down=el[np.where(dir==-1)]

                el_write=np.linspace(15,90,20)

                T0_coef_up=np.polyfit(el_up, T0_up, 7)
                T0_fit_up=np.poly1d(T0_coef_up)
                T0_write_up=T0_fit_up(el_write)
                T1_coef_up=np.polyfit(el_up, T1_up, 7)
                T1_fit_up=np.poly1d(T1_coef_up)
                T1_write_up=T1_fit_up(el_write)
                T2_coef_up=np.polyfit(el_up, T2_up, 7)
                T2_fit_up=np.poly1d(T2_coef_up)
                T2_write_up=T2_fit_up(el_write)
                T3_coef_up=np.polyfit(el_up, T3_up, 7)
                T3_fit_up=np.poly1d(T3_coef_up)
                T3_write_up=T3_fit_up(el_write)

                T0_coef_down=np.polyfit(el_down, T0_down, 7)
                T0_fit_down=np.poly1d(T0_coef_down)
                T0_write_down=T0_fit_down(el_write)
                T1_coef_down=np.polyfit(el_down, T1_down, 7)
                T1_fit_down=np.poly1d(T1_coef_down)
                T1_write_down=T1_fit_down(el_write)
                T2_coef_down=np.polyfit(el_down, T2_down, 7)
                T2_fit_down=np.poly1d(T2_coef_down)
                T2_write_down=T2_fit_down(el_write)
                T3_coef_down=np.polyfit(el_down, T3_down, 7)
                T3_fit_down=np.poly1d(T3_coef_down)
                T3_write_down=T3_fit_down(el_write)

                fig, axes = plt.subplots(2, 2, figsize=(18, 10))
                axes[0,0].scatter(el_up, T0_up, marker='^', s=3, alpha=0.5, c='red', label='Moving up')
                axes[0,0].plot(sorted(el_up), T0_fit_up(sorted(el_up)), c='orange', label='Fit up')
                axes[0,0].scatter(el_down, T0_down, marker='v', s=3, alpha=0.5, c='blue', label='Moving down')
                axes[0,0].plot(sorted(el_down), T0_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[0,0].legend()
                axes[0,0].set_title('Channel 0')

                axes[0,1].scatter(el_up, T1_up, marker='^', s=3, alpha=0.5, c='red', label='moving up')
                axes[0,1].plot(sorted(el_up), T1_fit_up(sorted(el_up)),  c='orange', label='Fit up')
                axes[0,1].scatter(el_down, T1_down, marker='v', s=3, alpha=0.5, c='blue', label='moving down')
                axes[0,1].plot(sorted(el_down), T1_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[0,1].legend()
                axes[0,1].set_title('Channel 1')

                axes[1,0].scatter(el_up, T2_up, marker='^', s=3, alpha=0.5, c='red', label='moving up')
                axes[1,0].plot(sorted(el_up), T2_fit_up(sorted(el_up)), c='orange', label='Fit up')
                axes[1,0].scatter(el_down, T2_down, marker='v', s=3, alpha=0.5, c='blue', label='moving down')
                axes[1,0].plot(sorted(el_down), T2_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2')

                axes[1,1].scatter(el_up, T3_up, marker='^', s=3, alpha=0.5, c='red', label='moving up')
                axes[1,1].plot(sorted(el_up), T3_fit_up(sorted(el_up)), c='orange', label='Fit up')
                axes[1,1].scatter(el_down, T3_down, marker='v', s=3, alpha=0.5, c='blue', label='moving down')
                axes[1,1].plot(sorted(el_down), T3_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3')

                fig.suptitle('Raw T readings - '+filename[:-4]+'\nFitted with 7th order Polynomial')
                pl.savefig('Output/Backlash_Fit/'+filename[:-4]+'_Raw_T_Readings'+'_'+str(Npoints)+'.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #to fit the backlash
                el_smooth=np.linspace(15, 90, 200)

                def poly_diff(x, backlash, fit_up, fit_down):
                    print('iter')
                    return fit_up(x)-fit_up(x+backlash)
                def err(backlash, x, y, fit_up, fit_down):
                    return poly_diff(x, backlash, fit_up, fit_down) - y
                def err_global(p, x, y1, y2, y3, y4):
                    err1 = err(p, x, y1, T0_fit_up, T0_fit_down)
                    err2 = err(p, x, y2, T1_fit_up, T1_fit_down)
                    err3 = err(p, x, y3, T2_fit_up, T2_fit_down)
                    err4 = err(p, x, y4, T3_fit_up, T3_fit_down)

                    return np.concatenate((err1, err2, err3, err4))

                p=[3]

                ydata1=T0_fit_up(el_smooth)-T0_fit_down(el_smooth)
                ydata2=T1_fit_up(el_smooth)-T1_fit_down(el_smooth)
                ydata3=T2_fit_up(el_smooth)-T2_fit_down(el_smooth)
                ydata4=T3_fit_up(el_smooth)-T3_fit_down(el_smooth)

                p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global, p, args=(el_smooth, ydata1, ydata2, ydata3, ydata4), full_output=1)

                print('p_best=',p_best)
                p_err=np.sqrt(np.diag(cov_x))

                fig, axes = plt.subplots(2, 2, figsize=(18, 10))
                axes[0, 0].plot(el_smooth, T0_fit_up(el_smooth)-T0_fit_down(el_smooth), label='diff_data')
                axes[0, 0].plot(el_smooth, poly_diff(el_smooth, *p_best, T0_fit_up, T0_fit_down), label='diff_fit')
                axes[0, 0].axhline(y=0, color='k', linestyle='--')
                axes[0, 0].set_title('Channel 0')
                axes[0, 0].legend(loc='upper right')

                axes[0, 1].plot(el_smooth, T1_fit_up(el_smooth)-T1_fit_down(el_smooth), label='diff_data')
                axes[0, 1].plot(el_smooth, poly_diff(el_smooth, *p_best, T1_fit_up, T1_fit_down), label='diff_fit')
                axes[0, 1].axhline(y=0, color='k', linestyle='--')
                axes[0, 1].set_title('Channel 1')
                axes[0, 1].legend(loc='upper right')

                axes[1, 0].plot(el_smooth, T2_fit_up(el_smooth)-T2_fit_down(el_smooth), label='diff_data')
                axes[1, 0].plot(el_smooth, poly_diff(el_smooth, *p_best, T2_fit_up, T2_fit_down), label='diff_fit')
                axes[1, 0].axhline(y=0, color='k', linestyle='--')
                axes[1, 0].set_title('Channel 2')
                axes[1, 0].legend(loc='upper right')

                axes[1, 1].plot(el_smooth, T3_fit_up(el_smooth)-T3_fit_down(el_smooth), label='diff_data')
                axes[1, 1].plot(el_smooth, poly_diff(el_smooth, *p_best, T3_fit_up, T3_fit_down), label='diff_fit')
                axes[1, 1].axhline(y=0, color='k', linestyle='--')
                axes[1, 1].set_title('Channel 3')
                axes[1, 1].legend(loc='upper right')

                fig.suptitle(filename[:-4]+' Backlash Fit \n Backlash[deg]='+str(round(p_best[0],3))+'±'+str(round(p_err[0],3)))
                pl.savefig('Output/Backlash_Fit/'+filename[:-4]+'_line_diff_fit'+'_'+str(Npoints)+'.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #Backlash correction
                fig, axes = plt.subplots(2, 2, figsize=(18, 10))
                axes[0,0].scatter(el_smooth, T0_fit_up(el_smooth), s=3, c='r', alpha=0.5, label='moving up')
                axes[0,0].scatter(el_smooth, T0_fit_down(el_smooth), s=3, c='b', alpha=0.5, label='moving down')
                axes[0,0].plot(el_smooth, T0_fit_up(el_smooth), c='r', label='moving up')
                axes[0,0].plot(el_smooth, T0_fit_down(el_smooth-p_best),c='b', label='moving down - Corrected Backlash[deg]={:.3f}'.format(p_best[0]))
                axes[0,0].legend()
                axes[0,0].set_title('Channel 0')

                axes[0,1].scatter(el_smooth, T1_fit_up(el_smooth), s=3, c='r', alpha=0.5, label='moving up')
                axes[0,1].scatter(el_smooth, T1_fit_down(el_smooth), s=3, c='b',alpha=0.5, label='moving down')
                axes[0,1].plot(el_smooth, T1_fit_up(el_smooth), c='r', label='moving up')
                axes[0,1].plot(el_smooth, T1_fit_down(el_smooth-p_best),c='b', label='moving down - Corrected Backlash[deg]={:.3f}'.format(p_best[0]))
                axes[0,1].legend()
                axes[0,1].set_title('Channel 1')

                axes[1,0].scatter(el_smooth, T2_fit_up(el_smooth), s=3, c='r', alpha=0.5, label='moving up')
                axes[1,0].scatter(el_smooth, T2_fit_down(el_smooth), s=3, c='b', alpha=0.5, label='moving down')
                axes[1,0].plot(el_smooth, T2_fit_up(el_smooth), c='r', label='moving up')
                axes[1,0].plot(el_smooth, T2_fit_down(el_smooth-p_best),c='b', label='moving down - Corrected Backlash[deg]={:.3f}'.format(p_best[0]))
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2')

                axes[1,1].scatter(el_smooth, T3_fit_up(el_smooth), s=3, c='r', alpha=0.5, label='moving up')
                axes[1,1].scatter(el_smooth, T3_fit_down(el_smooth), s=3, c='b', alpha=0.5, label='moving down')
                axes[1,1].plot(el_smooth, T3_fit_up(el_smooth), c='r', label='moving up')
                axes[1,1].plot(el_smooth, T3_fit_down(el_smooth-p_best),c='b', label='moving down - Corrected Backlash[deg]={:.3f}'.format(p_best[0]))
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3')

                pl.suptitle('Backlash-corrected data - ' + filename[:-4])
                pl.savefig('Output/Backlash_Fit/'+filename[:-4]+'_Backlash-corrected_data'+'_'+str(Npoints)+'.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #Creating a single array for rad moving up&down after correcting Backlash
                el_bl_corrected=np.concatenate((el_smooth, el_smooth), axis=None)
                T0_bl_corrected=np.concatenate((T0_fit_up(el_smooth), T0_fit_down(el_smooth-p_best)), axis=None)
                T1_bl_corrected=np.concatenate((T1_fit_up(el_smooth), T1_fit_down(el_smooth-p_best)), axis=None)
                T2_bl_corrected=np.concatenate((T2_fit_up(el_smooth), T2_fit_down(el_smooth-p_best)), axis=None)
                T3_bl_corrected=np.concatenate((T3_fit_up(el_smooth), T3_fit_down(el_smooth-p_best)), axis=None)


                el_write=np.linspace(15,90, Npoints)

                zenith_ch=np.zeros(4)
                zenith_corr=np.zeros(4)

                #Finding the Zenith
                def line_shape_1lay(x, p):
                    tau=np.exp(-p[1]/np.sin((x+p[2])*np.pi/180.))
                    T_meas=p[0]*(1.-tau)
                    return T_meas

                def err_lineshape(p, x, y):
                    return line_shape_1lay(x, p) - y

                def err_global_lineshape(p, x, y1, y2, y3, y4):
                    p1=p[0], p[1], p[5]
                    p2=p[0], p[2], p[5]
                    p3=p[0], p[3], p[5]
                    p4=p[0], p[4], p[5]

                    err1 = err_lineshape(p1, x, y1)
                    err2 = err_lineshape(p2, x, y2)
                    err3 = err_lineshape(p3, x, y3)
                    err4 = err_lineshape(p4, x, y4)

                    return np.concatenate((err1, err2, err3, err4))

                p=[200, 0.6, 0.25, 0.12, 0.05, 3]

                p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global_lineshape, p, args=(el_bl_corrected, T0_bl_corrected, T1_bl_corrected, T2_bl_corrected, T3_bl_corrected), full_output=1)

                print('p_best=',p_best)
                p_err=np.sqrt(np.diag(cov_x))

                p1=p_best[0], p_best[1], p_best[5]
                p2=p_best[0], p_best[2], p_best[5]
                p3=p_best[0], p_best[3], p_best[5]
                p4=p_best[0], p_best[4], p_best[5]

                fig, axes = plt.subplots(2, 2, figsize=(18, 10))
                string0='T_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nBeta={:.2f}'.format(round(p_best[1],2))+'±{:.2f}'.format(round(p_err[1],2))+'\nZenith_corr={:.2f}'.format(round(p[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                print(string0)
                pos0y=T0_bl_corrected.max()-(T0_bl_corrected.max()-T0_bl_corrected.min())/3.
                axes[0,0].text(60.,  pos0y, string0, fontsize=10)
                axes[0,0].scatter(el_bl_corrected, T0_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[0,0].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p1), s=3, c='r', label='Fit - Single slab atmospheric model')
                axes[0,0].set_ylim()
                axes[0,0].legend()
                axes[0,0].xlabel('El[deg]')
                axes[0,0].set_title('Channel 0')

                string1='T_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nBeta={:.2f}'.format(round(p_best[2],2))+'±{:.2f}'.format(round(p_err[2],2))+'\nZenith_corr={:.2f}'.format(round(p[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                print(string1)
                pos1y=T1_bl_corrected.max()-(T1_bl_corrected.max()-T1_bl_corrected.min())/3.
                axes[0,1].text(60.,  pos1y, string1, fontsize=10)
                axes[0,1].scatter(el_bl_corrected, T1_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[0,1].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p2), s=3, c='r', label='Fit - Single slab atmospheric model')
                axes[0,1].legend()
                axes[0,1].xlabel('El[deg]')
                axes[0,1].set_title('Channel 1')

                string2='T_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nBeta={:.2f}'.format(round(p_best[3],2))+'±{:.2f}'.format(round(p_err[3],2))+'\nZenith_corr={:.2f}'.format(round(p[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                print(string2)
                pos2y=T2_bl_corrected.max()-(T2_bl_corrected.max()-T2_bl_corrected.min())/3.
                axes[1,0].text(60.,  pos2y, string2, fontsize=10)
                axes[1,0].scatter(el_bl_corrected, T2_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[1,0].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p3), s=3, c='r', label='Fit - Single slab atmospheric model')
                axes[1,0].xlabel('El[deg]')
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2')

                string3='T_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nBeta={:.2f}'.format(round(p_best[4],2))+'±{:.2f}'.format(round(p_err[4],2))+'\nZenith_corr={:.2f}'.format(round(p[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                print(string3)
                pos3y=T3_bl_corrected.max()-(T3_bl_corrected.max()-T3_bl_corrected.min())/3.
                axes[1,1].text(60.,  pos3y, string3, fontsize=10)
                axes[1,1].scatter(el_bl_corrected, T3_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[1,1].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p4), s=3, c='r', label='Fit - Single slab atmospheric model')
                axes[1,1].xlabel('El[deg]')
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3')

                pl.suptitle('Fit to one slab atmospheric model - on backlash-corrected data\n' + filename[:-4])
                pl.savefig('Output/Backlash_Fit/'+filename[:-4]+'_SingleSlabModelFit_onBacklash-corrected_data'+'_'+str(Npoints)+'.png')

                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                p1_corr=p_best[0], p_best[1], 0.
                p2_corr=p_best[0], p_best[2], 0.
                p3_corr=p_best[0], p_best[3], 0.
                p4_corr=p_best[0], p_best[4], 0.

                betas=[p_best[1], p_best[2], p_best[3], p_best[4]]

                T0_singslab=line_shape_1lay(el_bl_corrected, p1)
                T0_singslab_corr=line_shape_1lay(el_bl_corrected, p1_corr)
                T0_write=line_shape_1lay(el_write, p1_corr)
                T1_singslab=line_shape_1lay(el_bl_corrected, p2)
                T1_singslab_corr=line_shape_1lay(el_bl_corrected, p2_corr)
                T1_write=line_shape_1lay(el_write, p2_corr)
                T2_singslab=line_shape_1lay(el_bl_corrected, p3)
                T2_singslab_corr=line_shape_1lay(el_bl_corrected, p3_corr)
                T2_write=line_shape_1lay(el_write, p3_corr)
                T3_singslab=line_shape_1lay(el_bl_corrected, p4)
                T3_singslab_corr=line_shape_1lay(el_bl_corrected, p4_corr)
                T3_write=line_shape_1lay(el_write, p4_corr)



                fig, axes = plt.subplots(2, 2, figsize=(18, 10))
                axes[0,0].scatter(el_bl_corrected, T0_singslab, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[0,0].scatter(el_bl_corrected, T0_singslab_corr, s=1, c='c', alpha=0.5, label='After zenith correction')
                axes[0,0].scatter(el_write, T0_write, s=5, c='b', label='Final set of Data Ponts')
                axes[0,0].legend()
                axes[0,0].set_title('Channel 0')

                axes[0,1].scatter(el_bl_corrected, T1_singslab, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[0,1].scatter(el_bl_corrected, T1_singslab_corr, s=1, c='c', alpha=0.5, label='After zenith correction')
                axes[0,1].scatter(el_write, T1_write, s=5, c='b', label='Final set of Data Ponts')
                axes[0,1].legend()
                axes[0,1].set_title('Channel 1')

                axes[1,0].scatter(el_bl_corrected, T2_singslab, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[1,0].scatter(el_bl_corrected, T2_singslab_corr, s=1, c='c', alpha=0.5, label='After zenith correction')
                axes[1,0].scatter(el_write, T2_write, s=5, c='b', label='Final set of Data Ponts')
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2')

                axes[1,1].scatter(el_bl_corrected, T3_singslab, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[1,1].scatter(el_bl_corrected, T3_singslab_corr, s=1, c='c', alpha=0.5, label='After zenith correction')
                axes[1,1].scatter(el_write, T3_write, s=5, c='b', label='Final set of Data Ponts')
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3')

                pl.suptitle('Single Slab Atm model data\nbefore and after Zenith correction Z_corr='+str(p_best[5])+' deg\n'+filename[:-4])
                pl.savefig('Output/Backlash_Fit/'+filename[:-4]+'_Zenith-corrected_data'+'_'+str(Npoints)+'.png')

                if showplots==1:
                    pl.show()
                else:
                    pl.close()

            dir=np.zeros(len(dir))

            if write_dat==1:

                for file in os.scandir('am_datfiles/'+template[:-4]+'/'+filename[:-4]):
                    if file.name.endswith("_red.dat"):
                        os.unlink(file.path)

                for i in range(len(el_write)):

                    with open('am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el_write[i],1))+'_am'+open_endstring,'w+') as f:
                        print('writing on'+' am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el_write[i],1))+'_am'+open_endstring)
                        f.write('1.25 1.50 {0}\n'.format(T0_write[i]))
                        f.write('3.25 2.50 {0}\n'.format(T1_write[i]))
                        f.write('5.5 2.00 {0}\n'.format(T2_write[i]))
                        f.write('7.25 1.50 {0}'.format(T3_write[i]))
                        print('writing done.')
                    f.close()

                pickle_Temps['am_dat_filename'].append('am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el_write[i],1))+'_'+dir_string[i]+'_am'+open_endstring)
                pickle_Temps['El'].append(round(el_write[i],1))
                pickle_Temps['dir'].append('NA')
                pickle_Temps['T0'].append(T0_write[i])
                pickle_Temps['T1'].append(T1_write[i])
                pickle_Temps['T2'].append(T2_write[i])
                pickle_Temps['T3'].append(T3_write[i])

            if os.path.exists(pickle_fn):
                os.remove(pickle_fn)
                print('Deleting old pickle file '+pickle_fn)


            f = open(pickle_fn,'wb')
            pickle.dump(pickle_Temps, f)
            f.close()

            f2 = open(pickle_betas,'wb')
            pickle.dump(betas, f2)
            f2.close()

        elif spline==2: #correct Temperature files for BL and Z_off but doesn't reduce the # points

            open_endstring='_corr.dat'
            with open(path_to_data+filename[:-4]+open_endstring,'w+') as f:
                time=data[:,0]
                T0=data[:,1]
                T1=data[:,2]
                T2=data[:,3]
                T3=data[:,4]
                el=data[:,5]


                time=time[el>0]
                T0=T0[el>0]
                T1=T1[el>0]
                T2=T2[el>0]
                T3=T3[el>0]
                el=el[el>0]

                T0[np.where(T0 == np.nan)]=0.
                T1[np.where(T0 == np.nan)]=0.
                T2[np.where(T0 == np.nan)]=0.
                T3[np.where(T0 == np.nan)]=0.

                previous_el=0
                dir=np.zeros(len(time))
                for j in range (len(time)):
                    if (el[j]>=previous_el):
                        dir[j]=1
                    elif (el[j]<previous_el):
                        dir[j]=-1
                    previous_el=el[j]

                up_index=np.where(dir==1)
                dn_index=np.where(dir==-1)

                T0_up=T0[np.where(dir==1)]
                T0_down=T0[np.where(dir==-1)]
                T1_up=T1[np.where(dir==1)]
                T1_down=T1[np.where(dir==-1)]
                T2_up=T2[np.where(dir==1)]
                T2_down=T2[np.where(dir==-1)]
                T3_up=T3[np.where(dir==1)]
                T3_down=T3[np.where(dir==-1)]

                el_up=el[np.where(dir==1)]
                el_down=el[np.where(dir==-1)]

                T0_coef_up=np.polyfit(el_up, T0_up, 7)
                T0_fit_up=np.poly1d(T0_coef_up)
                T1_coef_up=np.polyfit(el_up, T1_up, 7)
                T1_fit_up=np.poly1d(T1_coef_up)
                T2_coef_up=np.polyfit(el_up, T2_up, 7)
                T2_fit_up=np.poly1d(T2_coef_up)
                T3_coef_up=np.polyfit(el_up, T3_up, 7)
                T3_fit_up=np.poly1d(T3_coef_up)

                T0_coef_down=np.polyfit(el_down, T0_down, 7)
                T0_fit_down=np.poly1d(T0_coef_down)
                T1_coef_down=np.polyfit(el_down, T1_down, 7)
                T1_fit_down=np.poly1d(T1_coef_down)
                T2_coef_down=np.polyfit(el_down, T2_down, 7)
                T2_fit_down=np.poly1d(T2_coef_down)
                T3_coef_down=np.polyfit(el_down, T3_down, 7)
                T3_fit_down=np.poly1d(T3_coef_down)

                my_figsize=(12, 8)

                fig, axes = plt.subplots(2, 2, figsize=my_figsize)
                axes[0,0].scatter(el_up, T0_up, marker='^', s=3, alpha=0.5, c='red', label='Moving up')
                axes[0,0].plot(sorted(el_up), T0_fit_up(sorted(el_up)), c='orange', label='Fit up')
                axes[0,0].scatter(el_down, T0_down, marker='v', s=3, alpha=0.5, c='blue', label='Moving down')
                axes[0,0].plot(sorted(el_down), T0_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[0,0].legend(fontsize='large')
                axes[0,0].set_title('Channel 0', fontsize='large')
                axes[0,0].tick_params(axis='both', labelsize=14)
                axes[0,0].set_ylabel('T[K]', fontsize='large')

                axes[0,1].scatter(el_up, T1_up, marker='^', s=3, alpha=0.5, c='red', label='moving up')
                axes[0,1].plot(sorted(el_up), T1_fit_up(sorted(el_up)),  c='orange', label='Fit up')
                axes[0,1].scatter(el_down, T1_down, marker='v', s=3, alpha=0.5, c='blue', label='moving down')
                axes[0,1].plot(sorted(el_down), T1_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[0,1].legend(fontsize='large')
                axes[0,1].set_title('Channel 1', fontsize='large')
                axes[0,1].tick_params(axis='both', labelsize=14)
                axes[0,1].set_ylabel('T[K]', fontsize='large')

                axes[1,0].scatter(el_up, T2_up, marker='^', s=3, alpha=0.5, c='red', label='moving up')
                axes[1,0].plot(sorted(el_up), T2_fit_up(sorted(el_up)), c='orange', label='Fit up')
                axes[1,0].scatter(el_down, T2_down, marker='v', s=3, alpha=0.5, c='blue', label='moving down')
                axes[1,0].plot(sorted(el_down), T2_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[1,0].legend(fontsize='large')
                axes[1,0].set_title('Channel 2', fontsize='large')
                axes[1,0].tick_params(axis='both', labelsize=14)
                axes[1,0].set_ylabel('T[K]', fontsize='large')
                axes[1,0].set_xlabel('El[deg]', fontsize='large')

                axes[1,1].scatter(el_up, T3_up, marker='^', s=3, alpha=0.5, c='red', label='moving up')
                axes[1,1].plot(sorted(el_up), T3_fit_up(sorted(el_up)), c='orange', label='Fit up')
                axes[1,1].scatter(el_down, T3_down, marker='v', s=3, alpha=0.5, c='blue', label='moving down')
                axes[1,1].plot(sorted(el_down), T3_fit_down(sorted(el_down)), c='c', label='Fit down')
                axes[1,1].legend(fontsize='large')
                axes[1,1].set_title('Channel 3', fontsize='large')
                axes[1,1].tick_params(axis='both', labelsize=14)
                axes[1,1].set_ylabel('T[K]', fontsize='large')
                axes[1,1].set_xlabel('El[deg]', fontsize='large')

                fig.suptitle('T readings - '+filename[:-4]+'\nFitted with 7th order Polynomial', fontsize='large')
                pl.savefig('../Output/'+filename[:-4]+'_T_Readings.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #to fit the backlash
                el_smooth=np.linspace(15, 90, 200)

                def poly_diff(x, backlash, fit_up, fit_down):
                    return fit_up(x)-fit_up(x+backlash)
                def err(backlash, x, y, fit_up, fit_down):
                    return poly_diff(x, backlash, fit_up, fit_down) - y
                def err_global(p, x, y1, y2, y3, y4):
                    err1 = err(p, x, y1, T0_fit_up, T0_fit_down)
                    err2 = err(p, x, y2, T1_fit_up, T1_fit_down)
                    err3 = err(p, x, y3, T2_fit_up, T2_fit_down)
                    err4 = err(p, x, y4, T3_fit_up, T3_fit_down)

                    return np.concatenate((err1, err2, err3, err4))

                p=[3]

                ydata1=T0_fit_up(el_smooth)-T0_fit_down(el_smooth)
                ydata2=T1_fit_up(el_smooth)-T1_fit_down(el_smooth)
                ydata3=T2_fit_up(el_smooth)-T2_fit_down(el_smooth)
                ydata4=T3_fit_up(el_smooth)-T3_fit_down(el_smooth)

                p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global, p, args=(el_smooth, ydata1, ydata2, ydata3, ydata4), full_output=1)

                p_err=np.sqrt(np.diag(cov_x))

                fig, axes = plt.subplots(2, 2, figsize=my_figsize)
                axes[0, 0].plot(el_smooth, T0_fit_up(el_smooth)-T0_fit_down(el_smooth), label='diff_data')
                axes[0, 0].plot(el_smooth, poly_diff(el_smooth, *p_best, T0_fit_up, T0_fit_down), label='diff_fit')
                axes[0, 0].axhline(y=0, color='k', linestyle='--')
                axes[0, 0].set_title('Channel 0', fontsize='large')
                axes[0, 0].legend(loc='upper right', fontsize='large')
                axes[0, 0].set_ylabel('T[K]', fontsize='large')
                axes[0, 0].tick_params(axis='both', labelsize=14)

                axes[0, 1].plot(el_smooth, T1_fit_up(el_smooth)-T1_fit_down(el_smooth), label='diff_data')
                axes[0, 1].plot(el_smooth, poly_diff(el_smooth, *p_best, T1_fit_up, T1_fit_down), label='diff_fit')
                axes[0, 1].axhline(y=0, color='k', linestyle='--')
                axes[0, 1].set_title('Channel 1', fontsize='large')
                axes[0, 1].legend(loc='upper right', fontsize='large')
                axes[0, 1].set_ylabel('T[K]', fontsize='large')
                axes[0, 1].tick_params(axis='both', labelsize=14)

                axes[1, 0].plot(el_smooth, T2_fit_up(el_smooth)-T2_fit_down(el_smooth), label='diff_data')
                axes[1, 0].plot(el_smooth, poly_diff(el_smooth, *p_best, T2_fit_up, T2_fit_down), label='diff_fit')
                axes[1, 0].axhline(y=0, color='k', linestyle='--')
                axes[1, 0].set_title('Channel 2', fontsize='large')
                axes[1, 0].legend(loc='upper right', fontsize='large')
                axes[1, 0].set_ylabel('T[K]', fontsize='large')
                axes[1, 0].set_xlabel('El[deg]', fontsize='large')
                axes[1, 0].tick_params(axis='both', labelsize=14)

                axes[1, 1].plot(el_smooth, T3_fit_up(el_smooth)-T3_fit_down(el_smooth), label='diff_data')
                axes[1, 1].plot(el_smooth, poly_diff(el_smooth, *p_best, T3_fit_up, T3_fit_down), label='diff_fit')
                axes[1, 1].axhline(y=0, color='k', linestyle='--')
                axes[1, 1].set_title('Channel 3', fontsize='large')
                axes[1, 1].legend(loc='upper right', fontsize='large')
                axes[1, 1].set_ylabel('T[K]', fontsize='large')
                axes[1, 1].set_xlabel('El[deg]', fontsize='large')
                axes[1, 1].tick_params(axis='both', labelsize=14)

                fig.suptitle(filename[:-4]+' Backlash Fit \n Backlash[deg]='+str(round(p_best[0],3))+'±'+str(round(p_err[0],3)), fontsize='large')
                pl.savefig('../Output/'+filename[:-4]+'_line_diff_fit.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #Backlash correction

                el_down_prime=el_down+p_best

                fig, axes = plt.subplots(2, 2, figsize=my_figsize)
                axes[0,0].scatter(el_up, T0_up, marker='^', s=3, alpha=0.5, c='red', label='Moving up')
                axes[0,0].scatter(el_down_prime, T0_down, marker='v', s=3, alpha=0.5, c='blue', label='Moving down - BL corrected')
                axes[0,0].legend()
                axes[0,0].set_title('Channel 0', fontsize='large')
                axes[0,0].tick_params(axis='both', labelsize=14)
                axes[0,0].set_ylabel('T[K]', fontsize='large')

                axes[0,1].scatter(el_up, T1_up, marker='^', s=3, alpha=0.5, c='red', label='Moving up')
                axes[0,1].scatter(el_down_prime, T1_down, marker='v', s=3, alpha=0.5, c='blue', label='Moving down - BL corrected')
                axes[0,1].legend()
                axes[0,1].set_title('Channel 1', fontsize='large')
                axes[0,1].tick_params(axis='both', labelsize=14)
                axes[0,1].set_ylabel('T[K]', fontsize='large')

                axes[1,0].scatter(el_up, T2_up, marker='^', s=3, alpha=0.5, c='red', label='Moving up')
                axes[1,0].scatter(el_down_prime, T2_down, marker='v', s=3, alpha=0.5, c='blue', label='Moving down - BL corrected')
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2', fontsize='large')
                axes[1,0].tick_params(axis='both', labelsize=14)
                axes[1,0].set_ylabel('T[K]', fontsize='large')
                axes[1,0].set_xlabel('El[deg]', fontsize='large')

                axes[1,1].scatter(el_up, T3_up, marker='^', s=3, alpha=0.5, c='red', label='Moving up')
                axes[1,1].scatter(el_down_prime, T3_down, marker='v', s=3, alpha=0.5, c='blue', label='Moving down - BL corrected')
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3', fontsize='large')
                axes[1,1].tick_params(axis='both', labelsize=14)
                axes[1,1].set_ylabel('T[K]', fontsize='large')
                axes[1,1].set_xlabel('El[deg]', fontsize='large')

                pl.suptitle('Backlash-corrected data - ' + filename[:-4]+'\nBL correction[deg] = '+ str(p_best[0]), fontsize='large')
                pl.savefig('../Output/'+filename[:-4]+'_Backlash-corrected_data.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #Creating a single array for rad moving up&down after correcting Backlash
                dir_up=dir[np.where(dir==1)]
                dir_down=dir[np.where(dir==-1)]
                # dir_bl_corrected=np.concatenate((dir_up, dir_down), axis=None)
                #
                # el_bl_corrected=np.concatenate((el_up, el_down_prime), axis=None)
                # T0_bl_corrected=np.concatenate((T0_up, T0_down), axis=None)
                # T1_bl_corrected=np.concatenate((T1_up, T1_down), axis=None)
                # T2_bl_corrected=np.concatenate((T2_up, T2_down), axis=None)
                # T3_bl_corrected=np.concatenate((T3_up, T3_down), axis=None)

                dir_bl_corrected=np.zeros(len(dir))
                el_bl_corrected=np.zeros(len(dir))
                T0_bl_corrected=np.zeros(len(dir))
                T1_bl_corrected=np.zeros(len(dir))
                T2_bl_corrected=np.zeros(len(dir))
                T3_bl_corrected=np.zeros(len(dir))

                dir_bl_corrected[up_index]=dir_up
                dir_bl_corrected[dn_index]=dir_down

                el_bl_corrected[up_index]=el_up
                el_bl_corrected[dn_index]=el_down_prime

                T0_bl_corrected[up_index]=T0_up
                T0_bl_corrected[dn_index]=T0_down

                T1_bl_corrected[up_index]=T1_up
                T1_bl_corrected[dn_index]=T1_down

                T2_bl_corrected[up_index]=T2_up
                T2_bl_corrected[dn_index]=T2_down

                T3_bl_corrected[up_index]=T3_up
                T3_bl_corrected[dn_index]=T3_down


                zenith_corr=np.zeros(4)

                #Finding the Zenith
                def line_shape_1lay(x, p):
                    tau=np.exp(-p[1]/np.sin((x+p[2])*np.pi/180.))
                    T_meas=p[0]*(1.-tau)
                    return T_meas

                def err_lineshape(p, x, y):
                    return line_shape_1lay(x, p) - y

                def err_global_lineshape(p, x, y1, y2, y3, y4):
                    p1=p[0], p[1], p[5]
                    p2=p[0], p[2], p[5]
                    p3=p[0], p[3], p[5]
                    p4=p[0], p[4], p[5]

                    err1 = err_lineshape(p1, x, y1)
                    err2 = err_lineshape(p2, x, y2)
                    err3 = err_lineshape(p3, x, y3)
                    err4 = err_lineshape(p4, x, y4)

                    return np.concatenate((err1, err2, err3, err4))

                p=[200, 0.6, 0.25, 0.12, 0.05, 0.]

                p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global_lineshape, p, args=(el_bl_corrected, T0_bl_corrected, T1_bl_corrected, T2_bl_corrected, T3_bl_corrected), full_output=1)

                p_err=np.sqrt(np.diag(cov_x))
                z_corr_err=p_err[5]

                p1=p_best[0], p_best[1], p_best[5]
                p2=p_best[0], p_best[2], p_best[5]
                p3=p_best[0], p_best[3], p_best[5]
                p4=p_best[0], p_best[4], p_best[5]

                fig, axes = plt.subplots(2, 2, figsize=my_figsize)
                model_string='\nAtmospheric Model:\n T[ch]=T_atm(1-tau[ch])\ntau[ch]=exp(-beta[ch]/sin(el+z_corr))'
                string0='Output parameters:\n\nT_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nbeta={:.2f}'.format(round(p_best[1],2))+'±{:.2f}'.format(round(p_err[1],2))+'\nz_corr={:.2f}'.format(round(p_best[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                pos0y=T0_bl_corrected.max()-(T0_bl_corrected.max()-T0_bl_corrected.min())/3.
                axes[0,0].text(60.,  pos0y-(T0_bl_corrected.max()/10.), string0, fontsize=18)
                axes[0,0].scatter(el_bl_corrected, T0_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[0,0].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p1), s=2, c='r', label='Fit - Single slab atmospheric model')
                axes[0,0].legend()
                axes[0,0].set_title('Channel 0')

                string1='Output parameters:\n\nT_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nbeta={:.2f}'.format(round(p_best[2],2))+'±{:.2f}'.format(round(p_err[2],2))+'\nz_corr={:.2f}'.format(round(p_best[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                pos1y=T1_bl_corrected.max()-(T1_bl_corrected.max()-T1_bl_corrected.min())/3.
                axes[0,1].text(60.,  pos1y-(T1_bl_corrected.max()/10.), string1, fontsize=18)
                axes[0,1].scatter(el_bl_corrected, T1_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[0,1].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p2), s=2, c='r', label='Fit - Single slab atmospheric model')
                axes[0,1].legend()
                axes[0,1].set_title('Channel 1')

                string2='Output parameters:\n\nT_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nbeta={:.2f}'.format(round(p_best[3],2))+'±{:.2f}'.format(round(p_err[3],2))+'\nz_corr={:.2f}'.format(round(p_best[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                pos2y=T2_bl_corrected.max()-(T2_bl_corrected.max()-T2_bl_corrected.min())/3.
                axes[1,0].text(60.,  pos2y-(T2_bl_corrected.max()/10.), string2, fontsize=18)
                axes[1,0].scatter(el_bl_corrected, T2_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[1,0].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p3), s=2, c='r', label='Fit - Single slab atmospheric model')
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2')

                string3='Output parameters:\n\nT_atm ={:.2f}'.format(round(p_best[0],2))+'±{:.2f}'.format(round(p_err[0],2))+'\nbeta={:.2f}'.format(round(p_best[4],2))+'±{:.2f}'.format(round(p_err[4],2))+'\nz_corr={:.2f}'.format(round(p_best[5],2))+'±{:.2f}'.format(round(p_err[5],2))
                pos3y=T3_bl_corrected.max()-(T3_bl_corrected.max()-T3_bl_corrected.min())/3.
                axes[1,1].text(60.,  pos3y-(T3_bl_corrected.max()/10.), string3, fontsize=18)
                axes[1,1].scatter(el_bl_corrected, T3_bl_corrected, s=3, c='k', label='T data after BL correction')
                axes[1,1].scatter(el_bl_corrected, line_shape_1lay(el_bl_corrected, p4), s=2, c='r', label='Fit - Single slab atmospheric model')
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3')

                pl.suptitle('Fit to one slab atmospheric model - on backlash-corrected data\n' + filename[:-4]+'\n'+model_string)
                pl.savefig('../Output/'+filename[:-4]+'_SingleSlabModelFit.png')

                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                #To plot the tilt var
                fig, axes = plt.subplots(2, 2)

                el_tilt=np.linspace(30., 80., 500)
                T0_tilt=line_shape_1lay(el_tilt, p1)
                T1_tilt=line_shape_1lay(el_tilt, p2)
                T2_tilt=line_shape_1lay(el_tilt, p3)
                T3_tilt=line_shape_1lay(el_tilt, p4)



                delta=np.zeros(4) #DT/D(el) for each of the 4 channels around el=55deg

                el55_list=np.where(el_tilt>=55.)
                el55_list=el55_list[0]

                center0=T0_tilt[el55_list[0]]
                up0=center0+(0.01*center0)
                down0=center0-(0.01*center0)
                ind0_up=np.where(T0_tilt<=up0)
                ind0_up=ind0_up[0]
                if down0>=T0_tilt[len(T0_tilt)-1]:
                    ind0_dn=np.where(T0_tilt<=down0)
                    ind0_dn=ind0_dn[0]
                    i0_dn=ind0_dn[0]
                else:
                    i0_dn=len(T0_tilt)-1

                T0_tilt=T0_tilt[ind0_up[0]:i0_dn] #ind_up is at lower el
                el_tilt0=el_tilt[ind0_up[0]:i0_dn]
                delta[0]=(np.max(T0_tilt)-np.min(T0_tilt))/(np.max(el_tilt0)-np.min(el_tilt0))

                center1=T1_tilt[el55_list[0]]
                up1=center1+(0.01*center1)
                down1=center1-(0.01*center1)
                ind1_up=np.where(T1_tilt<=up1)
                ind1_up=ind1_up[0]
                ind1_dn=np.where(T1_tilt<=down1)
                ind1_dn=ind1_dn[0]
                T1_tilt=T1_tilt[ind1_up[0]:ind1_dn[0]] #ind_up is at lower el
                el_tilt1=el_tilt[ind1_up[0]:ind1_dn[0]]
                delta[1]=(np.max(T1_tilt)-np.min(T1_tilt))/(np.max(el_tilt1)-np.min(el_tilt1))

                center2=T2_tilt[el55_list[0]]
                up2=center2+(0.01*center2)
                down2=center2-(0.01*center2)
                ind2_up=np.where(T2_tilt<=up2)
                ind2_up=ind2_up[0]
                ind2_dn=np.where(T2_tilt<=down2)
                ind2_dn=ind2_dn[0]
                T2_tilt=T2_tilt[ind2_up[0]:ind2_dn[0]] #ind_up is at lower el
                el_tilt2=el_tilt[ind2_up[0]:ind2_dn[0]]
                delta[2]=(np.max(T2_tilt)-np.min(T2_tilt))/(np.max(el_tilt2)-np.min(el_tilt2))

                center3=T3_tilt[el55_list[0]]
                up3=center3+(0.01*center3)
                down3=center3-(0.01*center3)
                ind3_up=np.where(T3_tilt<=up3)
                ind3_up=ind3_up[0]
                ind3_dn=np.where(T3_tilt<=down3)
                ind3_dn=ind3_dn[0]
                T3_tilt=T3_tilt[ind3_up[0]:ind3_dn[0]] #ind_up is at lower el
                el_tilt3=el_tilt[ind3_up[0]:ind3_dn[0]]
                delta[3]=(np.max(T3_tilt)-np.min(T3_tilt))/(np.max(el_tilt3)-np.min(el_tilt3))

                axes[0,0].scatter(el_tilt0, T0_tilt, s=3, c='r', label='ΔT/Δ(el)[K/deg]='+str(round(delta[0],2)))
                #axes[0,0].set_ylim(center0-(0.1*center0), center0+(0.1*center0))
                axes[0,0].set_title('Channel 0')
                leg=axes[0,0].legend(loc='upper right')
                for item in leg.legendHandles:
                    item.set_visible(False)

                axes[0,1].scatter(el_tilt1, T1_tilt, s=3, c='r', label='ΔT/Δ(el)[K/deg]='+str(round(delta[1],2)))
                axes[0,1].set_title('Channel 1')
                leg=axes[0,1].legend(loc='upper right')
                for item in leg.legendHandles:
                    item.set_visible(False)

                axes[1,0].scatter(el_tilt2, T2_tilt, s=3, c='r', label='ΔT/Δ(el)[K/deg]='+str(round(delta[2],2)))
                axes[1,0].set_title('Channel 2')
                leg=axes[1,0].legend(loc='upper right')
                for item in leg.legendHandles:
                    item.set_visible(False)

                axes[1,1].scatter(el_tilt3, T3_tilt, s=3, c='r', label='ΔT/Δ(el)[K/deg]='+str(round(delta[3],2)))
                axes[1,1].set_title('Channel 3')
                leg=axes[1,1].legend(loc='upper right')
                for item in leg.legendHandles:
                    item.set_visible(False)

                pl.suptitle('T vs El Angle variation for WVR Tilt')
                pl.savefig('../Output/'+filename[:-4]+'_dTdEl.png')

                pl.close()


                betas=[p_best[1], p_best[2], p_best[3], p_best[4]]

                el_bl_z_corr = el_bl_corrected - p_best[5]

                fig, axes = plt.subplots(2, 2, figsize=(18, 10))
                #axes[0,0].scatter(el_bl_corrected, T0_bl_corrected, s=1, c='k', alpha=0.5, label='Before Zenith Corr')
                axes[0,0].scatter(el_bl_corrected, T0_bl_corrected, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[0,0].scatter(el_bl_z_corr, T0_bl_corrected, s=1, c='c', alpha=0.5, label='After zenith correction')
                #axes[0,0].axvline(x=zenith_ch[0], color='r', linestyle='--', label='Ch Zenith[deg]={:.2f}'.format(zenith_ch[0]))
                axes[0,0].legend()
                axes[0,0].set_title('Channel 0')

                #axes[0,1].scatter(el_bl_corrected, T1_bl_corrected, s=1, c='k', alpha=0.5, label='Before Zenith Corr')
                axes[0,1].scatter(el_bl_corrected, T1_bl_corrected, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[0,1].scatter(el_bl_z_corr, T1_bl_corrected, s=1, c='c', alpha=0.5, label='After zenith correction')
                #axes[0,1].axvline(x=zenith_ch[1], color='r', linestyle='--', label='Ch Zenith[deg]={:.2f}'.format(zenith_ch[1]))
                axes[0,1].legend()
                axes[0,1].set_title('Channel 1')

                #axes[1,0].scatter(el_bl_corrected, T2_bl_corrected, s=1, c='k', alpha=0.5, label='Before Zenith Corr')
                axes[1,0].scatter(el_bl_corrected, T2_bl_corrected, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[1,0].scatter(el_bl_z_corr, T2_bl_corrected, s=1, c='c', alpha=0.5, label='After zenith correction')
                #axes[1,0].axvline(x=zenith_ch[2], color='r', linestyle='--', label='Ch Zenith[deg]={:.2f}'.format(zenith_ch[2]))
                axes[1,0].legend()
                axes[1,0].set_title('Channel 2')

                #axes[1,1].scatter(el_bl_corrected, T3_bl_corrected, s=1, c='k', alpha=0.5, label='Before Zenith Corr')
                axes[1,1].scatter(el_bl_corrected, T3_bl_corrected, s=1, c='k', alpha=0.5, label='Before zenith correction')
                axes[1,1].scatter(el_bl_z_corr, T3_bl_corrected, s=1, c='c', alpha=0.5, label='After zenith correction')
                #axes[1,1].axvline(x=zenith_ch[3], color='r', linestyle='--', label='Ch Zenith[deg]={:.2f}'.format(zenith_ch[3]))
                axes[1,1].legend()
                axes[1,1].set_title('Channel 3')

                pl.suptitle('Single Slab Atm model data\nbefore and after Zenith correction Z_corr='+str(round(p_best[5],2))+' deg\n'+filename[:-4])
                pl.savefig('../Output/'+filename[:-4]+'_Zenith-corrected_data.png')

                if showplots==1:
                    pl.show()
                else:
                    pl.close()

            el_write=el_bl_z_corr
            T0_write=T0_bl_corrected
            T1_write=T1_bl_corrected
            T2_write=T2_bl_corrected
            T3_write=T3_bl_corrected
            dir_write=dir_bl_corrected

            print('T0_write=', T0_write)

            if write_dat==1:

                for file in os.scandir('am_datfiles/'+template[:-4]+'/'+filename[:-4]):
                    if file.name.endswith("_red.dat"):
                        os.unlink(file.path)

                dir_string=['None']*len(dir_write)

                for i in range(len(el_write)):
                    if dir_write[i]==1:
                        dir_string[i]='up'
                    elif dir_write[i]==-1:
                        dir_string[i]='dn'
                    with open('am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el_write[i],1))+'_'+dir_string[i]+'_am'+open_endstring,'w+') as f:
                        print('writing on'+'am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el_write[i],1))+'_'+dir_string[i]+'_am'+open_endstring)
                        f.write('1.25 1.50 {0}\n'.format(T0_write[i]))
                        f.write('3.25 2.50 {0}\n'.format(T1_write[i]))
                        f.write('5.5 2.00 {0}\n'.format(T2_write[i]))
                        f.write('7.25 1.50 {0}'.format(T3_write[i]))
                        print('writing done.')
                    f.close()

                    pickle_Temps['am_dat_filename'].append('am_datfiles/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_El'+str(round(el_write[i],1))+'_'+dir_string[i]+'_am'+open_endstring)
                    pickle_Temps['El'].append(round(el_write[i],1))
                    pickle_Temps['dir'].append(dir_string[i])
                    pickle_Temps['T0'].append(T0_write[i])
                    pickle_Temps['T1'].append(T1_write[i])
                    pickle_Temps['T2'].append(T2_write[i])
                    pickle_Temps['T3'].append(T3_write[i])



            if os.path.exists(pickle_fn):
                os.remove(pickle_fn)
                print('Deleting old pickle file '+pickle_fn)

            print('pickle_Temps=', pickle_Temps)
            f = open(pickle_fn,'wb')
            pickle.dump(pickle_Temps, f)
            f.close()

            f2 = open(pickle_betas,'wb')
            pickle.dump(betas, f2)
            f2.close()

            #zenith_eff=90.-p_best[5]
            #print('zenith_eff=', zenith_eff)

            z_offs_elnod['z_offs']=p_best[5]
            z_offs_elnod['z_offs_err']=z_corr_err

            f3=open(pickle_z, 'wb')
            pickle.dump(z_offs_elnod, f3)
            f3.close()

            f4=open(pickle_deltas, 'wb')
            pickle.dump(delta, f4)
            f4.close()

            return p_best[5], z_corr_err, delta

    def create_am_datfile_Az(self,filename, pathtofn='wvr1_data/', clean_mod=3, clean_method='import_model', path_to_data='', template='SPole_annual_50.amc'):
        if not os.path.exists(self.path_to_all+'am_datfiles_Az/'+template[:-4]):
            os.makedirs(self.path_to_all+'am_datfiles_Az/'+template[:-4])
        if not os.path.exists(self.path_to_all+'am_datfiles_Az/'+template[:-4]+'/'+filename[:-4]):
            os.makedirs(self.path_to_all+'am_datfiles_Az/'+template[:-4]+'/'+filename[:-4])

        if clean_mod==0:
            D_clean, waz, mod_removed_data=raz.read_Az_fast(filename, pathtofn=pathtofn, clean_mod=clean_mod, clean_method=clean_method)
        else:
            D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single, calib_data_Gavg, FH =raz.read_Az_fast(filename, pathtofn=pathtofn, clean_mod=clean_mod, clean_method=clean_method)

        path=self.path_to_all+'am_datfiles_Az/'+template[:-4]+'/'+filename[:-4]
        pickle_fn=path+'/'+filename[:-4]+'_clean_mod'+str(clean_mod)+'_method_'+clean_method+'_pickle_temps.txt'

        pickle_Temps= {'am_dat_filename': [], 'wAz':[],'scanN':[], 'T0':[], 'T1':[], 'T2':[], 'T3':[]}
        open_endstring='.dat'
        wwaz,fs =  raz.findScans(waz)
        nscans = len(fs.s)

        print('wwaz=', wwaz)
        print(nscans)

        T0=mod_removed_data[:,0]
        T1=mod_removed_data[:,1]
        T2=mod_removed_data[:,2]
        T3=mod_removed_data[:,3]


        pickle_Temps['T0']=T0
        pickle_Temps['T1']=T1
        pickle_Temps['T2']=T2
        pickle_Temps['T3']=T3

        pickle_Temps['wAz']=waz

        for azn in range(len(waz)):
            (scanN, az) = divmod(waz[azn],360)
            scanN=int(scanN)
            az=round(az,2)
            print('azn=', azn)
            print('scanN=',scanN)
            print('az=', az)
            dat_fn=self.path_to_all+'am_datfiles_Az/'+template[:-4]+'/'+filename[:-4]+'/'+filename[:-4]+'_scanN'+str(scanN)+'_Az'+str(az)+'_clean_mod'+str(clean_mod)+'_method_'+clean_method+'_am'+open_endstring
            with open(dat_fn,'w+') as f:
                print('writing on'+dat_fn)
                f.write('1.25 1.50 {0}\n'.format(T0[azn]))
                f.write('3.25 2.50 {0}\n'.format(T1[azn]))
                f.write('5.5 2.00 {0}\n'.format(T2[azn]))
                f.write('7.25 1.50 {0}'.format(T3[azn]))
                print('writing done.')
                f.close()

                pickle_Temps['am_dat_filename'].append(dat_fn)
                pickle_Temps['scanN'].append(scanN)

        if os.path.exists(pickle_fn):
            os.remove(pickle_fn)
            print('Deleting old pickle file '+pickle_fn)

        f = open(pickle_fn,'wb')
        pickle.dump(pickle_Temps, f)
        f.close()

        return pickle_Temps


    def create_template(self, date, time_h, starting_temp='SPole_annual_50.amc', path_to_temp='Templates/'):

        YY=date[:4]
        mm=date[5:7]
        dd=date[8:10]

        d=datetime.datetime(int(YY), int(mm), int(dd), int(time_h), 0, 0)

        f = open('SP_MERRA-2.pickle','rb')
        p = pickle.load(f)
        f.close()

        t_unix=p['t']
        t_iso = [datetime.datetime.fromtimestamp(t).isoformat() for t in t_unix]

        t_datetime = [datetime.datetime.fromisoformat(t) for t in t_iso]
        t_date = [t.date() for t in t_datetime]
        t_time = [t.time() for t in t_datetime]

        d_index=np.where(np.array(t_datetime)==d)

        h=p['h'] #m
        T=np.array(p['T'][:,d_index]).reshape(np.shape(h)) #K
        P=np.array(p['PL'][:,d_index]).reshape(np.shape(h)) #Pa
        mmr=np.array(p['QV'][:,d_index]).reshape(np.shape(h)) #mass mixing ratio [g/cm^3]?
        # pl.plot(h, T)
        # pl.xlabel('h[m]')
        # pl.ylabel('T[K]')
        # pl.close()
        # pl.plot(h, P)
        # pl.xlabel('h[m]')
        # pl.ylabel('P[Pa]')
        # pl.close()
        # pl.plot(h, mmr)
        # pl.xlabel('h[m]')
        # pl.ylabel('mmr')
        # pl.close()

        file_variable = open(path_to_temp+starting_temp)
        all_lines_variable = file_variable.readlines()
        out_lines=[]
        #for i in range (50):
        for i in range (175):
            out_lines.append(all_lines_variable[i])

        out_temp = open(path_to_temp+'MERRA_'+date+'_'+time_h+'.amc', "w")

        out_lines.append('')
        out_lines.append('\n')
        h_top=100
        for j in range (len(h)):
            h_top+=100
            i=len(h)-j-1
            vmr_h2o=(28.9644/18.01528)*mmr[i]
            out_lines.append('layer troposphere\n')
            #out_lines.append('h 100.0 m\n') #am requires layer thickness -- Layers will be stacked on top of each other
            out_lines.append('Pbase '+str(P[i])+' Pa\n')
            out_lines.append('Tbase '+str(T[i])+' K\n')
            out_lines.append('column dry_air vmr\n')
            out_lines.append('column h2o vmr '+str(vmr_h2o)+'\n')
            out_lines.append('')
            out_lines.append('\n')
        print('h_top=', h_top)
        #
        # while (h_top <= 70000.): #adding other dry layers
        #     out_lines.append('layer\n')
        #     out_lines.append('h 100.0 m\n')
        #     out_lines.append('T 240 K\n')
        #     out_lines.append('P 10 Pa\n')
        #     out_lines.append('column dry_air vmr\n')
        #     out_lines.append('column h2o vmr 4.76e-06\n')
        #     out_lines.append('column o3 vmr 6.42e-07\n')
        #     out_lines.append('')
        #     out_lines.append('\n')
        #     h_top+=100

        print('h_top=', h_top)

        out_temp.writelines(out_lines)
        out_temp.close()

    def fit_w_am(self, filename, path_to_data='wvr1_data/', template= 'SPole_annual_50.amc', path_to_temp='Templates/', el_range=(0,90.), spline=2):
        Npoints=0

        if not os.path.exists('am_datfiles/'+template[:-4]):
            os.makedirs('am_datfiles/'+template[:-4])

        path='am_datfiles/'+template[:-4]+'/'+filename[:-4]

        if not os.path.exists(path):
            os.makedirs(path)

        if spline==0:
            pickle_fn=path+'/'+filename[:-4]+'_fitoutput.txt'
            pickle_fn_read=path+'/'+filename[:-4]+'_pickle_temps.txt'
        elif spline==1:
            pickle_fn=path+'/'+filename[:-4]+'_fitoutput_red.txt'
            pickle_fn_read=path+'/'+filename[:-4]+'_pickle_temps_red.txt'
        elif spline==2:
            pickle_fn=path+'/'+filename[:-4]+'_fitoutput_corr.txt'
            pickle_fn_read=path+'/'+filename[:-4]+'_pickle_temps_corr.txt'
        low_el, high_el = el_range
        fit_output = {'filename': [], 'El': [], 'dir':[], 'Nscale': [], 'Nscale_err': [],
        'pwv_meso':[], 'pwv_strato':[], 'pwv_tropo':[], 'pwv_total':[], 'pwv_los_total':[]}
        #try:
        self.create_am_datfile(filename, path_to_data=path_to_data, template=template, spline=spline, showplots=0)
        print('dat file created')
        if spline==1:
            end_string='ed.dat'
            end_string_fn='_red'

        elif spline==0:
            end_string='am.dat'
            end_string_fn='_full'

        elif spline==2:
            end_string='rr.dat'
            end_string_fn='_corr'


        t = open(pickle_fn_read, 'rb')
        pickle_Temps = pickle.load(t)
        t.close()

        am_dat_list=pickle_Temps['am_dat_filename']
        el_list=pickle_Temps['El']
        dir_list=pickle_Temps['dir']


        #for flist in os.listdir(path):
        #    if flist[-6:]==end_string:
        T0_check=[]
        for findex in range (len(am_dat_list)):
            flist=am_dat_list[findex]
            print(flist)
            #print('flist[-6:]=',flist[-6:])
            #print('end_string=', end_string)
            if flist[-6:]==end_string:
                path_to_dat=flist
                #print(flist[-6:])
                el= el_list[findex]
                zen_ang=90-el
                print('zen_ang=', zen_ang)
                #print('El from list:', el)
                dir=flist[35:37]
                dir=dir_list[findex]
                if (low_el<= el <=high_el):
                    try:
                        am_template=path_to_temp+template
                        if not os.path.exists(path_to_dat+'.amc'):
                            os.system(f"am {am_template}  {zen_ang}  {path_to_dat}")
                        print(path_to_dat)

                        dat_file_variable = open(path_to_dat)
                        dat_lines = dat_file_variable.readlines()
                        T0_line=dat_lines[0]
                        T0_check.append(T0_line[10:15])
                        print('dat_file:', dat_lines)
                        file_variable = open(path_to_dat+'.amc')
                        all_lines_variable = file_variable.readlines()
                        nscale_line=all_lines_variable[44]
                        print('nscale_line=', nscale_line)
                        std_line=all_lines_variable[10]
                        print('std_line=', std_line)
                        if template[0:5]=='MERRA':
                            #pwv_line=all_lines_variable[1062] #wo dry layers
                            #pwv_line=all_lines_variable[1632]#[1633]
                            pwv_line_meso=all_lines_variable[1614]
                            pwv_line_strato=all_lines_variable[1621]#[1633]
                            pwv_line_tropo=all_lines_variable[1628]#[1633]
                            pwv_line_total=all_lines_variable[1633]#[1633]
                        else:
                            #pwv_line=all_lines_variable[503]#[504]
                            pwv_line_meso=all_lines_variable[483]
                            pwv_line_strato=all_lines_variable[490]#[1633]
                            pwv_line_tropo=all_lines_variable[497]#[1633]
                            pwv_line_total=all_lines_variable[504]#[1633]
                        print('pwv_line_tropo=', pwv_line_tropo)
                        Nscale=float(nscale_line[23:30])
                        print('nscale=', Nscale)
                        Nscale_std=float(std_line[34:38])
                        print('nscale_std=', Nscale_std)
                        pwv_meso=float(pwv_line_meso[30:36])
                        pwv_strato=float(pwv_line_strato[30:36])
                        pwv_tropo=float(pwv_line_tropo[30:36])
                        pwv_total=float(pwv_line_total[30:36])
                        print('pwv_total=', pwv_total)
                        pwv_los_total=float(pwv_line_total[49:56])
                        print('pwv_los=', pwv_los_total)
                        #print('T0_check=', T0_line[10:15])
                        fit_output['El'].append(el)
                        fit_output['dir'].append(dir)
                        fit_output['filename'].append(flist)
                        fit_output['Nscale'].append(Nscale)
                        fit_output['Nscale_err'].append(Nscale_std)
                        fit_output['pwv_meso'].append(pwv_meso)
                        fit_output['pwv_strato'].append(pwv_strato)
                        fit_output['pwv_tropo'].append(pwv_tropo)
                        fit_output['pwv_total'].append(pwv_total)
                        fit_output['pwv_los_total'].append(pwv_los_total)
                        Npoints+=1
                    except Exception as e:
                        print(f'{flist} failed: {e}')
        # except:
        #     zenith=np.nan
        #     print(filename+' failed.')


        if os.path.exists(pickle_fn):
            os.remove(pickle_fn)
            print('Deleting old pickle file '+pickle_fn)

        print('Writing on Pickle File '+pickle_fn)
        #print('Fit_output=', fit_output)

        f = open(pickle_fn,'wb')
        pickle.dump(fit_output, f)
        f.close()

        return(fit_output)


    def fit_w_am_Az(self, filename, clean_mod=3, clean_method='import_model', out_path='output_plots/', path_to_temp='Templates/', template='SPole_annual_50.amc'):

        #clean_mod=1--> removes single mod
        #clean_mod=2--> removes double mod
        #clean_mod=3--> removes single and double mod

        #clean_mod=1 (or 3) has to be used as an alternative to input the tilt correction par

        El=55.0

        if not os.path.exists(self.path_to_all+'am_datfiles_Az/'+template[:-4]):
            os.makedirs(self.path_to_all+'am_datfiles_Az/'+template[:-4])

        path=self.path_to_all+'am_datfiles_Az/'+template[:-4]+'/'+filename[:-4]
        #pathtofn='wvr1_data'+filename[:-9]+'/'
        pickle_fn=path+'/'+filename[:-4]+'_clean_mod'+str(clean_mod)+'_method_'+str(clean_method)+'_fitoutput.txt' #just double mod removed+tilt_correction
        pickle_fn_temps=path+'/'+filename[:-4]+'_clean_mod'+str(clean_mod)+'_method_'+str(clean_method)+'_pickle_temps.txt'
        fit_output = {'filename': [], 'Az': [], 'scanN':[], 'Nscale': [], 'Nscale_err': [], 'pwv':[], 'pwv_err':[], 'el_correction':[]}

        if not os.path.exists(pickle_fn_temps):
            pickle_Temps = self.create_am_datfile_Az(filename, pathtofn='wvr1_data/', clean_mod=clean_mod, clean_method=clean_method)
        else:
            print(pickle_fn_temps+' already exists.\n')
            answer = input("Do you want to overwrite it? ")
            if answer == "y":
                pickle_Temps = self.create_am_datfile_Az(filename, pathtofn='wvr1_data/', clean_mod=clean_mod, clean_method=clean_method)
            elif answer == "n":
                print('Loading pickle temps.')
                f = open(pickle_fn_temps,'rb')
                pickle_Temps = pickle.load(f)
                f.close()
            else:
                print("Please enter y or n.")


        end_string='am.dat'
        am_dat_list=pickle_Temps['am_dat_filename']
        waz_list=pickle_Temps['wAz']
        N_list=pickle_Temps['scanN']

        waz=pickle_Temps['wAz']
        scanN=pickle_Temps['scanN']

        pwv_ts=np.full(len(waz), np.nan)
        pwv_ts_mod3=np.full(len(waz), np.nan)
        nscale_ts=np.full(len(waz), np.nan)
        nscale_err_ts=np.full(len(waz), np.nan)

        D_pwv=np.full(np.shape(pickle_Temps['T0']), np.nan)

        el_correction=np.zeros(len(waz))

        az,fs =  raz.findScans(waz)

        pl.plot(waz[0:int(len(waz)/10)], label='waz')
        pl.plot(az[0:int(len(waz)/10)], label='az')
        pl.title('waz vs az')
        pl.legend()
        pl.close()

        am_template=path_to_temp+template

        for azn in range(len(waz)):
            print('\n\n\n\n\niteration=', azn)
            print('az=', az)
            print('pwv_ts=', pwv_ts)
            print('fs=', fs)
            D_pwv = raz.interpToImage(az, pwv_ts, fs)
            (scanN, azN) = divmod(waz[azn],360)
            scanN=int(scanN)
            print('az[azn]=', az[azn])
            print('azN=', azN)
            azN=round(azN,2)
            print('azN_rounded=', azN)

            fig=plt.figure()
            im=plt.imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
            plt.suptitle('PWV Atmogram tmp - '+filename[:-4]+'\naz='+str(azN)+' - scanN='+str(scanN)+'\n iteration = '+str(azn))
            cbar = fig.colorbar(im, extend='both')
            cbar.set_label('PWV[um]')
            plt.savefig(out_path+filename[:-4]+'_pwvatmo'+'_clean_mod'+str(clean_mod)+'_clean_method_'+str(clean_method)+'_tmp.png')
            plt.close()
            path_to_dat=path+'/'+filename[:-4]+'_scanN'+str(scanN)+'_Az'+str(azN)+'_clean_mod'+str(clean_mod)+'_method_'+str(clean_method)+'_am.dat'
            print('path_to_dat=', path_to_dat)
            fit_output['Az'].append(az[azn])
            fit_output['scanN'].append(scanN)
            fit_output['filename'].append(filename[:-4]+'_scanN'+str(scanN)+'_Az'+str(azN)+'_clean_mod'+str(clean_mod)+'_method_'+str(clean_method)+'_am.dat')


            try:
                if not os.path.exists(path_to_dat+'.amc'):
                    os.system(f"am {am_template}  {El}  {path_to_dat}")
                    print('am fit done for file '+path_to_dat)
                file_variable = open(path_to_dat+'.amc')
                all_lines_variable = file_variable.readlines()
                nscale_line=all_lines_variable[8]
                std_line=all_lines_variable[10]
                pwv_line=all_lines_variable[504]
                print('pwv_line=', pwv_line)
                Nscale=float(nscale_line[33:40])
                Nscale_std=float(std_line[34:38])
                print('pwv=', pwv_line[30:36])
                try:
                    pwv=float(pwv_line[30:37])
                except:
                    pwv=float(pwv_line[30:36])
                else:
                    pwv=float(pwv_line[30:35])
                print('pwv=', pwv)
                pwv_ts[azn]=pwv
                nscale_ts[azn]=Nscale
                nscale_err_ts[azn]=Nscale_std
            except Exception as e:
               print(f'{path_to_dat} failed: {e}')

        fit_output['Nscale']=nscale_ts
        fit_output['Nscale_err']=nscale_err_ts
        fit_output['pwv']=pwv_ts
        fit_output['el_correction']=el_correction

        if os.path.exists(pickle_fn):
            os.remove(pickle_fn)
            print('Deleting old pickle file '+pickle_fn)

        print('Writing on Pickle File '+pickle_fn)

        f = open(pickle_fn,'wb')
        pickle.dump(fit_output, f)
        f.close()

        fig=plt.figure()
        im=plt.imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
        plt.suptitle('PWV Atmogram\n'+filename[:-4])
        cbar = fig.colorbar(im, extend='both')
        cbar.set_label('PWV[um]')
        plt.savefig(out_path+filename[:-4]+'_pwvatmo_clean_mod'+str(clean_mod)+'_clean_method_'+str(clean_method)+'.png')
        plt.close()

        return(fit_output)

    def fit_w_am_zenith(self, filename, path_to_data='wvr1_data/', template= 'SPole_annual_50.amc', path_to_temp='Templates/'):
        Npoints=0

        if not os.path.exists('am_datfiles/'+template[:-4]):
            os.makedirs('am_datfiles/'+template[:-4])

        path='am_datfiles/'+template[:-4]+'/'+filename[:-4]


        if not os.path.exists(path):
            os.makedirs(path)

        pickle_fn=path+'/'+filename[:-4]+'_fitoutput_corr.txt'
        pickle_fn_read=path+'/'+filename[:-4]+'_pickle_temps_corr.txt'

        fit_output = {'filename': [], 'El': [], 'dir':[], 'Nscale': [], 'Nscale_err': [],
        'pwv_meso':[], 'pwv_strato':[], 'pwv_tropo':[], 'pwv_total':[], 'pwv_los_total':[]}

        end_string='rr.dat'
        end_string_fn='_corr'


        # try:
        #     t = open(pickle_fn_read, 'rb')
        #     print(pickle_fn_read+' exists.')
        #     pickle_Temps = pickle.load(t)
        #     t.close()
        #     print('pickle_Temps =', pickle_Temps)
        # except:
        self.create_am_datfile(filename, path_to_data=path_to_data, template=template, spline=2)
        print('dat file created')


        t = open(pickle_fn_read, 'rb')
        pickle_Temps = pickle.load(t)
        t.close()

        am_dat_list=pickle_Temps['am_dat_filename']
        el_list=pickle_Temps['El']
        dir_list=pickle_Temps['dir']



        T0_check=[]
        for findex in range (len(am_dat_list)):
            flist=am_dat_list[findex]
            if flist[-6:]==end_string:
                path_to_dat=flist
                print('El_max =', np.max(el_list))
                el= el_list[findex]
                print('El=',el)
                zen_ang=90-el
                dir=flist[35:37]
                dir=dir_list[findex]
                if (el==np.max(el_list)):
                    try:
                        print('fitting.')
                        am_template=path_to_temp+template
                        if not os.path.exists(path_to_dat+'.amc'):
                            os.system(f"am {am_template}  {zen_ang}  {path_to_dat}")
                        dat_file_variable = open(path_to_dat)
                        dat_lines = dat_file_variable.readlines()
                        T0_line=dat_lines[0]
                        T0_check.append(T0_line[10:15])
                        file_variable = open(path_to_dat+'.amc')
                        all_lines_variable = file_variable.readlines()
                        nscale_line=all_lines_variable[44]
                        std_line=all_lines_variable[10]
                        if template[0:5]=='MERRA':
                            #pwv_line=all_lines_variable[1062] #wo dry layers
                            #pwv_line=all_lines_variable[1632]#[1633]
                            pwv_line_meso=all_lines_variable[1614]
                            pwv_line_strato=all_lines_variable[1621]#[1633]
                            pwv_line_tropo=all_lines_variable[1628]#[1633]
                            pwv_line_total=all_lines_variable[1633]#[1633]
                        else:
                            #pwv_line=all_lines_variable[503]#[504]
                            pwv_line_meso=all_lines_variable[483]
                            pwv_line_strato=all_lines_variable[490]#[1633]
                            pwv_line_tropo=all_lines_variable[497]#[1633]
                            pwv_line_total=all_lines_variable[504]#[1633]
                        print('pwv_line_tropo=', pwv_line_tropo)
                        try:
                            Nscale=float(nscale_line[23:30])
                        except:
                            Nscale=float(nscale_line[23:28])
                        print('nscale=', Nscale)
                        try:
                            Nscale_std=float(std_line[34:38])
                        except:
                            Nscale_std=float(std_line[34:36])
                        print('nscale_std=', Nscale_std)
                        try:
                            pwv_meso=float(pwv_line_meso[30:36])
                        except:
                            pwv_meso=float(pwv_line_meso[30:34])
                        try:
                            pwv_strato=float(pwv_line_strato[30:36])
                        except:
                            pwv_strato=float(pwv_line_strato[30:34])
                        try:
                            pwv_tropo=float(pwv_line_tropo[30:36])
                        except:
                            pwv_tropo=float(pwv_line_tropo[30:34])
                        try:
                            pwv_total=float(pwv_line_total[30:36])
                        except:
                            pwv_total=float(pwv_line_total[30:34])
                        print('pwv_total=', pwv_total)
                        try:
                            pwv_los_total=float(pwv_line_total[49:56])
                        except:
                            pwv_los_total=float(pwv_line_total[49:54])
                        print('pwv_los=', pwv_los_total)
                        #print('T0_check=', T0_line[10:15])
                        fit_output['El'].append(el)
                        fit_output['dir'].append(dir)
                        fit_output['filename'].append(flist)
                        fit_output['Nscale'].append(Nscale)
                        fit_output['Nscale_err'].append(Nscale_std)
                        fit_output['pwv_meso'].append(pwv_meso)
                        fit_output['pwv_strato'].append(pwv_strato)
                        fit_output['pwv_tropo'].append(pwv_tropo)
                        fit_output['pwv_total'].append(pwv_total)
                        fit_output['pwv_los_total'].append(pwv_los_total)
                        Npoints+=1
                    except Exception as e:
                        print(f'{flist} failed: {e}')


        if os.path.exists(pickle_fn):
            os.remove(pickle_fn)
            print('Deleting old pickle file '+pickle_fn)

        print('Writing on Pickle File '+pickle_fn)

        f = open(pickle_fn,'wb')
        pickle.dump(fit_output, f)
        f.close()

        return(fit_output)


    def plot_am_fit_2(self, filename, var='pwv', path_to_amc='', pwv_layer='total', template= 'SPole_annual_50.amc', el_range=(0,90.), spline=0, pf='None', pf2='None', show=0):
        if not os.path.exists(f'fit_output_plots/'+template[:-4]):
            os.makedirs(f'fit_output_plots/'+template[:-4])
        dirpath='am_datfiles/'+path_to_amc+template[:-4]+'/'+filename[:-4]+'/'

        if spline==0:

            pickle_fn=dirpath+filename[:-4]+'_fitoutput.txt'
            pickle_fn_temps=dirpath+filename[:-4]+'_pickle_temps.txt'

            paths = sorted(Path(dirpath).iterdir(), key=os.path.getmtime)
            check_path=paths[0:10]
            Nscale_list=np.zeros(len(paths))
            Nscale_std_list=np.zeros(len(paths))
            pwv_list=np.zeros(len(paths))
            T0_list=np.zeros(len(paths))
            T1_list=np.zeros(len(paths))
            T2_list=np.zeros(len(paths))
            T3_list=np.zeros(len(paths))
            El_list=np.zeros(len(paths))
            dir_list=['None']*len(paths)
            dir_list=np.array(dir_list)

            i_list=0
            el=0

            f = open(pickle_fn, 'rb')
            fit_output = pickle.load(f)
            f.close()

            t = open(pickle_fn_temps, 'rb')
            pickle_Temps = pickle.load(t)
            t.close()

            Nscale_list=fit_output['Nscale']
            Nscale_std_list=fit_output['Nscale_err']
            pwv_list=fit_output['pwv']
            El_list=fit_output['El']
            dir_list=np.array(fit_output['dir'])

            El_list_T=pickle_Temps['El']
            dir_list_T=np.array(pickle_Temps['dir'])

            T0_list=np.array(pickle_Temps['T0'])
            T1_list=np.array(pickle_Temps['T1'])
            T2_list=np.array(pickle_Temps['T2'])
            T3_list=np.array(pickle_Temps['T3'])

            dir_col=np.array([None]*len(dir_list))
            dir_mkr=np.array([None]*len(dir_list))
            pwv_list=np.array(pwv_list)
            El_list=np.array(El_list)
            El_list_T=np.array(El_list_T)


            def mscatter(x,y,ax=None, m=None, **kw):
                import matplotlib.markers as mmarkers
                if not ax: ax=plt.gca()
                sc = ax.scatter(x,y,**kw)
                if (m is not None) and (len(m)==len(x)):
                    paths = []
                    for marker in m:
                        if isinstance(marker, mmarkers.MarkerStyle):
                            marker_obj = marker
                        else:
                            marker_obj = mmarkers.MarkerStyle(marker)
                        path = marker_obj.get_path().transformed(
                            marker_obj.get_transform())
                        paths.append(path)
                    sc.set_paths(paths)
                return sc

            pwv_up=pwv_list[np.where(dir_list=='up')]
            pwv_dn=pwv_list[np.where(dir_list=='dn')]

            T0_up=T0_list[np.where(dir_list_T=='up')]
            T0_dn=T0_list[np.where(dir_list_T=='dn')]

            El_up=El_list[np.where(dir_list=='up')]
            El_dn=El_list[np.where(dir_list=='dn')]

            El_up_T=El_list_T[np.where(dir_list_T=='up')]
            El_dn_T=El_list_T[np.where(dir_list_T=='dn')]

            print('El_list_T_shape=', np.shape(El_list_T))
            print('dir_list_T=', dir_list_T)

            print('pwv=', pwv_list)
            print('T0_list=', T0_list)

            print('pwv_up=', pwv_up)
            print('pwv_dn=', pwv_dn)


            pl.figure(figsize=(16,10))
            pl.scatter(El_up, pwv_up, color='r', s=55, marker='^', label='WVR moving upwards')
            pl.scatter(El_dn, pwv_dn, color='b', s=55, marker='v', label='WVR moving downwards')
            pl.xlabel('El[deg]', fontsize=16)
            pl.ylabel('pwv[um]', fontsize=16)
            pl.text(np.max(El_list)-30.,  np.max(pwv_list)-3., 'pwv_mean[um]={:.2f}'.format(np.mean(pwv_list))+'±{:.2f}'.format(np.std(pwv_list)), fontsize=16)
            pl.text(np.max(El_list)-30.,  np.max(pwv_list)-6., 'pwv_mean_up[um]={:.2f}'.format(np.mean(pwv_up))+'±{:.2f}'.format(np.std(pwv_up)), color='r', fontsize=16)
            pl.text(np.max(El_list)-30.,  np.max(pwv_list)-9., 'pwv_mean_dn[um]={:.2f}'.format(np.mean(pwv_dn))+'±{:.2f}'.format(np.std(pwv_dn)), color='b', fontsize=16)
            pl.axhline(y=np.mean(pwv_list), c='k',alpha=0.6, linestyle='--')
            pl.axhline(y=np.mean(pwv_up), c='r',alpha=0.6, linestyle='--')
            pl.axhline(y=np.mean(pwv_dn), c='b',alpha=0.6, linestyle='--')
            pl.legend(loc='lower left', fontsize=14)
            pl.suptitle(filename[:-4]+'\nTemplate: '+template[:-4], fontsize=16)
            pl.title('Before BL+Z_offset correction', fontsize=16)

            pl.savefig('fit_output_plots/'+template[:-4]+'/'+filename[:-4]+'_fitoutput_pwv_wlineshape.png')
            if pf != 'None':
                pl.savefig(pf+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_wlineshape.png')
            if pf2 != 'None':
                pl.savefig(pf2+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_wlineshape.png')

            #pl.close()
            pl.close()

        if spline==1:

            pickle_fn=dirpath+filename[:-4]+'_fitoutput_red.txt'
            pickle_fn_temps=dirpath+filename[:-4]+'_pickle_temps_red.txt'
            pickle_fn_betas=dirpath+filename[:-4]+'_pickle_betas.txt'

            paths = sorted(Path(dirpath).iterdir(), key=os.path.getmtime)
            check_path=paths[0:10]
            Nscale_list=np.zeros(len(paths))
            Nscale_std_list=np.zeros(len(paths))
            pwv_list=np.zeros(len(paths))
            T0_list=np.zeros(len(paths))
            T1_list=np.zeros(len(paths))
            T2_list=np.zeros(len(paths))
            T3_list=np.zeros(len(paths))
            El_list=np.zeros(len(paths))
            dir_list=['None']*len(paths)
            dir_list=np.array(dir_list)


            i_list=0
            el=0

            f = open(pickle_fn, 'rb')
            fit_output = pickle.load(f)
            f.close()

            t = open(pickle_fn_temps, 'rb')
            pickle_Temps = pickle.load(t)
            t.close()

            t2 = open(pickle_fn_betas, 'rb')
            betas = pickle.load(t2)
            t2.close()


            Nscale_list=fit_output['Nscale']
            Nscale_std_list=fit_output['Nscale_err']
            pwv_list=fit_output[var]
            El_list=fit_output['El']
            dir_list=np.array(fit_output['dir'])

            El_list_T=pickle_Temps['El']
            dir_list_T=np.array(pickle_Temps['dir'])

            T0_list=np.array(pickle_Temps['T0'])
            T1_list=np.array(pickle_Temps['T1'])
            T2_list=np.array(pickle_Temps['T2'])
            T3_list=np.array(pickle_Temps['T3'])

            pwv_list=np.array(pwv_list)
            El_list=np.array(El_list)
            El_list_T=np.array(El_list_T)
            print('El_list=', El_list)
            print('pwv_list=', pwv_list)


            def line_gauss(x, c, a, sigma):
                return c + a * np.exp(-(x - 183.3)**2 / (2 * sigma**2))  #to understand why it doesn't work if i define it out of here

            freqs_fit=np.linspace(170, 190, 2000)

            for i in range(len(El_list)):
                T_dualside=[T3_list[i], T2_list[i], T1_list[i], T0_list[i], T0_list[i], T1_list[i], T2_list[i], T3_list[i]]
                popt,pcov = sp.curve_fit(line_gauss, self.freqs, T_dualside, p0 = [1, 10, 10], maxfev=1000000)
                fig = plt.figure(figsize=(16, 12))
                pl.scatter(self.freqs, T_dualside, c='r')
                pl.xlabel('f[GHz]')
                pl.ylabel('T[K]')
                pl.ylim(np.min(line_gauss(freqs_fit, *popt))-5., np.max(line_gauss(freqs_fit, *popt))+5.)
                pl.plot(freqs_fit, line_gauss(freqs_fit, *popt), c='k', linewidth=1, label='Gaussian Fit to Lineshape in T')
                pl.legend()
                pl.xlim(175,190)
                pl.suptitle('Lineshape in T', fontsize=25)
                pl.title(filename[:-4]+' El='+str(El_list[i]), fontsize=24)

                pl.savefig('fit_output_plots/'+template[:-4]+'/'+filename[:-4]+'_El'+str(El_list[i])+'_T_wlineshape.png')

                #pl.close()
                pl.close()


            El55_up=np.where(El_list>=54.5)[0]
            print('El55_up=', El55_up)
            print('El55=', El55_up[0])


            beta_dualside=[betas[3], betas[2], betas[1], betas[0], betas[0], betas[1], betas[2], betas[3]]
            popt,pcov = sp.curve_fit(line_gauss, self.freqs, beta_dualside, p0 = [1, 10, 10], maxfev=1000000)

            fig, ax = plt.subplots(2, 1, figsize=(16, 12))

            ax[0].scatter(self.freqs, beta_dualside, c='r', marker='*', label='Opacities from Single Slab Model')
            ax[0].set_xlabel('f[GHz]')
            ax[0].set_ylabel('Betas_at_Zenith')
            ax[0].set_ylim(np.min(betas)-0.2, np.max(betas)+0.2)
            #ax[0].plot(freqs_fit, line_gauss(freqs_fit, *popt), linewidth=1, label='Gaussian Fit to Lineshape')
            ax[0].legend()
            ax[0].set_xlim(175,190)

            ax[1].scatter(El_list, pwv_list, color='b', s=60, marker='o')
            ax[1].set_xlabel('El[deg]', fontsize=22)
            ax[1].set_ylabel('pwv[um]', fontsize=22)
            #ax[1].set_xticks(fontsize=18)
            #ax[1].set_yticks(fontsize=18)
            ax[1].text(np.max(El_list)-20.,  np.max(pwv_list)-3., 'pwv_mean[um]={:.2f}'.format(np.mean(pwv_list)), fontsize=20)
            ax[1].text(np.max(El_list)-20.,  np.max(pwv_list)-6., 'pwv_std[um]={:.2f}'.format(np.std(pwv_list)), fontsize=20)
            ax[1].text(np.max(El_list)-20.,  np.max(pwv_list)-9., 'pwv_El55[um]={:.2f}'.format(pwv_list[El55_up[0]]), fontsize=20)
            ax[1].axhline(y=np.mean(pwv_list), c='k', alpha=0.6, linestyle='--')

            pl.suptitle(filename[:-4]+'\nTemplate: '+template[:-4], fontsize=20)
            pl.title('After BL+Z_offset correction and Data Reduction', fontsize=20)

            pl.savefig('fit_output_plots/'+template[:-4]+'/'+filename[:-4]+'_fitoutput_pwv_red_20_wlineshape.png')
            if pf != 'None':
                pl.savefig(pf+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_red_20_wlineshape.png')
            if pf2 != 'None':
                pl.savefig(pf2+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_red_20_wlineshape.png')

            #pl.close()
            pl.close()


        if spline==2:

            pl.figure(figsize=(16,10))
            #Modified for template fit comparison
            #temp_list=['SPole_annual_50.amc', 'MERRA_'+filename[:4]+'-'+filename[4:6]+'-'+filename[6:8]+'_'+filename[9:11]+'.amc']
            # col_list=['r', 'b']
            # col_star_list=['y', 'g']

            temp_list=[template]
            col_list=['r']
            col_star_list=['y']

            for i in range(len(temp_list)):
                temp=temp_list[i]
                dirpath='am_datfiles/'+path_to_amc+temp[:-4]+'/'+filename[:-4]+'/'
                col=col_list[i]
                col_star=col_star_list[i]

                pickle_fn=dirpath+filename[:-4]+'_fitoutput_corr.txt'
                pickle_fn_temps=dirpath+filename[:-4]+'_pickle_temps_corr.txt'
                pickle_fn_betas=dirpath+filename[:-4]+'_pickle_betas.txt'

                paths = sorted(Path(dirpath).iterdir(), key=os.path.getmtime)
                check_path=paths[0:10]
                Nscale_list=np.zeros(len(paths))
                Nscale_std_list=np.zeros(len(paths))
                pwv_list=np.zeros(len(paths))
                T0_list=np.zeros(len(paths))
                T1_list=np.zeros(len(paths))
                T2_list=np.zeros(len(paths))
                T3_list=np.zeros(len(paths))
                El_list=np.zeros(len(paths))
                dir_list=['None']*len(paths)
                dir_list=np.array(dir_list)

                i_list=0
                el=0

                f = open(pickle_fn, 'rb')
                fit_output = pickle.load(f)
                f.close()

                t = open(pickle_fn_temps, 'rb')
                pickle_Temps = pickle.load(t)
                t.close()

                Nscale_list=fit_output['Nscale']
                Nscale_std_list=fit_output['Nscale_err']
                try:
                    pwv_list=fit_output['pwv_'+pwv_layer]
                    pwv_los_list=fit_output['pwv_los_total']
                except:
                    pwv_list=fit_output['pwv']

                El_list=fit_output['El']
                dir_list=np.array(fit_output['dir'])

                El_list_T=pickle_Temps['El']
                dir_list_T=np.array(pickle_Temps['dir'])

                T0_list=np.array(pickle_Temps['T0'])
                T1_list=np.array(pickle_Temps['T1'])
                T2_list=np.array(pickle_Temps['T2'])
                T3_list=np.array(pickle_Temps['T3'])

                dir_col=np.array([None]*len(dir_list))
                dir_mkr=np.array([None]*len(dir_list))
                pwv_list=np.array(pwv_list)
            #     pwv_los_list=np.array(pwv_los_list)
                El_list=np.array(El_list)
                El_list_T=np.array(El_list_T)
                print('El_list_T=', El_list_T)

                pwv_up=pwv_list[np.where(dir_list=='up')]
                pwv_dn=pwv_list[np.where(dir_list=='dn')]

                # pwv_los_up=pwv_los_list[np.where(dir_list=='up')]
                # pwv_los_dn=pwv_los_list[np.where(dir_list=='dn')]

                T0_up=T0_list[np.where(dir_list_T=='up')]
                T0_dn=T0_list[np.where(dir_list_T=='dn')]

                El_up=El_list[np.where(dir_list=='up')]
                El_dn=El_list[np.where(dir_list=='dn')]

                El_up_T=El_list_T[np.where(dir_list_T=='up')]
                El_dn_T=El_list_T[np.where(dir_list_T=='dn')]

                print('El_list_T_shape=', np.shape(El_list_T))
                print('dir_list_T_shape=', np.shape(dir_list_T))
                print('T0_list_shape=', np.shape(T0_list))

                pwv_start=[pwv_list[0]]
                pwv_end=[pwv_list[len(El_list)-1]]
                el_start=[El_list[0]]
                el_end=[El_list[len(El_list)-1]]

                if var=='nscale':

                    pl.scatter(El_list, Nscale_list, color=col, s=45, label=temp)
                    pl.plot(El_list, Nscale_list, color=col, alpha=0.6)
                    pl.axhline(y=np.mean(Nscale_list), c=col, alpha=0.6, linestyle='--', label=('nscale_avg='+str(round(np.mean(Nscale_list),2))))

                if var=='pwv':

                    pl.scatter(El_up, pwv_up, color=col, s=55, marker='^', label=temp)
                    #pl.scatter(El_up, pwv_los_up, color='r', s=55, marker='^', label='WVR moving upwards - PWV_los')
                    pl.plot(El_list, pwv_list, color=col, alpha=0.5)
                    pl.scatter(El_dn, pwv_dn, color=col, s=55, marker='v')
                    #pl.scatter(El_dn, pwv_los_dn, color='r', s=55, marker='v', label='WVR moving downwards - PWV_los')
                    pl.scatter(el_start, pwv_start, color=col_star, marker='*', s=55, label='Start')
                    pl.scatter(el_end, pwv_end, color=col_star, s=55, marker='X', label='End')

                    #pl.text(np.max(El_list)-30.,  np.max(pwv_list)-3., 'pwv_mean[um]={:.2f}'.format(np.mean(pwv_list))+'±{:.2f}'.format(np.std(pwv_list)), fontsize=16)
                    #pl.text(np.max(El_list)-30.,  np.max(pwv_list)-6., 'pwv_mean_up[um]={:.2f}'.format(np.mean(pwv_up))+'±{:.2f}'.format(np.std(pwv_up)), color='r', fontsize=16)
                    #pl.text(np.max(El_list)-30.,  np.max(pwv_list)-9., 'pwv_mean_dn[um]={:.2f}'.format(np.mean(pwv_dn))+'±{:.2f}'.format(np.std(pwv_dn)), color='b', fontsize=16)
                    #pl.axhline(y=np.mean(pwv_list), c=col,alpha=0.6, linestyle='--', label=('pwv_avg='+str(round(np.mean(pwv_list),2))))


            if var=='pwv':
                pl.xlabel('El[deg]', fontsize=16)
                pl.ylabel('pwv[um]', fontsize=16)
                pl.legend(loc='upper right', fontsize=14)
                pl.suptitle(filename[:-4]+'\nTemplate: '+template[:-4], fontsize=16)
                pl.title('Layer: '+pwv_layer, fontsize=16)
                pl.savefig('fit_output_plots/'+template[:-4]+'/'+filename[:-4]+'_fitoutput_pwv_wlineshape_corr.png')
                if pf != 'None':
                    pl.savefig(pf+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_wlineshape_corr.png')
                if pf2 != 'None':
                    pl.savefig(pf2+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_wlineshape_corr.png')
                if show==1:
                    pl.show()
                else:
                    pl.close()

            if var=='nscale':
                pl.xlabel('El[deg]', fontsize=16)
                pl.ylabel('Nscale[um]', fontsize=16)
                pl.legend(loc='upper right', fontsize=14)
                pl.suptitle(filename[:-4]+'\nTemplate: '+template[:-4], fontsize=16)
                pl.title('Layer: '+pwv_layer, fontsize=16)
                pl.savefig('fit_output_plots/'+template[:-4]+'/'+filename[:-4]+'_fitoutput_nscale_wlineshape_corr.png')
                if pf != 'None':
                    pl.savefig(pf+filename[:-4]+'_fitoutput_nscale_'+pwv_layer+'_wlineshape_corr.png')
                if pf2 != 'None':
                    pl.savefig(pf2+filename[:-4]+'_fitoutput_nscale_'+pwv_layer+'_wlineshape_corr.png')
                if show==1:
                    pl.show()
                else:
                    pl.close()


        return El_list, pwv_list



    def read_amc(self, filename, pathtofn='', full_spectrum=1, plot=1):

        if full_spectrum==1:

            famc=open(pathtofn+filename,"r")
            amclines=famc.readlines()

            amclines[27] = "f 10 GHz  350 GHz  200 MHz\n"
            amclines[28]=''
            amclines[29]=''

            newfile = open(pathtofn+filename[:-4]+'_fullspectrum.amc', "w")
            newfile.writelines(amclines)
            newfile.close()

            filename=filename[:-4]+'_fullspectrum.amc'

        os.system('am '+pathtofn+filename+' >'+pathtofn+filename[:-4]+'.out')
        outfile=pathtofn+filename[:-4]+'.out'
        f=open(outfile,"r")
        lines=f.readlines()

        out_dict={'freq': [], 'tau':[], 'Trj':[]}

        for x in lines:
            out_dict['freq'].append(float(x.split(' ')[0]))
            out_dict['tau'].append(float(x.split(' ')[1]))
            out_dict['Trj'].append(float(x.split(' ')[2]))
        f.close()

        f=open(outfile[:-4]+'.txt',"wb")
        pickle.dump(out_dict, f)
        f.close()

        if plot==1:
            plt.plot(out_dict['freq'], out_dict['tau'])
            plt.ylabel('tau')
            plt.xlabel('Freq[GHz]')
            plt.title(filename[:-4])
            plt.savefig(outfile[:-4]+'_tau.png')
            plt.show()

            plt.plot(out_dict['freq'], out_dict['Trj'])
            plt.ylabel('Trj[K]')
            plt.xlabel('Freq[GHz]')
            plt.title(filename[:-4])
            plt.savefig(outfile[:-4]+'_Trj.png')
            plt.show()

        return out_dict

    def plot_Trj_el(self, day, remake_post_fig=1, rewrite_txt=1, datafolder='', posting_folder=''):
            path='am_datfiles/SPole_annual_50/'
            fold_day=path+str(day)+'_bk_spectrum/'

            def write_txt():
                scan_list=[]
                for fn in os.listdir(path):
                    if fn[-11:]=='skyDip_fast':
                        scan_list.append(fn)

                T_spectrum_El55={'day':[], 'time':[], '30GHz':[], '40GHz':[], '90GHz':[], '150GHz':[], '220GHz':[], '270GHz':[]}

                for scan in scan_list:
                    folder=path+scan
                    dat_exist=1

                    pk_pwv=path+scan+'/'+scan+'_fitoutput_corr.txt'
                    # t = open(pk_pwv, 'rb')
                    # fit_output = pickle.load(t)
                    # t.close()

                    if not os.path.exists(folder+'/spectra/plots'):
                        os.makedirs(folder+'/spectra/plots')

                    test_data=scan+'.txt'
                    path_to_test=''

                    if not os.path.exists(pk_pwv):
                        try:
                            self.create_am_datfile(test_data, path_to_data=path_to_test, spline=2)
                            self.fit_w_am(test_data, path_to_data=path_to_test, spline=2)
                        except:
                            print('dat file '+test_data+' not found.')
                            dat_exist=0
                    print('dat exist=', dat_exist)
                    if dat_exist==1:
                        T_spectrum_El55['day'].append(scan[:8])
                        T_spectrum_El55['time'].append(scan[9:15])

                        Trj_dict={'fn': [], 'el':[], 'Trj':[], 'pwv': []}
                        for flist in os.listdir(folder):
                            print('flist:', flist)
                            if (flist[-4:]=='.amc') & (flist[-12:-8]=='corr'):
                                if flist[30:32]=='55':
                                    print(flist)
                                    print(flist[-11:-8])
                                    amcfile=flist

                                    a_file = open(folder+'/'+amcfile, "r")
                                    list_of_lines = a_file.readlines()
                                    print(list_of_lines[27])
                                    list_of_lines[27] = "f 10 GHz  350 GHz  200 MHz\n"
                                    list_of_lines[28]=''
                                    list_of_lines[29]=''
                                    print(list_of_lines[27])
                                    b_file = open(folder+'/spectra/'+amcfile[:-4]+'_bk.amc', "w")
                                    b_file.writelines(list_of_lines)
                                    b_file.close()

                                    newamcfile=amcfile[:-4]+'_bk.amc'

                                    os.system('am '+folder+'/spectra/'+newamcfile+' >'+folder+'/spectra/'+newamcfile[:-4]+'.out')

                                    outfile=folder+'/spectra/'+newamcfile[:-4]+'.out'
                                    f=open(outfile,"r")
                                    lines=f.readlines()
                                    print(lines[0])

                                    freq=[]
                                    tau=[]
                                    Trj=[]
                                    for x in lines:
                                        freq.append(float(x.split(' ')[0]))
                                        tau.append(float(x.split(' ')[1]))
                                        Trj.append(float(x.split(' ')[2]))
                                    f.close()

                                    Trj_f=np.zeros(6)

                                    fig, ax = pl.subplots(6,1)
                                    avg30=np.mean(np.array(Trj)[np.where((np.array(freq)>25)&(np.array(freq)<35))])
                                    ax[0].plot(freq, Trj, label='Trj[30GHz]='+str(round(avg30,2)))
                                    Trj_f[0]=round(avg30,2)
                                    ax[0].set_ylabel('T_rj[K]')
                                    ax[0].set_xlim(25,35)
                                    int0=np.where((np.array(freq)>25)&(np.array(freq)<35))
                                    print(np.array(Trj))
                                    ax[0].set_ylim(np.min(np.array(Trj)[int0[0]]), np.max(np.array(Trj)[int0[0]]))
                                    ax[0].legend()

                                    avg40=np.mean(np.array(Trj)[np.where((np.array(freq)>35)&(np.array(freq)<=45))])
                                    ax[1].plot(freq, Trj, label='Trj[40GHz]='+str(round(avg40,2)))
                                    Trj_f[1]=round(avg40,2)
                                    ax[1].set_ylabel('T_rj[K]')
                                    ax[1].set_xlim(35,45)
                                    int1=np.where((np.array(freq)>35)&(np.array(freq)<45))
                                    print(np.array(Trj)[int1])
                                    ax[1].set_ylim(np.min(np.array(Trj)[int1[0]]), np.max(np.array(Trj)[int1[0]]))
                                    ax[1].legend()

                                    avg90=np.mean(np.array(Trj)[np.where((np.array(freq)>80.9)&(np.array(freq)<107.5))]) #BC=94.2; BW=26.7
                                    ax[2].plot(freq, Trj, label='Trj[90GHz]='+str(round(avg90,2)))
                                    Trj_f[2]=round(avg90,2)
                                    ax[2].set_ylabel('T_rj[K]')
                                    ax[2].set_xlim(80,108)
                                    int2=np.where((np.array(freq)>80.9)&(np.array(freq)<107.5))
                                    print(np.array(Trj)[int2])
                                    ax[2].set_ylim(np.min(np.array(Trj)[int2[0]]), np.max(np.array(Trj)[int2[0]]))
                                    ax[2].legend()

                                    avg150=np.mean(np.array(Trj)[np.where((np.array(freq)>127.3)&(np.array(freq)<170.7))]) #BC=149;BW=43.4
                                    ax[3].plot(freq, Trj, label='Trj[150GHz]='+str(round(avg150,2)))
                                    Trj_f[3]=round(avg150,2)
                                    ax[3].set_xlabel('freq[GHz]')
                                    ax[3].set_ylabel('T_rj[K]')
                                    ax[3].set_xlim(127, 171)
                                    int3=np.where((np.array(freq)>127.3)&(np.array(freq)<170.7))
                                    print(np.array(Trj)[int3])
                                    ax[3].set_ylim(np.min(np.array(Trj)[int3[0]]), np.max(np.array(Trj)[int3[0]]))
                                    ax[3].legend()

                                    avg220=np.mean(np.array(Trj)[np.where((np.array(freq)>205.6)&(np.array(freq)<256.8))]) #BC=231.2; BW=51.2
                                    ax[4].plot(freq, Trj, label='Trj[220GHz]='+str(round(avg220,2)))
                                    Trj_f[4]=round(avg220,2)
                                    ax[4].set_ylabel('T_rj[K]')
                                    ax[4].set_xlim(205,257)
                                    int4=np.where((np.array(freq)>205.6)&(np.array(freq)<256.8))
                                    print(np.array(Trj)[int4])
                                    ax[4].set_ylim(np.min(np.array(Trj)[int4[0]]), np.max(np.array(Trj)[int4[0]]))
                                    ax[4].legend()

                                    avg270=np.mean(np.array(Trj)[np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))])#BC=275.5;BW=69.8
                                    ax[5].plot(freq, Trj, label='Trj[270GHz]='+str(round(avg270,2)))
                                    Trj_f[5]=round(avg270,2)
                                    ax[5].set_xlabel('freq[GHz]')
                                    ax[5].set_ylabel('T_rj[K]')
                                    ax[5].set_xlim(240,311)
                                    int5=np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))
                                    print(np.array(Trj)[int5])
                                    ax[5].set_ylim(np.min(np.array(Trj)[int5[0]]), np.max(np.array(Trj)[int5[0]]))
                                    ax[5].legend()

                                    if not os.path.exists(folder+'/spectra/plots/'):
                                        os.makedirs(folder+'/spectra/plots/')
                                    if remake_post_fig==1:
                                        plt.savefig(posting_folder+'spectra/TrjinBKbands_'+newamcfile[:-4]+'.png')
                                    plt.close()

                                    plt.plot(freq, Trj)
                                    plt.ylabel('T_rj[K]')
                                    plt.xlabel('Freq[GHz]')
                                    plt.title(newamcfile[:-4])
                                    if remake_post_fig==1:
                                        plt.savefig(posting_folder+'spectra/TrjinBKbands_fullspectrum_'+newamcfile[:-4]+'.png')
                                    plt.savefig(folder+'/spectra/plots/'+newamcfile[:-4]+'_fullspectrum.png')
                                    plt.close()

                                    Trj_dict['fn'].append(newamcfile[:-4]+'_Trj_pk.txt')
                                    Trj_dict['el'].append(float(newamcfile[30:34]))
                                    Trj_dict['Trj'].append(Trj_f)

                                print(Trj_dict)
                                f=open(fold_day+'/spectra/Trj_el_'+day+'_pk.txt',"wb")
                                pickle.dump(Trj_dict, f)
                                f.close()

            def read_txt():
                f=open(fold_day+'/spectra/Trj_el_'+day+'_pk.txt',"rb")
                Trj_dict=pickle.load(f)
                f.close()

                print(Trj_dict)
                print(Trj_dict['el'][0])

                El_list=[]
                T30=[]
                T40=[]
                T90=[]
                T150=[]
                T210=[]
                T270=[]

                for i in range (len(Trj_dict['fn'])):
                    El_list.append(Trj_dict['el'][i])
                    T=Trj_dict['Trj'][i]
                    print(T)
                    T30.append(T[0])
                    T40.append(T[1])
                    T90.append(T[2])
                    T150.append(T[3])
                    T210.append(T[4])
                    T270.append(T[5])

                # fig, ax = pl.subplots(2,1, sharex=True)
                #
                # ax[0].plot(El_list, T30, label='30GHz')
                # ax[0].plot(El_list, T40, label='40GHz')
                # ax[0].plot(El_list, T90, label='90GHz')
                # ax[0].plot(El_list, T150, label='150GHz')
                # ax[0].plot(El_list, T210, label='220GHz')
                # ax[0].plot(El_list, T270, label='270GHz')
                # # ax[0].xlabel('El[deg]')
                # ax[0].set_ylabel('T_rj[K]')
                # ax[0].legend()
                #
                # ax[1].scatter(fit_output['El'], fit_output['pwv'], s=5)
                # ax[1].plot(fit_output['El'], fit_output['pwv'], c='k', alpha=0.5)
                # ax[1].set_xlabel('El[deg]')
                # ax[1].set_ylabel('PWV[um]')
                # ax[1].set_xlim(55.,55.9)
                #
                # plt.suptitle(scan)
                # if remake_post_fig==1:
                #     plt.savefig(posting_folder+'TrjandPWV_'+scan+'.png')
                # plt.savefig(folder+'/'+scan+'.png')
                # #plt.show()
                # plt.close()

                T_spectrum_El55['30GHz'].append(T30[0])
                T_spectrum_El55['40GHz'].append(T40[0])
                T_spectrum_El55['90GHz'].append(T90[0])
                T_spectrum_El55['150GHz'].append(T150[0])
                T_spectrum_El55['220GHz'].append(T210[0])
                T_spectrum_El55['270GHz'].append(T270[0])

                if not os.path.exists(path+day+'_bk_spectrum'):
                    os.makedirs(path+day+'_bk_spectrum')

                f=open(path+day+'_bk_spectrum'+'/T_spectrum_El55.txt',"wb")
                pickle.dump(T_spectrum_El55, f)
                f.close()

                f=open(path+day+'_bk_spectrum'+'/T_spectrum_El55.txt',"rb")
                T_spectrum_El55=pickle.load(f)
                f.close()

                time_axis=[]
                for j in range (len(T_spectrum_El55['time'])):
                    t=T_spectrum_El55['time'][j]
                    time_axis.append(t[:2])

                print('time_axis=', time_axis)
                print('shape=', np.shape(time_axis))
                print('T_spectrum_El55[\'30GHz\']=', T_spectrum_El55['30GHz'])
                print('shape=', np.shape(T_spectrum_El55['30GHz']))

                plt.scatter(time_axis, T_spectrum_El55['30GHz'], label='30GHz')
                plt.plot(time_axis, T_spectrum_El55['30GHz'], c='k', alpha=0.5)
                plt.scatter(time_axis, T_spectrum_El55['40GHz'], label='40GHz')
                plt.plot(time_axis, T_spectrum_El55['40GHz'], c='k', alpha=0.5)
                plt.scatter(time_axis, T_spectrum_El55['90GHz'], label='90GHz')
                plt.plot(time_axis, T_spectrum_El55['90GHz'], c='k', alpha=0.5)
                plt.scatter(time_axis, T_spectrum_El55['150GHz'], label='150GHz')
                plt.plot(time_axis, T_spectrum_El55['150GHz'], c='k', alpha=0.5)
                plt.scatter(time_axis, T_spectrum_El55['220GHz'], label='220GHz')
                plt.plot(time_axis, T_spectrum_El55['220GHz'], c='k', alpha=0.5)
                plt.scatter(time_axis, T_spectrum_El55['270GHz'], label='270GHz')
                plt.plot(time_axis, T_spectrum_El55['270GHz'], c='k', alpha=0.5)
                plt.xlabel('UTC time')
                plt.ylabel('Trj[K]')
                plt.suptitle('Trj Time Fluctuations at El=55deg')
                plt.title(day)
                plt.legend()
                plt.savefig(path+day+'_bk_spectrum'+'/T_spectrum_El55.png')
                if remake_post_fig==1:
                    plt.savefig(posting_folder+'T_spectrum_El55_'+day+'.png')
                plt.close()

            if rewrite_txt==1:
                write_txt()
                read_txt()
            else:
                read_txt()



    def plot_Trj_skydip(self, wvr_scan, remake_post_fig=1, rewrite_txt=1, datafolder='', posting_folder='/Volumes/Data Analysis/BICEP/Postings/WVR_postings/20210421_TrjSpectralExtrapolation/plots/'):
        path='am_datfiles/SPole_annual_50/'+wvr_scan[:-4]+'/'

        def write_txt():

            T_spectrum_El55={'day':[], 'time':[], '30GHz':[], '40GHz':[], '90GHz':[], '150GHz':[], '220GHz':[], '270GHz':[]}

            dat_exist=1

            pk_pwv=path+wvr_scan[:-4]+'_fitoutput_corr.txt'
            # t = open(pk_pwv, 'rb')
            # fit_output = pickle.load(t)
            # t.close()

            if not os.path.exists(path+'/spectra/plots'):
                os.makedirs(path+'/spectra/plots')

            test_data=wvr_scan


            if not os.path.exists(pk_pwv):
                try:
                    self.create_am_datfile(test_data, path_to_data='wvr1_data/', spline=2)
                    self.fit_w_am(test_data, path_to_data='wvr1_data/', spline=2)
                except Exception as e:
                    print('dat file '+test_data+' not found.')
                    print(e)
                    dat_exist=0

            if dat_exist==1:
                T_spectrum_El55['day'].append(wvr_scan[:8])
                T_spectrum_El55['time'].append(wvr_scan[9:15])

                Trj_dict={'fn': [], 'el':[], 'Trj':[]}
                for amcfile in os.listdir(path):
                    if (amcfile[-4:]=='.amc') & (amcfile[-12:-8]=='corr'):
                        el=amcfile[30:32]

                        a_file = open(path+'/'+amcfile, "r")
                        list_of_lines = a_file.readlines()
                        print(list_of_lines[27])
                        list_of_lines[27] = "f 10 GHz  350 GHz  200 MHz\n"
                        list_of_lines[28]=''
                        list_of_lines[29]=''
                        print(list_of_lines[27])
                        b_file = open(path+'/spectra/'+amcfile[:-4]+'_bk.amc', "w")
                        b_file.writelines(list_of_lines)
                        b_file.close()

                        newamcfile=amcfile[:-4]+'_bk.amc'

                        os.system('am '+path+'/spectra/'+newamcfile+' >'+path+'/spectra/'+newamcfile[:-4]+'.out')

                        outfile=path+'/spectra/'+newamcfile[:-4]+'.out'
                        f=open(outfile,"r")
                        lines=f.readlines()
                        print(lines[0])

                        freq=[]
                        tau=[]
                        Trj=[]
                        for x in lines:
                            freq.append(float(x.split(' ')[0]))
                            tau.append(float(x.split(' ')[1]))
                            Trj.append(float(x.split(' ')[2]))
                        f.close()

                        Trj_f=np.zeros(6)

                        fig, ax = pl.subplots(6,1)
                        avg30=np.mean(np.array(Trj)[np.where((np.array(freq)>25)&(np.array(freq)<35))])
                        ax[0].plot(freq, Trj, label='Trj[30GHz]='+str(round(avg30,2)))
                        Trj_f[0]=round(avg30,2)
                        ax[0].set_ylabel('T_rj[K]')
                        ax[0].set_xlim(25,35)
                        int0=np.where((np.array(freq)>25)&(np.array(freq)<35))
                        print(np.array(Trj))
                        ax[0].set_ylim(np.min(np.array(Trj)[int0[0]]), np.max(np.array(Trj)[int0[0]]))
                        ax[0].legend()

                        avg40=np.mean(np.array(Trj)[np.where((np.array(freq)>35)&(np.array(freq)<=45))])
                        ax[1].plot(freq, Trj, label='Trj[40GHz]='+str(round(avg40,2)))
                        Trj_f[1]=round(avg40,2)
                        ax[1].set_ylabel('T_rj[K]')
                        ax[1].set_xlim(35,45)
                        int1=np.where((np.array(freq)>35)&(np.array(freq)<45))
                        print(np.array(Trj)[int1])
                        ax[1].set_ylim(np.min(np.array(Trj)[int1[0]]), np.max(np.array(Trj)[int1[0]]))
                        ax[1].legend()

                        avg90=np.mean(np.array(Trj)[np.where((np.array(freq)>80.9)&(np.array(freq)<107.5))]) #BC=94.2; BW=26.7
                        ax[2].plot(freq, Trj, label='Trj[90GHz]='+str(round(avg90,2)))
                        Trj_f[2]=round(avg90,2)
                        ax[2].set_ylabel('T_rj[K]')
                        ax[2].set_xlim(80,108)
                        int2=np.where((np.array(freq)>80.9)&(np.array(freq)<107.5))
                        print(np.array(Trj)[int2])
                        ax[2].set_ylim(np.min(np.array(Trj)[int2[0]]), np.max(np.array(Trj)[int2[0]]))
                        ax[2].legend()

                        avg150=np.mean(np.array(Trj)[np.where((np.array(freq)>127.3)&(np.array(freq)<170.7))]) #BC=149;BW=43.4
                        ax[3].plot(freq, Trj, label='Trj[150GHz]='+str(round(avg150,2)))
                        Trj_f[3]=round(avg150,2)
                        ax[3].set_xlabel('freq[GHz]')
                        ax[3].set_ylabel('T_rj[K]')
                        ax[3].set_xlim(127, 171)
                        int3=np.where((np.array(freq)>127.3)&(np.array(freq)<170.7))
                        print(np.array(Trj)[int3])
                        ax[3].set_ylim(np.min(np.array(Trj)[int3[0]]), np.max(np.array(Trj)[int3[0]]))
                        ax[3].legend()

                        avg220=np.mean(np.array(Trj)[np.where((np.array(freq)>205.6)&(np.array(freq)<256.8))]) #BC=231.2; BW=51.2
                        ax[4].plot(freq, Trj, label='Trj[220GHz]='+str(round(avg220,2)))
                        Trj_f[4]=round(avg220,2)
                        ax[4].set_ylabel('T_rj[K]')
                        ax[4].set_xlim(205,257)
                        int4=np.where((np.array(freq)>205.6)&(np.array(freq)<256.8))
                        print(np.array(Trj)[int4])
                        ax[4].set_ylim(np.min(np.array(Trj)[int4[0]]), np.max(np.array(Trj)[int4[0]]))
                        ax[4].legend()

                        avg270=np.mean(np.array(Trj)[np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))])#BC=275.5;BW=69.8
                        ax[5].plot(freq, Trj, label='Trj[270GHz]='+str(round(avg270,2)))
                        Trj_f[5]=round(avg270,2)
                        ax[5].set_xlabel('freq[GHz]')
                        ax[5].set_ylabel('T_rj[K]')
                        ax[5].set_xlim(240,311)
                        int5=np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))
                        print(np.array(Trj)[int5])
                        ax[5].set_ylim(np.min(np.array(Trj)[int5[0]]), np.max(np.array(Trj)[int5[0]]))
                        ax[5].legend()

                        if not os.path.exists(path+'/spectra/plots/'):
                            os.makedirs(path+'/spectra/plots/')
                        # if remake_post_fig==1:
                        #     plt.savefig(posting_path+'spectra/TrjinBKbands_'+newamcfile[:-4]+'.png')
                        pl.title('El='+str(el))
                        plt.close()

                        plt.plot(freq, Trj)
                        plt.ylabel('T_rj[K]')
                        plt.xlabel('Freq[GHz]')

                        plt.title(newamcfile[:-4])
                        # if remake_post_fig==1:
                        #     plt.savefig(posting_path+'spectra/TrjinBKbands_fullspectrum_'+newamcfile[:-4]+'.png')
                        plt.savefig(path+'/spectra/plots/'+newamcfile[:-4]+'_fullspectrum.png')
                        plt.close()

                        Trj_dict['fn'].append(newamcfile[:-4]+'_Trj_pk.txt')
                        Trj_dict['el'].append(float(newamcfile[30:34]))
                        Trj_dict['Trj'].append(Trj_f)

                        pl.plot(Trj_f)
                        pl.title('el='+str(el))
                        pl.close()


                f=open(path+'Trj_el_'+wvr_scan[:-9]+'_pk.txt',"wb")
                pickle.dump(Trj_dict, f)
                f.close()

                pl.plot(Trj_dict['el'], Trj_dict['Trj'][0])
                pl.show()

                return Trj_dict



        def read_txt():
            T_spectrum_El55={'day':[], 'time':[], '30GHz':[], '40GHz':[], '90GHz':[], '150GHz':[], '220GHz':[], '270GHz':[]}

            f=open(path+'Trj_el_'+wvr_scan[:-9]+'_pk.txt',"rb")
            Trj_dict=pickle.load(f)
            f.close()
            #
            # print(Trj_dict)
            # print(Trj_dict['el'][0])

            El_list=[]
            T30=[]
            T40=[]
            T90=[]
            T150=[]
            T210=[]
            T270=[]

            for i in range (len(Trj_dict['fn'])):
                El_list.append(Trj_dict['el'][i])
                T=Trj_dict['Trj'][i]
                print(T)
                T30.append(T[0])
                T40.append(T[1])
                T90.append(T[2])
                T150.append(T[3])
                T210.append(T[4])
                T270.append(T[5])

            # fig, ax = pl.subplots(2,1, sharex=True)
            #
            # ax[0].plot(El_list, T30, label='30GHz')
            # ax[0].plot(El_list, T40, label='40GHz')
            # ax[0].plot(El_list, T90, label='90GHz')
            # ax[0].plot(El_list, T150, label='150GHz')
            # ax[0].plot(El_list, T210, label='220GHz')
            # ax[0].plot(El_list, T270, label='270GHz')
            # # ax[0].xlabel('El[deg]')
            # ax[0].set_ylabel('T_rj[K]')
            # ax[0].legend()
            #
            # ax[1].scatter(fit_output['El'], fit_output['pwv'], s=5)
            # ax[1].plot(fit_output['El'], fit_output['pwv'], c='k', alpha=0.5)
            # ax[1].set_xlabel('El[deg]')
            # ax[1].set_ylabel('PWV[um]')
            # ax[1].set_xlim(55.,55.9)
            #
            # plt.suptitle(wvr_scan)
            # if remake_post_fig==1:
            #     plt.savefig(posting_path+'TrjandPWV_'+wvr_scan+'.png')
            # plt.savefig(path+'/'+wvr_scan+'.png')
            # #plt.show()
            # plt.close()

            T_spectrum_El55['30GHz'].append(T30[0])
            T_spectrum_El55['40GHz'].append(T40[0])
            T_spectrum_El55['90GHz'].append(T90[0])
            T_spectrum_El55['150GHz'].append(T150[0])
            T_spectrum_El55['220GHz'].append(T210[0])
            T_spectrum_El55['270GHz'].append(T270[0])

            if not os.path.exists(path+wvr_scan+'_bk_spectrum'):
                os.makedirs(path+wvr_scan+'_bk_spectrum')

            f=open(path+wvr_scan+'_bk_spectrum'+'/T_spectrum_El55.txt',"wb")
            pickle.dump(T_spectrum_El55, f)
            f.close()

            f=open(path+wvr_scan+'_bk_spectrum'+'/T_spectrum_El55.txt',"rb")
            T_spectrum_El55=pickle.load(f)
            f.close()

            time_axis=[]
            for j in range (len(T_spectrum_El55['time'])):
                t=T_spectrum_El55['time'][j]
                time_axis.append(t[:2])

            print('time_axis=', time_axis)
            print('shape=', np.shape(time_axis))
            print('T_spectrum_El55[\'30GHz\']=', T_spectrum_El55['30GHz'])
            print('shape=', np.shape(T_spectrum_El55['30GHz']))

            plt.scatter(time_axis, T_spectrum_El55['30GHz'], label='30GHz')
            plt.plot(time_axis, T_spectrum_El55['30GHz'], c='k', alpha=0.5)
            plt.scatter(time_axis, T_spectrum_El55['40GHz'], label='40GHz')
            plt.plot(time_axis, T_spectrum_El55['40GHz'], c='k', alpha=0.5)
            plt.scatter(time_axis, T_spectrum_El55['90GHz'], label='90GHz')
            plt.plot(time_axis, T_spectrum_El55['90GHz'], c='k', alpha=0.5)
            plt.scatter(time_axis, T_spectrum_El55['150GHz'], label='150GHz')
            plt.plot(time_axis, T_spectrum_El55['150GHz'], c='k', alpha=0.5)
            plt.scatter(time_axis, T_spectrum_El55['220GHz'], label='220GHz')
            plt.plot(time_axis, T_spectrum_El55['220GHz'], c='k', alpha=0.5)
            plt.scatter(time_axis, T_spectrum_El55['270GHz'], label='270GHz')
            plt.plot(time_axis, T_spectrum_El55['270GHz'], c='k', alpha=0.5)
            plt.xlabel('UTC time')
            plt.ylabel('Trj[K]')
            plt.suptitle('Trj Time Fluctuations at El=55deg')
            plt.title(day)
            plt.legend()
            plt.savefig(path+day+'_bk_spectrum'+'/T_spectrum_El55.png')
            # if remake_post_fig==1:
            #     plt.savefig(posting_path+'T_spectrum_El55_'+day+'.png')
            plt.show()

        if rewrite_txt==1:
            write_txt()
            read_txt()
        else:
            read_txt()

    def plot_Trj_az(self, scan, Az_min=0, Az_max=360, remake_post_fig=1, rewrite_txt=0, rewrite_amc=0, out_pk='Trj_pk_BAK.txt', datafolder='BAElnod_data/', posting_folder='/Volumes/LaCie/BICEP/Postings/WVR_postings/20210421_TrjSpectralExtrapolation/plots/'):

        def plot_atmo(waz, pwv_list, title='', show=0):
            az, fs =  raz.findScans(waz)
            D_pwv = raz.interpToImage(az, pwv_list, fs)

            fig=plt.figure()
            im=plt.imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
            plt.suptitle('PWV Atmogram\n'+title)
            cbar = fig.colorbar(im, extend='both')
            cbar.set_label('PWV[um]')
            if show==0:
                pl.show()
            else:
                pl.close()

            return az, fs, D_pwv


        def Trjatmo_from_pwvatmo(fn):

            if not os.path.exists('am_datfiles_Az/SPole_annual_50/'+fn[:-4]+'/'):
                os.makedirs('am_datfiles_Az/SPole_annual_50/'+fn[:-4]+'/')

            pickle_fn=self.path_to_all+'am_datfiles_Az/SPole_annual_50/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'
            pickle_fn_temps=self.path_to_all+'am_datfiles_Az/SPole_annual_50/'+fn[:-4]+'/'+fn[:-4]+'_clean_mod3_method_import_model_pickle_temps.txt'
            #using fit or import_model should be the same here as I am just using it to extract waz

            f = open(pickle_fn_temps,'rb')
            pickle_Temps = pickle.load(f)
            f.close()
            waz=pickle_Temps['wAz']

            f = open(pickle_fn,'rb')
            fit_output = pickle.load(f)
            f.close()
            pwv_ts=fit_output['pwv']

            az, fs, D_pwv = plot_atmo(waz, pwv_ts)

            return az, fs, pwv_ts, D_pwv

        az, fs, pwv, D_pwv = Trjatmo_from_pwvatmo(scan)

        path=self.path_to_all+'am_datfiles_Az/SPole_annual_50/'
        #path='../../pipeline/am_datfiles_Az/SPole_annual_50/'

        folder=path+scan[:-4]


        def write_txt():
            Trj_dict={'fn': [], 'scanN':[], 'Az':[], 'waz':[], 'Trj':[]}
            Az_infolder_list=[]
            for flist in os.listdir(folder):
                if flist[-4:]=='.amc':
                    if flist[-16:-11]=='model':
                        amcfile=flist
                        scanN_f_list= r.findall(r'\d+', flist[30:36])
                        scanN_f=scanN_f_list[0]
                        Az_f_list=r.findall(r'\d+', flist[37:47])
                        Az_f=Az_f_list[0]+'.'+Az_f_list[1]
                        waz_f=float(Az_f)+(int(scanN_f)*360.)
                        if ((float(Az_f)<Az_max)&(float(Az_f)>Az_min)):
                            #print(flist)
                            Trj_dict['scanN'].append(float(scanN_f))
                            Trj_dict['Az'].append(float(Az_f))
                            Trj_dict['waz'].append(waz_f)
                            newamcfile=amcfile[:-4]+'_bk.amc'
                            print('scanN, Az, waz =', scanN_f, Az_f, waz_f)

                            os.system('mv '+folder+'/spectra/'+amcfile[:-4]+'_bk.out '+folder+'/spectra/'+amcfile[:-4]+'_bk_old.out')
                            #temporary: I want to remake the am files but can't in one shot (killed 9) and I don't want to restart from zero
                            #every time

                            if rewrite_amc == 0:
                                condition = os.path.exists(folder+'/spectra/'+amcfile[:-4]+'_bk.out')
                            else:
                                condition = False #Always like if amc doesn't exist

                            if not condition:
                                a_file = open(folder+'/'+amcfile, "r")
                                list_of_lines = a_file.readlines()
                                # list_of_lines[27] = "f 10 GHz  300 GHz  200 MHz\n"
                                list_of_lines[27] = "f 10 GHz  350 GHz  200 MHz\n"
                                list_of_lines[28]=''
                                list_of_lines[29]=''
                                if not os.path.exists(folder+'/spectra/'):
                                    os.makedirs(folder+'/spectra/')
                                b_file = open(folder+'/spectra/'+amcfile[:-4]+'_bk.amc', "w")
                                b_file.writelines(list_of_lines)
                                b_file.close()
                                os.system('am '+folder+'/spectra/'+newamcfile+' >'+folder+'/spectra/'+newamcfile[:-4]+'.out')
                            outfile=folder+'/spectra/'+newamcfile[:-4]+'.out'
                            f=open(outfile,"r")
                            lines=f.readlines()

                            freq=[]
                            tau=[]
                            Trj=[]
                            for x in lines:
                                freq.append(float(x.split(' ')[0]))
                                tau.append(float(x.split(' ')[1]))
                                Trj.append(float(x.split(' ')[2]))
                            f.close()

                            Trj = np.array(Trj)
                            Trj_f=np.zeros(6)

                            try:
                                fig, ax = pl.subplots(6,1)
                                avg30=np.mean(np.array(Trj)[np.where((np.array(freq)>25)&(np.array(freq)<35))])
                                ax[0].plot(freq, Trj, label='Trj[30GHz]='+str(round(avg30,2)))
                                Trj_f[0]=avg30
                                ax[0].set_ylabel('T_rj[K]')
                                ax[0].set_xlim(25,35)
                                int0=np.where((np.array(freq)>25)&(np.array(freq)<35))
                                ax[0].set_ylim(np.min(np.array(Trj)[int0[0]]), np.max(np.array(Trj)[int0[0]]))
                                ax[0].legend()

                                avg40=np.mean(np.array(Trj)[np.where((np.array(freq)>35)&(np.array(freq)<=45))])
                                ax[1].plot(freq, Trj, label='Trj[40GHz]='+str(round(avg40,2)))
                                Trj_f[1]=avg40
                                ax[1].set_ylabel('T_rj[K]')
                                ax[1].set_xlim(35,45)
                                int1=np.where((np.array(freq)>35)&(np.array(freq)<45))
                                ax[1].set_ylim(np.min(np.array(Trj)[int1[0]]), np.max(np.array(Trj)[int1[0]]))
                                ax[1].legend()

                                avg90=np.mean(np.array(Trj)[np.where((np.array(freq)>80.9)&(np.array(freq)<107.5))]) #BC=94.2; BW=26.7
                                ax[2].plot(freq, Trj, label='Trj[90GHz]='+str(round(avg90,2)))
                                Trj_f[2]=avg90
                                ax[2].set_ylabel('T_rj[K]')
                                ax[2].set_xlim(80,108)
                                int2=np.where((np.array(freq)>80.9)&(np.array(freq)<107.5))
                                ax[2].set_ylim(np.min(np.array(Trj)[int2[0]]), np.max(np.array(Trj)[int2[0]]))
                                ax[2].legend()

                                avg150=np.mean(np.array(Trj)[np.where((np.array(freq)>127.3)&(np.array(freq)<170.7))]) #BC=149;BW=43.4
                                ax[3].plot(freq, Trj, label='Trj[150GHz]='+str(round(avg150,2)))
                                Trj_f[3]=avg150
                                ax[3].set_xlabel('freq[GHz]')
                                ax[3].set_ylabel('T_rj[K]')
                                ax[3].set_xlim(127, 171)
                                int3=np.where((np.array(freq)>127.3)&(np.array(freq)<170.7))
                                ax[3].set_ylim(np.min(np.array(Trj)[int3[0]]), np.max(np.array(Trj)[int3[0]]))
                                ax[3].legend()


                                passband_k210 = bk.keck_bandpasses(bp_file='K210_frequency_spectrum_20160101.txt')
                                freq_pb, response_unform, response_rj = passband_k210
                                f = sint.interp1d(freq_pb, response_rj)
                                freq = np.array(freq)
                                mask_210 = np.where((np.array(freq)>np.min(freq_pb))&(np.array(freq)<np.max(freq_pb)))
                                freq_210 = freq[mask_210]
                                band_210 = f(freq_210)

                                avg210=np.nanmean(np.array(Trj[mask_210])*band_210) #BC=231.2; BW=51.2
                                ax[4].plot(freq, Trj, label='Trj[210GHz]='+str(round(avg210,2)))
                                Trj_f[4]=avg210
                                ax[4].set_ylabel('T_rj[K]')
                                ax[4].set_xlim(205,257)
                                #int4=np.where((np.array(freq)>205.6)&(np.array(freq)<256.8))
                                ax[4].set_ylim(np.min(Trj[mask_210]), np.max(Trj[mask_210]))
                                ax[4].legend()

                                passband_k270 = bk.keck_bandpasses(bp_file='K270_frequency_spectrum_20170710.txt')
                                freq_pb, response_unform, response_rj = passband_k270
                                print('max_freq_270=',np.max(freq_pb))
                                f = sint.interp1d(freq_pb, response_rj)
                                freq = np.array(freq)
                                mask_270 = np.where((np.array(freq)>np.min(freq_pb))&(np.array(freq)<np.max(freq_pb)))
                                freq_270 = freq[mask_270]
                                band_270 = f(freq_270)
                                #avg270=np.mean(np.array(Trj)[np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))])#BC=275.5;BW=69.8
                                Trj_270 = np.array(Trj[mask_270])*band_270

                                avg270=np.nanmean(np.array(Trj[mask_270])*band_270)#BC=275.5;BW=69.8
                                ax[5].plot(freq, Trj, label='Trj[270GHz]='+str(round(avg270,2)))
                                Trj_f[5]=avg270
                                ax[5].set_xlabel('freq[GHz]')
                                ax[5].set_ylabel('T_rj[K]')
                                ax[5].set_xlim(240,311)
                                #int5=np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))
                                ax[5].set_ylim(np.min(Trj[mask_270]), np.max(Trj[mask_270]))
                                ax[5].legend()
                                plt.savefig(folder+'/spectra/plots/'+newamcfile[:-4]+'.png')
                                pl.close()


                                pl.figure()
                                mean_pwv = round(np.nanmean(D_pwv), 2)
                                #pl.plot(freq[np.where((np.array(freq)>240.6)&(np.array(freq)<310.4))], band_270, label='K270 bandpass')
                                pl.plot(freq_270, band_270*np.nanmax(Trj[mask_270]))
                                pl.fill_between(freq_270, band_270*np.nanmax(Trj[mask_270]), alpha=0.5, label='Keck270')
                                pl.plot(freq_210, band_210*np.nanmax(Trj[mask_210]))
                                pl.fill_between(freq_210, band_210*np.nanmax(Trj[mask_210]), alpha=0.5, label='Keck210')
                                pl.plot(freq, Trj, label='Trj', c='k')
                                pl.xlabel('freq[GHz]')
                                pl.ylabel('T_rj[K]')
                                pl.xlim(100, 350)
                                pl.suptitle('Data Set : '+scan[:-16])
                                pl.title('Nscan = '+scanN_f+' - Az = '+ Az_f)
                                pl.legend()
                                plt.savefig(folder+'/spectra/plots/'+newamcfile[:-4]+'_K210_K270_bandpasses.png')
                                pl.close()

                                if not os.path.exists(folder+'/spectra/plots/'):
                                    os.makedirs(folder+'/spectra/plots/')

                                Trj_dict['fn'].append(newamcfile)
                                Trj_dict['Trj'].append(Trj_f)

                            except:
                                Trj_dict['Trj'].append(np.full((6,), np.nan))


            f=open(folder+'/spectra/'+out_pk,"wb")
            pickle.dump(Trj_dict, f)
            f.close()

        def read_txt():
            print('path to pk =', folder+'/spectra/'+out_pk)
            if os.path.exists(folder+'/spectra/'+out_pk):
                f=open(folder+'/spectra/'+out_pk,"rb")
                Trj_dict=pickle.load(f)
                f.close()
            else:
                pk_input = input(out_pk+" does not exist. Which pk file do you want to use?\n")
                f=open(folder+'/spectra/'+pk_input,"rb")
                Trj_dict=pickle.load(f)
                f.close()

            El_list=[]
            T30=[]
            T40=[]
            T90=[]
            T150=[]
            T210=[]
            T270=[]
            Az_list=[]
            ScanN_list=[]
            waz_list=[]


            print(Trj_dict)

            for i in range (len(Trj_dict['fn'])):
                Az_list.append(Trj_dict['Az'][i])
                ScanN_list.append(int(Trj_dict['scanN'][i]))
                waz_list.append(Trj_dict['waz'][i])

                T=Trj_dict['Trj'][i]

                try:
                    T30.append(T[0])
                    T40.append(T[1])
                    T90.append(T[2])
                    T150.append(T[3])
                    T210.append(T[4])
                    T270.append(T[5])
                except:
                    T30.append(np.nan)
                    T40.append(np.nan)
                    T90.append(np.nan)
                    T150.append(np.nan)
                    T210.append(np.nan)
                    T270.append(np.nan)


            waz_list=np.array(waz_list)
            Az_list=np.array(Az_list)
            T30=np.array(T30)
            T40=np.array(T40)
            T90=np.array(T90)
            T150=np.array(T150)
            T210=np.array(T210)
            T270=np.array(T270)

            idx=np.argsort(waz_list)
            waz_list=waz_list[idx]
            Az_list=Az_list[idx]
            T30=T30[idx]
            T40=T40[idx]
            T90=T90[idx]
            T150=T150[idx]
            T210=T210[idx]
            T270=T270[idx]

            plt.figure()
            ax1 = plt.subplot(611)
            ax1.scatter(waz_list, T30, s=3)
            ax1.set_title('T30')
            ax2 = plt.subplot(612)
            ax2.scatter(waz_list, T40, s=3)
            ax2.set_title('T40')
            ax3 = plt.subplot(613)
            ax3.scatter(waz_list, T90, s=3)
            ax3.set_title('T90')
            ax4 = plt.subplot(614)
            ax4.scatter(waz_list, T150, s=3)
            ax4.set_title('T150')
            ax5 = plt.subplot(615)
            ax5.scatter(waz_list, T210, s=3)
            ax5.set_title('T210')
            ax6 = plt.subplot(616)
            ax6.scatter(waz_list, T270, s=3)
            ax6.set_title('T270')
            plt.suptitle('Trj fluctuations')
            ax6.set_xlabel('wAz')
            plt.savefig(posting_folder+'Trj_fluctuations.png')
            #plt.show()
            plt.close()

            # plt.figure()
            # ax1 = plt.subplot(211)
            # ax1.scatter(Az_list, T210, s=3, c='r')
            # ax1.plot(Az_list, T210 , c='k', alpha=0.5)
            # ax1.set_title('Az')
            # ax2 = plt.subplot(212)
            # ax2.scatter(waz_list, T210, s=3, c='r')
            # ax2.plot(waz_list, T210, c='k', alpha=0.5)
            # ax2.set_title('Waz')
            # plt.show()
            #
            # plt.scatter(Az_list, waz_list, s=3, c='r')
            # plt.plot(Az_list, waz_list, c='k', alpha=0.5)
            #
            # plt.show()


            newaz,newfs =  raz.findScans(waz_list)

            plt.figure()
            ax1 = plt.subplot(611)
            ax1.scatter(newaz, T30, c='r', s=2)
            ax1.plot(newaz, T30, c='k', alpha=0.5)
            ax1.set_ylabel('T[K]')
            ax1.set_title('T30')
            ax2 = plt.subplot(612)
            ax2.scatter(newaz, T40, c='r', s=2)
            ax2.plot(newaz, T40, c='k', alpha=0.5)
            ax2.set_ylabel('T[K]')
            ax2.set_title('T40')
            ax3 = plt.subplot(613)
            ax3.scatter(newaz, T90, c='r', s=2)
            ax3.plot(newaz, T90, c='k', alpha=0.5)
            ax3.set_ylabel('T[K]')
            ax3.set_title('T90')
            ax4 = plt.subplot(614)
            ax4.scatter(newaz, T150, c='r', s=2)
            ax4.plot(newaz, T150, c='k', alpha=0.5)
            ax4.set_ylabel('T[K]')
            ax4.set_title('T150')
            ax5 = plt.subplot(615)
            ax5.scatter(newaz, T210, c='r', s=2)
            ax5.plot(newaz, T210, c='k', alpha=0.5)
            ax5.set_ylabel('T[K]')
            ax5.set_title('T210')
            ax6 = plt.subplot(616)
            ax6.scatter(newaz, T270, c='r', s=2)
            ax6.plot(newaz, T270, c='k', alpha=0.5)
            ax6.set_ylabel('T[K]')
            ax6.set_title('T270')
            plt.suptitle('Trj fluctuations')
            ax6.set_xlabel('Az')
            plt.savefig(posting_folder+'Trj_fluctuations_Az_Subplots_WVR_Scan'+scan[:-4]+'.png')
            #plt.show()
            plt.close()




            D_30GHz= raz.interpToImage(newaz, T30, newfs)
            D_40GHz= raz.interpToImage(newaz, T40, newfs)
            D_90GHz= raz.interpToImage(newaz, T90, newfs)
            D_150GHz= raz.interpToImage(newaz, T150, newfs)
            D_210GHz= raz.interpToImage(newaz, T210, newfs)
            D_270GHz= raz.interpToImage(newaz, T270, newfs)

            figsize=(16,12)
            fs=18
            fs_ticks=16

            fig=plt.figure(figsize=figsize)
            ax1 = plt.subplot(111)
            cax1 = plt.imshow(D_30GHz,aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)

            #plt.title('Trj_30GHz[K]')
            cbar = fig.colorbar(cax1)
            cbar.set_label('30GHz', fontsize=fs)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)
            #plt.colorbar()

            title = 'T_rj[K]\n'+scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)
            plt.savefig(posting_folder+'Trj_30_atmo_clean_mod3_clean_method_monthly_tilt_model_tmp_zoom_off_WVR_Scan'+scan[:-4]+'.png')
            #plt.show()
            plt.close()

            fig=plt.figure(figsize=figsize)
            ax2 = plt.subplot(111, sharex=ax1)
            #cax2 = plt.imshow(D_40GHz[:,:66], aspect='auto',interpolation='nearest', origin='lower')
            cax2 = plt.imshow(D_40GHz, aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax2.get_xticklabels(), visible=False)
            ax2.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            cbar = fig.colorbar(cax2)
            cbar.set_label('40GHz', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)
            title = 'T_rj[K]\n'+scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)
            plt.savefig(posting_folder+'Trj_40_atmo_clean_mod3_clean_method_monthly_tilt_model_tmp_zoom_off_WVR_Scan_'+scan[:-4]+'.png')
            #plt.show()
            plt.close()

            fig=plt.figure(figsize=figsize)
            ax3 = plt.subplot(111, sharex=ax1)
            cax3=plt.imshow(D_90GHz, aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax3.get_xticklabels(), visible=False)
            ax3.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            #plt.title('Trj_90GHz[K]')
            #plt.colorbar()
            cbar = fig.colorbar(cax3)
            cbar.set_label('90GHz', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)
            title = 'T_rj[K]\n'+scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)
            plt.savefig(posting_folder+'Trj_90_atmo_clean_mod3_clean_method_monthly_tilt_model_tmp_zoom_off_WVR_Scan_'+scan[:-4]+'.png')
            #plt.show()
            plt.close()

            fig=plt.figure(figsize=figsize)
            ax4 = plt.subplot(111, sharex=ax1)
            cax4=plt.imshow(D_150GHz, aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax4.get_xticklabels(), visible=False)
            ax4.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)
            #plt.title('Trj_150GHz[K]')
            #plt.colorbar()
            cbar = fig.colorbar(cax4)
            cbar.set_label('150GHz', fontsize=fs)
            title = 'T_rj[K]\n'+scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)
            plt.savefig(posting_folder+'Trj_150_atmo_clean_mod3_clean_method_monthly_tilt_model_tmp_zoom_off_WVR_Scan_'+scan[:-4]+'.png')
            #plt.show()
            plt.close()

            fig=plt.figure(figsize=figsize)
            ax5 = plt.subplot(111, sharex=ax1)
            cax5=plt.imshow(D_210GHz, aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax5.get_xticklabels(), visible=False)
            ax5.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            #plt.title('Trj_210GHz[K]')
            #plt.colorbar()
            cbar = fig.colorbar(cax5)
            cbar.set_label('210GHz', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)
            title = 'T_rj[K]\n'+scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)
            plt.savefig(posting_folder+'Trj_210_atmo_clean_mod3_clean_method_monthly_tilt_model_tmp_zoom_off_WVR_Scan_'+scan[:-4]+'.png')
            #plt.show()
            plt.close()

            fig=plt.figure(figsize=figsize)
            ax6 = plt.subplot(111, sharex=ax1)
            cax6=plt.imshow(D_270GHz, aspect='auto',interpolation='nearest', origin='lower')
            ax6.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            plt.xlabel('Scan_n', fontsize=fs)
            #plt.title('Trj_270GHz[K]')
            #plt.colorbar()
            cbar = fig.colorbar(cax6)
            cbar.set_label('270GHz', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)

            title = 'T_rj[K]\n'+scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)
            #plt.savefig(folder+'/'+title+'.png')
            #plt.savefig('../Postings/WVR_postings/20210421_TrjSpectralExtrapolation/plots/'+scan[:-4]+'_Trjatmo_clean_mod3_clean_method_monthly_tilt_model_tmp.png')
            plt.savefig(posting_folder+'Trj_270_atmo_clean_mod3_clean_method_monthly_tilt_model_tmp_zoom_off_WVR_Scan_'+scan[:-4]+'.png')
            #plt.show()
            plt.close()


            plt.figure(figsize=(14,12))
            fs= 18
            ax1 = plt.subplot(311)

            cax1=plt.imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax5.get_xticklabels(), visible=False)
            ax1.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            plt.title('PWV', fontsize=16)
            plt.colorbar(label='um')
            cbar = fig.colorbar(cax1)
            cbar.set_label('um', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            plt.xticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)

            ax5 = plt.subplot(312, sharex=ax1)

            cax5=plt.imshow(D_210GHz, aspect='auto',interpolation='nearest', origin='lower')
            #plt.setp(ax5.get_xticklabels(), visible=False)
            ax5.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            plt.title('Trj_210GHz', fontsize=16)
            plt.colorbar(label='K')
            cbar = fig.colorbar(cax5)
            cbar.set_label('K', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            plt.xticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)

            ax6 = plt.subplot(313, sharex=ax1)

            cax6=plt.imshow(D_270GHz, aspect='auto',interpolation='nearest', origin='lower')
            ax6.set_ylim(Az_min, Az_max)
            plt.ylabel('Az [deg]', fontsize=fs)
            plt.xlabel('Scan_n', fontsize=fs)
            plt.title('Trj_270GHz', fontsize=16)
            plt.colorbar(label='K')
            cbar = fig.colorbar(cax6)
            cbar.set_label('K', fontsize=fs)
            plt.yticks(fontsize=fs_ticks)
            plt.xticks(fontsize=fs_ticks)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fs_ticks)

            title = scan.replace('.txt','')
            plt.suptitle(title,y=0.97, fontsize=20)

            #plt.show()
            plt.close()

            return D_30GHz, D_40GHz, D_90GHz, D_150GHz, D_210GHz, D_270GHz, newaz


        if rewrite_txt==1:
            print('Rewriting txt.')
            write_txt()
            D_30GHz, D_40GHz, D_90GHz, D_150GHz, D_210GHz, D_270GHz, newaz = read_txt()
        else:
            D_30GHz, D_40GHz, D_90GHz, D_150GHz, D_210GHz, D_270GHz, newaz = read_txt()

        return D_30GHz, D_40GHz, D_90GHz, D_150GHz, D_210GHz, D_270GHz, newaz
