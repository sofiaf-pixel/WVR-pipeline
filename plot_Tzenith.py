import plotElnod as pE
import AM_functions as am
import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time
from zipfile import ZipFile
import tarfile
import datetime


x_am=am.AM_functions()
x=pE.ReadElnod()


class plot_Tz(object):
    '''

    '''

    def __init__(self, unit=None, verb=True):

        '''


        '''

    def save_fn(self, data_set, f):
        f=open(f,"wb")
        pk.dump(data_set, f)
        f.close()

    def read_fn(self,f):
        f=open(f,"rb")
        data=pk.load(f)
        f.close()
        return data


    def unzip_data(self, path, filepath):
        if filepath[-7:]=='.tar.gz':

            outpath=path+filepath[:-7]
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            tar = tarfile.open(path+'/'+filepath, "r:gz")
            tar.extractall(outpath)
            tar.close()
            os.system(f"rm "+path+filepath)



    def read_Tz(self, month, path='wvr1_data/'):

        T0=[]
        T1=[]
        T2=[]
        T3=[]
        el=[]


        filelist=[]
        filelistcut=[]
        t=[]
        for file in os.listdir(path):
            if file[16:22]=='skyDip':
                if file[:6]==month:
                    if not os.path.exists(file[:16]+'_skyDip'):
                        self.unzip_data(path, file)
                    filepath=path+file+'/'+file+'_fast.txt'

                    try:

                        FH_fast=x.read_elnod_fast(filepath)
                        filelist.append(filepath)
                        yy=int(file[:4])
                        mm=int(file[4:6])
                        dd=int(file[6:8])
                        hh=int(file[9:11])
                        t.append(datetime.datetime(yy, mm, dd, hh))
                        filelistcut.append(file[:15])
                        T0_file=FH_fast[:,1]
                        T1_file=FH_fast[:,2]
                        T2_file=FH_fast[:,3]
                        T3_file=FH_fast[:,4]
                        el_file=FH_fast[:,0]

                        el_file=el_file[np.where(el_file>89.)]
                        T0_file=T0_file[np.where(el_file>89.)]
                        T1_file=T1_file[np.where(el_file>89.)]
                        T2_file=T2_file[np.where(el_file>89.)]
                        T3_file=T3_file[np.where(el_file>89.)]

                        T0.append(np.nanmean(T0_file))
                        T1.append(np.nanmean(T1_file))
                        T2.append(np.nanmean(T2_file))
                        T3.append(np.nanmean(T3_file))
                        el.append(np.nanmean(el_file))

                    except:
                        print('File '+filepath+' does not exist.')


        Tzen={'date':t, 'T0':T0, 'T1':T1, 'T2':T2, 'T3':T3}
        self.save_fn(Tzen, 'T_zenith_'+month)

        return Tzen


    def plot_Tz(self, fn, start='20230101', end='20230105'):
        #start/end in format '20230101'

        data_set = self.read_fn(fn)

        t=np.array(data_set['date'])
        T0=np.array(data_set['T0'])
        T1=np.array(data_set['T1'])
        T2=np.array(data_set['T2'])
        T3=np.array(data_set['T3'])

        if not (start==0):
            ys=int(start[:4])
            ms=int(start[4:6])
            ds=int(start[6:8])
            start_date=datetime.datetime(ys, ms, ds)
            ye=int(end[:4])
            me=int(end[4:6])
            de=int(end[6:8])
            end_date=datetime.datetime(ye, me, de)

            mask = (t >= start_date) & (t <= end_date)

            t=t[mask]
            T0=T0[mask]
            T1=T1[mask]
            T2=T2[mask]
            T3=T3[mask]


        fig, axes = pl.subplots(4, 1)
        axes[0].scatter(t, (T0-np.nanmean(T0))/np.nanmean(T0), s=3, c='r')
        axes[0].legend(loc='upper right')
        axes[0].set_ylabel('T_Ch0[K]')

        axes[1].scatter(t, (T1-np.nanmean(T1))/np.nanmean(T1), s=3, c='y')
        axes[1].legend(loc='upper right')
        axes[1].set_ylabel('T_Ch1[K]')

        axes[2].scatter(t, (T2-np.nanmean(T2))/np.nanmean(T2), s=3, c='c')
        axes[2].legend(loc='upper right')
        axes[2].set_ylabel('T_Ch2[K]')

        axes[3].scatter(t, (T3-np.nanmean(T3))/np.nanmean(T3), s=3, c='b')
        axes[3].legend(loc='upper right')
        axes[3].set_ylabel('T_Ch3[K]')
        axes[3].set_xlabel('day')

        if not os.path.exists('output'):
            os.makedirs('output')
        pl.savefig('output/')
        pl.suptitle('T_zenith '+fn[-6:])

        pl.show()
