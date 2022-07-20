import os, sys
#import zipfile
from zipfile import ZipFile
import tarfile

path='wvr1_data/BAElnod_data/'
month='201909'

os.system(f"rm "+path+month+"*.txt")

for filelist in os.listdir(path):
    #os.system(f"rm "+path+".DS_store")
    print(filelist)
    if filelist[:6]==month:
        if filelist[-3:]=='.gz':
            print(filelist[8:])
            if filelist[8:]=='.tar.gz':
                outpath=path
            else:
                outpath=path+filelist[:-7]
                if not os.path.exists(outpath):
                    os.makedirs(outpath)
            print(outpath)
            print('unzipping filelist=', filelist)
            tar = tarfile.open(path+'/'+filelist, "r:gz")
            tar.extractall(outpath)
            tar.close()
    os.system(f"rm "+path+".DS_store")



#
# day_str_list_full=[]
# for d in range(1, 30):
#     if d<10:
#         day_str=month+'0'+str(d)
#     else:
#         day_str=month+str(d)
#     day_str_list_full.append(day_str)
#
#
# for day in day_str_list_full:
#     os.system(f"rm "+path+".DS_store")
#     os.system(f"rm "+path+day+"/.DS_store")
#     #for filelist in os.listdir(path+day+'/'):
#     for filelist in os.listdir(path):
#         if filelist[-3:]=='.gz':
#             print('unzipping filelist=', filelist)
#             outpath=path+filelist[:-7]
#             if not os.path.exists(outpath):
#                 os.makedirs(outpath)
#                 print('unzipping filelist=', filelist)
#                 print('outpath=', outpath)
#                 tar = tarfile.open(path+day+'/'+filelist, "r:gz")
#                 tar.extractall(outpath)
#                 tar.close()
#     os.system(f"rm "+path+".DS_store")
#     os.system(f"rm "+path+day+"/.DS_store")
