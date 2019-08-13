"""
s2_select.py

select drifter files based on condition
@author: Vitalii Sheremet, FATE Project
"""
# -*- coding: utf-8 -*-
import numpy as np
import shutil

# use ls *.npz > FList.csv to create the file list

SOURCEDIR='driftfvcom_data1/'
DESTINDIR='driftfvcom_data2/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])
kdr=0

while kdr in range(len(FNs)):
#while kdr in range(0,284):
#while kdr in range(9,10):

    FN=FNs[kdr]
    FN1=SOURCEDIR+FN
#np.savez(FN2,tdh=tdh,londh=londh,latdh=latdh,udh=udh,vdh=vdh,umoh=u0,vmoh=v0)
    Z=np.load(FN1)
    tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
    udh=Z['udh'];vdh=Z['vdh'];
    umoh=Z['umoh'];vmoh=Z['vmoh'];
    Z.close()

    FND='drift_data/'+FN[0:-4]+'.csv'
#    FN='drift_data/ID_96101.csv'
    D = np.genfromtxt(FND,dtype=None,names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH_I'],delimiter=',')
    depd=np.median(np.array(D['DEPTH_I']))
# only surface drifters
    print kdr, FN, depd

    if np.isnan(umoh).all() or tdh[-1]-tdh[0] < 2. or np.abs(depd)>10.: # track must be longer than 2 days:
        pass
    else: # track must be longer than 2 days
        FN2=DESTINDIR+FN
        print FN2
        #        np.savez(FN2,tdh=tdh,londh=londh,latdh=latdh,udh=udh,vdh=vdh,umoh=umoh,vmoh=vmoh)
        shutil.copyfile(FN1,FN2)

    kdr=kdr+1
