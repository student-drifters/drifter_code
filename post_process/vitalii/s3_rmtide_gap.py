"""
s3_rmtide_gap.py

remove tide and find gaps in data

compare drifter dataset velocity with fvcom hincast
following the drifter trajectory.
fvcom data pydap access from monthly files.

To do:
1. spacial interpolation of fvcom velocity data
currently, the nearest neighbor
@author: Vitalii Sheremet, FATE Project
"""

# -*- coding: utf-8 -*-
import numpy as np
#from pydap.client import open_url
import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
from scipy import interpolate
import sys
from SeaHorseTide import *


SOURCEDIR='driftfvcom_data2/'
DESTINDIR='driftfvcom_data3/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])
kdr=0


while kdr in range(len(FNs)):
#while kdr in range(0,284):
#while kdr in range(9,10):

    FN=FNs[kdr]
    FN1=SOURCEDIR+FN
    Z=np.load(FN1)
    tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
    udh=Z['udh'];vdh=Z['vdh'];
    umoh=Z['umoh'];vmoh=Z['vmoh'];
    uwsh=Z['uwsh'];vwsh=Z['vwsh'];
    Z.close()

    print kdr, FN

# need the raw drifter tracks to find gaps
    FND=FN[0:-4]
    FND='drift_data/'+FND+'.csv'
#driftfvcom_data/ID_19965381_dm.npz
    D = np.genfromtxt(FND,dtype=None,names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH_I'],delimiter=',')

    td=np.array(D['TimeRD'])
    latd=np.array(D['LAT_DD'])
    lond=np.array(D['LON_DD'])

    #start and end times close to a whole hour
    t1=np.ceil(np.min(td)*24.)/24.
    t2=np.floor(np.max(td)*24.)/24.

# remove tidal signal
    udm=sh_rmtide(udh,ends=np.NaN)
    vdm=sh_rmtide(vdh,ends=np.NaN)
    udti=udh-udm
    vdti=vdh-vdm

    umom=sh_rmtide(umoh,ends=np.NaN)
    vmom=sh_rmtide(vmoh,ends=np.NaN)
    umoti=umoh-umom
    vmoti=vmoh-vmom

# it does not make much sense to remove tides from the wind stress
# but it might be helpful if comparing with the rmtide filtered velocities
    uwsm=sh_rmtide(uwsh,ends=np.NaN)
    vwsm=sh_rmtide(vwsh,ends=np.NaN)


# find gaps
    tgap=tdh*0.+1./24.
    for i in range(len(tdh)):
        ti=tdh[i]
        i1=max(np.argwhere(td<=ti))
        i2=min(np.argwhere(td>ti))
        tgap[i]=td[i2]-td[i1]

    flag_gap=tdh*0.+1.
    for i in range(len(tdh)):
        if (tgap[i] > 9./24.):
            flag_gap[i]=np.NaN

    flag=umoh*0.+1. # 1 inside GOM3 NaN outside
    flag=flag*flag_gap

    plt.figure(1)
    plt.plot(londh*flag,latdh*flag,'r.-')


    """
    plt.figure()
    plt.plot(tdhgap,udh,'b.-',tdhgap,umoh,'r.-',tdhgap,udm,'k-',tdhgap,umom,'g-')
    plt.title('U '+FN0)
    plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
    plt.show()

    plt.figure()
    plt.plot(tdhgap,vdh,'b.-',tdhgap,vmoh,'r.-',tdhgap,vdm,'k-',tdhgap,vmom,'g-')
    plt.title('V '+FN0)
    plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
    plt.show()
    """
    FN2=DESTINDIR+FN
    np.savez(FN2,tdh=tdh,tgap=tgap,flag=flag,londh=londh,latdh=latdh,udh=udh,vdh=vdh,udm=udm,vdm=vdm,udti=udti,vdti=vdti,umoh=umoh,vmoh=vmoh,umom=umom,vmom=vmom,umoti=umoti,vmoti=vmoti,uwsh=uwsh,vwsh=vwsh,uwsm=uwsm,vwsm=vwsm)

    kdr=kdr+1

plt.show()
