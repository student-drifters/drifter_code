"""
s5_binavguv_drifter.py
Author: Vitalii

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
#from scipy import interpolate
import sys
from SeaHorseTide import *
import shutil
import matplotlib.mlab as mlab
import matplotlib.cm as cm


def sh_bindata(x, y, z, xbins, ybins):
    """
    Bin irregularly spaced data on a rectangular grid.

    """
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_median=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_std=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_num=np.zeros((len(xbins)-1,len(ybins)-1),dtype=int)
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
#            k=np.where((ix==iix) and (iy==iiy)) # wrong syntax
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
            zb_median[iix-1,iiy-1]=np.median(z[k])
            zb_std[iix-1,iiy-1]=np.std(z[k])
            zb_num[iix-1,iiy-1]=len(z[k])

    return xb,yb,zb_mean,zb_median,zb_std,zb_num
#    return xb,yb,zb_median,zb_num

###########################################################

SOURCEDIR='driftfvcom_data4/'
DESTINDIR='driftfvcom_data4/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])

lath=np.array([])
lonh=np.array([])
th=np.array([])
flagh=np.array([])
u=np.array([])
v=np.array([])
#for k in range(10):
for k in range(len(FNs)):
    FN=FNs[k]
    #ID_19965381.npz
    FN1=SOURCEDIR+FN
    print k, FN1
    Z=np.load(FN1)
    tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
    udh=Z['udh'];vdh=Z['vdh'];
    umoh=Z['umoh'];vmoh=Z['vmoh'];
    tgap=Z['tgap'];flag=Z['flag'];
    udm=Z['udm'];vdm=Z['vdm'];
    udti=Z['udti'];vdti=Z['vdti'];
    umom=Z['umom'];vmom=Z['vmom'];
    umoti=Z['umoti'];vmoti=Z['vmoti'];
    Z.close()

    lath=np.append(lath,latdh)
    lonh=np.append(lonh,londh)
    th=np.append(th,tdh)
    flagh=np.append(flagh,flag)
#    ekem=((udm-umom)*(udm-umom)+(vdm-vmom)*(vdm-vmom))*0.5*flag
    u1=udm*flag
    v1=vdm*flag
    u=np.append(u,u1)
    v=np.append(v,v1)

i=np.argwhere(np.isnan(u)==False).flatten()
u=u[i]
v=v[i]
lath=lath[i]
lonh=lonh[i]
th=th[i]

x=lonh
y=lath

xi = np.arange(-76.,-56.000001,0.5)
yi = np.arange(35.,47.000001,0.25)

xb,yb,ub_mean,ub_median,ub_std,ub_num = sh_bindata(x, y, u, xi, yi)
xb,yb,vb_mean,vb_median,vb_std,vb_num = sh_bindata(x, y, v, xi, yi)
np.savez('binned_drifter.npz',xb=xb,yb=yb,ub_mean=ub_mean,ub_median=ub_median,ub_std=ub_std,ub_num=ub_num,vb_mean=vb_mean,vb_median=vb_median,vb_std=vb_std,vb_num=vb_num)
