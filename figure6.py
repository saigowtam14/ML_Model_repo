# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 09:50:58 2022

@author: SAI GOWTAM VALLURI
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import warnings
from datetime import timedelta
warnings.filterwarnings("ignore")
# =============================================================================
with open('20170907_2201.pkl','rb') as f:  # Python 3: open(..., 'rb')
    mlx, mly, mlz, dmx, dmy, dmz, wx, wy, wz, sdx, sdy, sdz, cpdmsp,\
        weimertm,tmpp,vmlats,vmlons = pickle.load(f)
f.close()

plt.rcParams['axes.grid'] = False
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.weight': 'bold'})
fig,axes = plt.subplots(figsize=(20,15),dpi=300)
plt.subplots_adjust(hspace=0, wspace=0)
   
cpcp =  np.round(np.max(mlz)-np.min(mlz))
cpcp = np.round(cpcp)
vmin = -75
vmax = 75
cmap = 'seismic'
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(3,3,1,projection='polar')
qcs1 = ax1.pcolormesh(mly, mlx, mlz, norm=norm, cmap = cmap)
cs = ax1.contour(mly, mlx, mlz, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(a) ML-Model',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
tmpp = pd.to_datetime(tmpp) + timedelta(minutes=1)
ax1.text(0.9, 1.1, tmpp, transform=ax1.transAxes, fontsize = 12, color='m')


wxc = wx[wx<=28]
wyc = wy
wzc = wz[0:29,:]

cpcp =  np.round(np.max(wz)-np.min(wz))
ax1 = plt.subplot(3,3,2,projection='polar')
qcs1 = ax1.pcolormesh(wy, wx, wz,norm=norm, cmap = cmap)
cs = ax1.contour(wyc, wxc, wzc, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
cpcp1 = 'Cross Polar Cap Potential = '+str(cpcp)+' kV'
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
# ax1.set_title(cpcp1,fontweight="bold")
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(b) Weimer 2005',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
cpcp1 = 'DMSP $\u03A6_{PC}$ = '+str(cpdmsp)+' kV'
ax1.text(0.9, 1.1, cpcp1, transform=ax1.transAxes, fontsize = 12, color='m')



cpcp = np.round(np.max(sdz)-np.min(sdz))
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(3,3,3,projection='polar')
qcs1 = ax1.pcolormesh(sdy, sdx, sdz, norm=norm, cmap=cmap)
cs = ax1.contour(sdy, sdx, sdz, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=7)
ax1.set_theta_zero_location("S")
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.scatter(vmlons,90-vmlats,s=3,color='g')
ax1.grid(linestyle='--', linewidth=0.8)
ax1.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(c) SuperDARN',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)



with open('20170907_2342.pkl','rb') as f:  # Python 3: open(..., 'rb')
    mlx, mly, mlz, dmx, dmy, dmz, wx, wy, wz, sdx, sdy, sdz, cpdmsp,\
        weimertm,tmpp,vmlats,vmlons = pickle.load(f)
f.close()

  
cpcp =  np.round(np.max(mlz)-np.min(mlz))
cpcp = np.round(cpcp)
vmin = -75
vmax = 75
cmap = 'seismic'
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(3,3,4,projection='polar')
qcs1 = ax1.pcolormesh(mly, mlx, mlz, norm=norm, cmap = cmap)
cs = ax1.contour(mly, mlx, mlz, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.grid(linestyle='--', linewidth=0.8)
ax1.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(d) ML-Model',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
ax1.text(0.9, 1.1, tmpp, transform=ax1.transAxes, fontsize = 12, color='m')


wxc = wx[wx<=28]
wyc = wy
wzc = wz[0:29,:]

cpcp =  np.round(np.max(wz)-np.min(wz))
ax1 = plt.subplot(3,3,5,projection='polar')
qcs1 = ax1.pcolormesh(wy, wx, wz,norm=norm, cmap = cmap)
cs = ax1.contour(wyc, wxc, wzc, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
cpcp1 = 'Cross Polar Cap Potential = '+str(cpcp)+' kV'
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
# ax1.set_title(cpcp1,fontweight="bold")
ax1.grid(linestyle='--', linewidth=0.8)
ax1.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(e) Weimer 2005',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
cpcp1 = 'DMSP $\u03A6_{PC}$ = '+str(cpdmsp)+' kV'
ax1.text(0.9, 1.1, cpcp1, transform=ax1.transAxes, fontsize = 12, color='m')



cpcp = np.round(np.max(sdz)-np.min(sdz))
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(3,3,6,projection='polar')
qcs1 = ax1.pcolormesh(sdy, sdx, sdz, norm=norm, cmap=cmap)
cs = ax1.contour(sdy, sdx, sdz, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=7)
ax1.set_theta_zero_location("S")
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.scatter(vmlons,90-vmlats,s=3,color='g')
ax1.grid(linestyle='--', linewidth=0.8)
ax1.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(f) SuperDARN',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
plt.tight_layout()




with open('20170908_1503.pkl','rb') as f:  # Python 3: open(..., 'rb')
    mlx, mly, mlz, dmx, dmy, dmz, wx, wy, wz, sdx, sdy, sdz, cpdmsp,\
        weimertm,tmpp,vmlats,vmlons = pickle.load(f)
f.close()

# plt.rcParams['axes.grid'] = False
# plt.rcParams.update({'font.size': 12})
# plt.rcParams.update({'font.weight': 'bold'})
# fig,axes = plt.subplots(figsize=(15,8),dpi=600)
# plt.subplots_adjust(hspace=0, wspace=0)
   
cpcp =  np.round(np.max(mlz)-np.min(mlz))
cpcp = np.round(cpcp)
vmin = -75
vmax = 75
cmap = 'seismic'
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(3,3,7,projection='polar')
qcs1 = ax1.pcolormesh(mly, mlx, mlz, norm=norm, cmap = cmap)
cs = ax1.contour(mly, mlx, mlz, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(g) ML-Model',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
tmpp = pd.to_datetime(tmpp) + timedelta(minutes=1)
ax1.text(0.9, 1.1, tmpp, transform=ax1.transAxes, fontsize = 12, color='m')


wxc = wx[wx<=24]
wyc = wy
wzc = wz[0:25,:]

cpcp =  np.round(np.max(wz)-np.min(wz))
ax1 = plt.subplot(3,3,8,projection='polar')
qcs1 = ax1.pcolormesh(wy, wx, wz,norm=norm, cmap = cmap)
cs = ax1.contour(wyc, wxc, wzc, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
cpcp1 = 'Cross Polar Cap Potential = '+str(cpcp)+' kV'
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
# ax1.set_title(cpcp1,fontweight="bold")
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(h) Weimer 2005',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
cpcp1 = 'DMSP $\u03A6_{PC}$ = '+str(cpdmsp)+' kV'
ax1.text(0.9, 1.1, cpcp1, transform=ax1.transAxes, fontsize = 12, color='m')



cpcp = np.round(np.max(sdz)-np.min(sdz))
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(3,3,9,projection='polar')
qcs1 = ax1.pcolormesh(sdy, sdx, sdz, norm=norm, cmap=cmap)
cs = ax1.contour(sdy, sdx, sdz, colors = 'k',levels=range(-60,61,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=7)
ax1.set_theta_zero_location("S")
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.scatter(vmlons,90-vmlats,s=3,color='g')
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
cbar1.set_label('Potential (kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
plt.title('(i) SuperDARN',fontweight="bold")
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)



# with open('20170908_1645.pkl','rb') as f:  # Python 3: open(..., 'rb')
#     mlx, mly, mlz, dmx, dmy, dmz, wx, wy, wz, sdx, sdy, sdz, cpdmsp,\
#         weimertm,tmpp,vmlats,vmlons = pickle.load(f)
# f.close()

  
# cpcp =  np.round(np.max(mlz)-np.min(mlz))
# cpcp = np.round(cpcp)
# vmin = -75
# vmax = 75
# cmap = 'seismic'
# norm = colors.Normalize(vmin=vmin, vmax=vmax)
# ax1 = plt.subplot(2,3,4,projection='polar')
# qcs1 = ax1.pcolormesh(mly, mlx, mlz, norm=norm, cmap = cmap)
# cs = ax1.contour(mly, mlx, mlz, colors = 'k',levels=range(-60,61,10))
# plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
# ax1.set_theta_zero_location("S")
# ax1.set_ylim(0, 40)
# ax1.xaxis.grid(linestyle='--', linewidth=0.8)
# ax1.yaxis.grid(linestyle='--', linewidth=0.8)
# cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
# cbar1.set_label('Potential (kV)',fontweight="bold")
# xtickpos = ax1.get_xticks()
# xticks = ['00 MLT', '03','06','09','12','15','18','21']
# ax1.set_xticks(xtickpos,xticks)
# ytickpos = ax1.get_yticks()
# yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
# ax1.set_yticks(ytickpos,yticks)
# qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
# plt.title('(d) ML-Model',fontweight="bold")
# cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
# ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
# tmpp = pd.to_datetime(tmpp) + timedelta(minutes=1)
# ax1.text(0.9, 1.1, tmpp, transform=ax1.transAxes, fontsize = 12, color='m')


# wxc = wx[wx<=27]
# wyc = wy
# wzc = wz[0:28,:]

# cpcp =  np.round(np.max(wz)-np.min(wz))
# ax1 = plt.subplot(2,3,5,projection='polar')
# qcs1 = ax1.pcolormesh(wy, wx, wz,norm=norm, cmap = cmap)
# cs = ax1.contour(wyc, wxc, wzc, colors = 'k',levels=range(-60,61,10))
# plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=6)
# cpcp1 = 'Cross Polar Cap Potential = '+str(cpcp)+' kV'
# ax1.set_theta_zero_location("S")
# ax1.set_ylim(0, 40)
# # ax1.set_title(cpcp1,fontweight="bold")
# ax1.xaxis.grid(linestyle='--', linewidth=0.8)
# ax1.yaxis.grid(linestyle='--', linewidth=0.8)
# cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
# cbar1.set_label('Potential (kV)',fontweight="bold")
# xtickpos = ax1.get_xticks()
# xticks = ['00 MLT', '03','06','09','12','15','18','21']
# ax1.set_xticks(xtickpos,xticks)
# ytickpos = ax1.get_yticks()
# yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
# ax1.set_yticks(ytickpos,yticks)
# qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
# plt.title('(e) Weimer 2005',fontweight="bold")
# cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
# ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
# cpcp1 = 'DMSP $\u03A6_{PC}$ = '+str(cpdmsp)+' kV'
# ax1.text(0.9, 1.1, cpcp1, transform=ax1.transAxes, fontsize = 12, color='m')



# cpcp = np.round(np.max(sdz)-np.min(sdz))
# norm = colors.Normalize(vmin=vmin, vmax=vmax)
# ax1 = plt.subplot(2,3,6,projection='polar')
# qcs1 = ax1.pcolormesh(sdy, sdx, sdz, norm=norm, cmap=cmap)
# cs = ax1.contour(sdy, sdx, sdz, colors = 'k',levels=range(-60,61,10))
# plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=7)
# ax1.set_theta_zero_location("S")
# ax1.set_theta_zero_location("S")
# ax1.set_ylim(0, 40)
# ax1.scatter(vmlons,90-vmlats,s=3,color='g')
# ax1.xaxis.grid(linestyle='--', linewidth=0.8)
# ax1.yaxis.grid(linestyle='--', linewidth=0.8)
# cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.1)
# cbar1.set_label('Potential (kV)',fontweight="bold")
# xtickpos = ax1.get_xticks()
# xticks = ['00 MLT', '03','06','09','12','15','18','21']
# ax1.set_xticks(xtickpos,xticks)
# ytickpos = ax1.get_yticks()
# yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
# ax1.set_yticks(ytickpos,yticks)
# qcs1 = ax1.quiver(dmy, dmx, dmy, dmx+dmz,scale=15000,color='m')
# plt.title('(f) SuperDARN',fontweight="bold")
# cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
# ax1.text(0.95, 0.8, cpcp1, transform=ax1.transAxes, fontsize = 12)
# plt.tight_layout()