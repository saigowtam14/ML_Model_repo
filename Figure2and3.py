# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:50:00 2022

@author: SAI GOWTAM VALLURI
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import pyIGRF
import pyproj
import math
import warnings
warnings.filterwarnings("ignore")

data = np.loadtxt('potelej.txt');
test1 = data;
mlon = test1[:,6];
colat = test1[:,5];
fac = test1[:,7];
pot = test1[:,10];
pot = pot/1000;


cpcp = np.round(np.max(pot)-np.min(pot))

et = test1[:,11]
ep = test1[:,12]
ctiot = test1[:,13]
ctiop = test1[:,14]
sigp = test1[:,8]
sigh = test1[:,9]

colat = np.array(np.reshape(colat,[50,25]))
mlon = np.array(np.reshape(mlon,[50,25]))
fac = np.array(np.reshape(fac,[50,25]))
pot = np.array(np.reshape(pot,[50,25]))
et = np.array(np.reshape(et,[50,25]))
ep = np.array(np.reshape(ep,[50,25]))
ctiot = np.array(np.reshape(ctiot,[50,25]))
ctiop = np.array(np.reshape(ctiop,[50,25]))
sigp = np.array(np.reshape(sigp,[50,25]))
sigh = np.array(np.reshape(sigh,[50,25]))
mlon = mlon + 180;

xjh=ctiot*et+ctiop*ep
etbar = et/(np.sqrt(et**2+ep**2))
epbar = ep/(np.sqrt(et**2+ep**2))
jpt = ctiot*etbar;
jpp = ctiop*epbar;
jht = ctiot-jpt;
jhp = ctiop-jpp;


plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.weight': 'bold'})
fig,axes = plt.subplots(figsize=(15,7),dpi=300)
plt.subplots_adjust(hspace=0.25, wspace=0.1)

x1 = colat;
y1 = np.deg2rad(mlon);
z1 = fac;
vmin = -0.6;
vmax = 0.6;
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cmap = 'seismic'
ax1 = plt.subplot(1,3,1,projection='polar')
qcs1 = ax1.pcolormesh(y1, x1, z1, norm=norm, cmap=cmap)
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
# ax1.set_title('(a) Field Aligned Currents '+ fname,fontweight="bold")
ax1.set_title('(a) Field Aligned Currents',fontweight="bold")
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.2,pad=0.1)
cbar1.set_label('($\u03BC A/m^2$)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
plt.tight_layout()

z1 = sigp;
vmin = 5;
vmax = 15;
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(1,3,2,projection='polar')
qcs1 = ax1.pcolormesh(y1, x1, z1, norm=norm, cmap='Reds')
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
ax1.set_title('(b) Pedersen conductance',fontweight="bold")
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.2,pad=0.1)
cbar1.set_label('(S)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
plt.tight_layout()


z1 = sigh;
vmin = 5;
vmax = 25;
ax1 = plt.subplot(1,3,3,projection='polar')
norm = colors.Normalize(vmin=vmin, vmax=vmax)
qcs1 = ax1.pcolormesh(y1, x1, z1, norm=norm, cmap='Reds')
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.set_title('(c) Hall conductance',fontweight="bold")
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.2,pad=0.1)
cbar1.set_label('(S)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
plt.tight_layout()



plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.weight': 'bold'})
fig,axes = plt.subplots(figsize=(15,7),dpi=300)
plt.subplots_adjust(hspace=0.25, wspace=0.1)
   
x1 = colat;
y1 = np.deg2rad(mlon);
z1 = pot;

vmin = -50;
vmax = 50;
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(1,3,1,projection='polar')
qcs1 = ax1.pcolormesh(y1, x1, z1, norm=norm, cmap=cmap)
cs = ax1.contour(y1, x1, z1, colors = 'k',levels=range(-40,41,10))
plt.clabel(cs, fmt = '%2.1d', colors = 'k', fontsize=9)
cpcp1 = 'Cross Polar Cap Potential = '+str(cpcp)+' kV'
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.set_title('(a) Potential',fontweight="bold")
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.15)
cbar1.set_label('(kV)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
plt.tight_layout()


file_name = 'dms_20130514_17s1.001.hdf5.txt'
start_tm = '20130514 12:00:00'
end_tm = '20130514 13:10:00'
start_tm = pd.to_datetime(start_tm)
end_tm = pd.to_datetime(end_tm)
tmpp = pd.to_datetime('20130514 12:44:00')

cols = ['year','month','day','hour','minute','second','RECNO','RECNO1','KINDAT','UT1_UNIX','UT2_UNIX',\
        'gdlat','glon','gdalt','SAT_ID','mlt','mlat','MLONG','NE','hor_ion_v','vert_ion_v','BD',\
            'B_FORWARD','B_PERP', 'DIFF_BD','DIFF_B_FOR','DIFF_B_PERP']
dmsp1 = pd.read_csv(file_name,names=cols,skiprows=1,delim_whitespace=True,low_memory=False)
dmsp1['timestamps'] = pd.to_datetime(dmsp1[['month','day','year','hour','minute']]);
dmsp1 = dmsp1[dmsp1['mlat'].between(40,90)]
dmsp1 = dmsp1.reset_index()
dmsp1 = dmsp1[dmsp1['hor_ion_v'].between(-3000,3000)]
dmsp1 = dmsp1.reset_index()

intind = (dmsp1['timestamps']>=start_tm) & (dmsp1['timestamps']<=end_tm)
dmsp1 = dmsp1.loc[intind]
dmsp1.reset_index(inplace=True,drop=True)


lat = np.array(dmsp1['gdlat'])
lon = np.array(dmsp1['glon'])
alt = np.array(dmsp1['gdalt'])
date = np.repeat(2013,np.size(alt))
igrf = pd.DataFrame([])
for i1 in range(0,np.size(alt)):
    ig =  pyIGRF.igrf_value(lat[i1], lon[i1], alt[i1], date[i1])
    ig = np.transpose(ig)
    tmpig = pd.DataFrame(ig)
    igrf = pd.concat([igrf,tmpig])

igrf = np.reshape(np.array(igrf),[int(np.size(igrf)/7),7])

Vy = np.array(dmsp1['hor_ion_v'])
Vz = np.array(dmsp1['vert_ion_v'])
Bx = igrf[:,3]*(10**-9)
By = igrf[:,4]*(10**-9)
Bz = igrf[:,5]*(10**-9)
dmsp1['Bx'] = Bx
dmsp1['By'] = By
dmsp1['Bz'] = Bz
 
Ex = -Vy*Bz + Vz*By
dmsp1['Ex'] = Ex 
Exd1 = Ex


x1 = 90-dmsp1['mlat'];
y1 = dmsp1['mlt']*15;
y1 = np.deg2rad(y1);
z1 = dmsp1['hor_ion_v'];
qcs1 = ax1.quiver(y1, x1, y1, x1+z1,scale=15000,color='m')

tmm = pd.to_datetime(tmpp)
title2 = str(tmm.year)+ "{:02d}".format(tmm.month) +\
    "{:02d}".format(tmm.day) +' '+"{:02d}".format(tmm.hour) + \
        ':'+ "{:02d}".format(tmm.minute) + 'UT'
# plt.title('(b) Weimer model')
cpcp1 = '$\u03A6_{PC}$ = '+str(cpcp)+' kV'
ax1.text(1, 1, cpcp1, transform=ax1.transAxes, fontsize = 12)




ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')


for i2 in range(0,np.size(alt)-1):
    a = pyproj.transform(lla, ecef, lon[i2], lat[i2], alt[i2], radians=False)
    b = pyproj.transform(lla, ecef, lon[i2+1], lat[i2+1], alt[i2+1], radians=False)
    d = math.dist(b, a)
    Exd1[i2] = Ex[i2]*(d)

Exd1[Exd1>750] = np.nan
Exd1[Exd1<-750] = np.nan

dmsp1['Ex*dl_igrf'] = Exd1

tmp1 = np.array(dmsp1['Ex*dl_igrf'].cumsum())

dl_cor = []
for i2 in range(0,np.size(alt)):
    a1 = pyproj.transform(lla, ecef, lon[0], lat[0], alt[0], radians=False)
    b1 = pyproj.transform(lla, ecef, lon[-1], lat[-1], alt[-1], radians=False)
    d1 = math.dist(b1, a1)
    
    a = pyproj.transform(lla, ecef, lon[0], lat[0], alt[0], radians=False)
    b = pyproj.transform(lla, ecef, lon[i2], lat[i2], alt[i2], radians=False)
    d = math.dist(b, a)
    
    dddd = d/d1
    dl_cor = np.append(dl_cor,dddd)
    # print(d/d1)
    # dl_cor.iloc[i1] = d/d1
    


tmp1 = tmp1 - tmp1[-1]*dl_cor
cpdmsp = np.round(max(tmp1/1000)-min(tmp1/1000))
cpcp1 = 'DMSP $\u03A6_{PC}$ = '+str(cpdmsp)+' kV'
ax1.text(1, 1.1, cpcp1, transform=ax1.transAxes, fontsize = 12, color='m')


# plt.rcParams.update({'font.size': 14})
# plt.rcParams.update({'font.weight': 'bold'})
# fig,axes = plt.subplots(figsize=(8,6),dpi=300)
# plt.subplots_adjust(hspace=0.25, wspace=0.1)

x1 = colat;
y1 = np.deg2rad(mlon);
z1 = xjh;
vmin = 0.0;
vmax = 0.1;
norm = colors.Normalize(vmin=vmin, vmax=vmax)
ax1 = plt.subplot(1,3,2,projection='polar')
qcs1 = ax1.pcolormesh(y1, x1, z1, norm=norm, cmap='Reds')
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.set_title('(b) Joule Heating rate',fontweight="bold")
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.15)
cbar1.set_label('$(W/m^2)$',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
plt.tight_layout()

# plt.rcParams.update({'font.size': 14})
# plt.rcParams.update({'font.weight': 'bold'})
# fig,axes = plt.subplots(figsize=(8,6),dpi=300)
# plt.subplots_adjust(hspace=0.25, wspace=0.1)

x1 = colat;
y1 = np.deg2rad(mlon);
z1 = jpp;
z2 = jpt
jpnorm = np.sqrt(z1**2 + z2**2)
# vmin = -2.5;
# vmax = 2.5;
# norm = colors.Normalize(vmin=vmin, vmax=vmax)
# cmap = 'seismic'
# ax1 = plt.subplot(1,2,1,projection='polar')
# qcs1 = ax1.pcolormesh(y1, x1, jpnorm, cmap='Reds')
# ax1.quiver(y1, x1, z1/jpnorm, z2/jpnorm, angles = 'xy', color='magenta')
# ax1.set_theta_zero_location("S")
# ax1.set_ylim(0, 30)
# ax1.set_title('Pedersen currents',fontweight="bold")
# ax1.xaxis.grid(linestyle='--', linewidth=0.8)
# ax1.yaxis.grid(linestyle='--', linewidth=0.8)
# xtickpos = ax1.get_xticks()
# cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.15)
# xticks = ['00', '03','06','09','12','15','18','21']
# ax1.set_xticks(xtickpos,xticks)
# plt.tight_layout()


x1 = colat;
y1 = np.deg2rad(mlon);
z1 = jhp;
z2 = jht
jhnorm = np.sqrt(z1**2 + z2**2)
vmin = 0;
vmax = 4;
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cmap = 'seismic'
ax1 = plt.subplot(1,3,3,projection='polar')
qcs1 = ax1.pcolormesh(y1, x1, jhnorm, norm = norm, cmap='Reds')
ax1.quiver(y1, x1, z1/jhnorm, z2/jhnorm, angles = 'xy', color='magenta')
ax1.set_theta_zero_location("S")
ax1.set_ylim(0, 40)
ax1.set_title('(c) Hall currents',fontweight="bold")
ax1.xaxis.grid(linestyle='--', linewidth=0.8)
ax1.yaxis.grid(linestyle='--', linewidth=0.8)
cbar1 = fig.colorbar(qcs1, ax=ax1,shrink=0.4,pad=0.15)
cbar1.set_label('(A/m)',fontweight="bold")
xtickpos = ax1.get_xticks()
xticks = ['00 MLT', '03','06','09','12','15','18','21']
ax1.set_xticks(xtickpos,xticks)
ytickpos = ax1.get_yticks()
yticks = ['N', '80$^\circ$','','70$^\circ$','','60$^\circ$','','50$^\circ$']
ax1.set_yticks(ytickpos,yticks)
plt.tight_layout()