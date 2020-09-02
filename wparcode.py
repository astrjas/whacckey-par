#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import wparfuncs
from astropy.io import fits
#import scipy.interpolate as sin
import matplotlib as mpl

from astropy.modeling import blackbody
#from astropy import units as u
#from astropy.visualization import quantity_support

def convert(s):
    new=""
    for x in s:
        new+=x
    return new

# This code block is for adjusting plot visualization. It is not necessary for the core functionality.
# NOTE: Individual visualization parameters can be adjusted on the fly (see below).
rcParams['figure.figsize'] = (12, 8)
rcParams['axes.titlesize'] = 16
rcParams['axes.labelsize'] = 16
rcParams['legend.fontsize'] = 16
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14

medmod=19.01
medmoderr=0.08

power=(medmod/5)+1
dsmc=10**power
dsmcerr=0.2*np.log(10)*(10**((0.2*medmod)+1))*medmoderr

dsmc_kpc=dsmc/1000.
dsmcerr_kpc=dsmcerr/1000

dsmc_cm=dsmc*(3.086e18) #cm
dsmcerr_cm=dsmcerr*(3.086e18)
#%%
#LIN 358 DATA
hdul358 = fits.open("lin358.fits")
headers358=hdul358[0].data
data358=hdul358[1].data

fluxjyall358=data358.field('_sed_flux') #in Jy
fluxejyall358=data358.field('_sed_eflux') #in Jy
freqGHzall358=data358.field('_sed_freq') #in GHz

fluxjy358=fluxjyall358[~np.isnan(fluxejyall358)]
fluxejy358=fluxejyall358[~np.isnan(fluxejyall358)]
freqGHz358=freqGHzall358[~np.isnan(fluxejyall358)]

freq358=freqGHz358*1.0e9
wave358=(3.0e18)/freq358 #in Angstroms
flux358=fluxjy358*(1e-23)*freq358/wave358
fluxe358=fluxejy358*(1e-23)*freq358/wave358

l=0
filter358=([])
filter358allarr=data358.field('_sed_filter')
filter358arr=filter358allarr[~np.isnan(fluxejyall358)]
for l in np.arange(len(filter358arr)): 
    filter358=np.append(filter358,convert(filter358arr[l]))

dat358={'flux':flux358,'ferror':fluxe358,'filter':filter358}
table358=pd.DataFrame(data=dat358,index=wave358,columns=['flux','ferror','filter'])

wu=list(np.unique(wave358))

spec358,wavetable358=wparfuncs.spec_vot(table358,wave358)
#noot=wavetable358.at['2MASS:Ks','wavelength']
#wu.index(noot)

hdulx1=fits.open("xrayflux1001.fits")
headersx1=hdulx1[0].data
datax1=hdulx1[1].data

channelallx1=datax1.field('CHANNEL')
fluxpallx1=datax1.field('FLUX') #in photons/
fluxpeallx1=datax1.field('ERROR')

channelx1=channelallx1[~np.isnan(fluxpallx1)]
fluxpx1=fluxpallx1[~np.isnan(fluxpallx1)]
fluxpex1=fluxpeallx1[~np.isnan(fluxpallx1)]

fluxx1=np.zeros(len(fluxpx1))
fluxex1=np.zeros(len(fluxpex1))

channelx1cm=channelx1*(1.0e-8)

p=0
for p in range(len(channelx1)):
    fluxx1[p]=fluxpx1[p]*(((6.626e-27)*(3.0e10))/channelx1cm[p])
    fluxex1[p]=fluxpex1[p]*(((6.626e-27)*(3.0e10))/channelx1cm[p])    


hdulx3=fits.open("xrayflux1003.fits")
headersx3=hdulx3[0].data
datax3=hdulx3[1].data

channelallx3=datax3.field('CHANNEL')
fluxpallx3=datax3.field('FLUX') #in photons/
fluxpeallx3=datax3.field('ERROR')

channelx3=channelallx3[~np.isnan(fluxpallx3)]
fluxpx3=fluxpallx3[~np.isnan(fluxpallx3)]
fluxpex3=fluxpeallx3[~np.isnan(fluxpallx3)]

fluxx3=np.zeros(len(fluxpx3))
fluxex3=np.zeros(len(fluxpex3))

channelx3cm=channelx3*(1.0e-8)

p=0
for p in range(len(channelx3)):
    fluxx3[p]=fluxpx3[p]*(((6.626e-27)*(3.0e10))/channelx3cm[p])
    fluxex3[p]=fluxpex3[p]*(((6.626e-27)*(3.0e10))/channelx3cm[p])    


#channel=np.zeros(len())
#((6.626e-27)*(3.0e10))/channel

#fig, axs = plt.subplots(figsize=(12,8))
#plt.rc('font',family='Times New Roman')
#axs.scatter(np.log10(channelx3),np.log10(fluxx3),color='blue',label='Scaled Model Kurucz',zorder=1)
#axs.set_xlabel('log(wavelength) [A]')  # This is the a the ordinate, or vertical axis.
#axs.set_ylim(bottom=-18.5,top=-13.6)#bscissa, or horizontal axis.
#axs.set_ylabel('log(flux) [erg/cm2/s/A]')  # This is
#axs.set_xlim(left=3.0,right=5.5)


#%%
#N73 DATA
hdul73 = fits.open("n73.fits")
headers73=hdul73[0].data
data73=hdul73[1].data

fluxjyall73=data73.field('_sed_flux') #in Jy
fluxejyall73=data73.field('_sed_eflux')
freqGHzall73=data73.field('_sed_freq') #in GHz

fluxjy73=fluxjyall73[~np.isnan(fluxejyall73)]
fluxejy73=fluxejyall73[~np.isnan(fluxejyall73)]
freqGHz73=freqGHzall73[~np.isnan(fluxejyall73)]

freq73=freqGHz73*1.0e9
wave73=(3.0e18)/freq73 #in Angstroms
flux73=fluxjy73*(1e-23)*freq73/wave73
fluxe73=fluxejy73*(1e-23)*freq73/wave73

n=0
filter73=([])
filter73allarr=data73.field('_sed_filter')
filter73arr=filter73allarr[~np.isnan(fluxejyall73)]
for n in np.arange(len(filter73arr)):
    filter73=np.append(filter73,convert(filter73arr[n]))

dat73={'flux':flux73,'ferror':fluxe73,'filter':filter73}
table73=pd.DataFrame(data=dat73,index=wave73,columns=['flux','ferror','filter'])

spec73,wavetable73=wparfuncs.spec_vot(table73,wave73)
#%%
#KURUCZ
met=np.array(["p00","m10"])
logg=np.array(["1.0","0.5"])
temp=3750 #goes to 6500 in incr of 250

#chii358,model_kur_name358,scale_kur358,flux_kur358=twoobs_func.fits_kur(spec358)
chii358,model_kur_name358,scale_kur358,flux_kur358,band358,funck358=wparfuncs.fits_kur_apogee(spec358,3750,"1.0","p00",wavetable358)
kshort358=model_kur_name358[7:47]
ktemp358=kshort358[22:26]
kg358=kshort358[33:36]
kmet358=kshort358[1:4]
wkur358,fkur358=wparfuncs.kurucz_getter(kshort358)

alpha358=scale_kur358 #=(r/d)^2
rms358=dsmc_cm*np.sqrt(alpha358) #in cm
rmserr358=dsmcerr_cm*np.sqrt(alpha358)
rms358_sol=rms358/(6.96e10)
rmserr358_sol=rmserr358/(6.96e10)

mms358=((rms358**2)/(6.67e-8))*10.**(float(kg358))
mmserr358=((2.*rms358)/(6.67e-8))*rmserr358*10.**(float(kg358))
mms358_sol=mms358/(1.99e33)
mmserr358_sol=mmserr358/(1.99e33)

#%%
#chii73,model_kur_name73,scale_kur73,flux_kur73=twoobs_func.fits_kur(spec73)
chii73,model_kur_name73,scale_kur73,flux_kur73,band73,funck73=wparfuncs.fits_kur_apogee(spec73,3750,"0.5","m10",wavetable73)
kshort73=model_kur_name73[7:47]
kg73=kshort73[33:36]
kmet73=kshort73[1:4]
ktemp73=kshort73[22:26]
wkur73,fkur73=wparfuncs.kurucz_getter(kshort73)

alpha73=scale_kur73 #=(r/d)^2
rms73=dsmc_cm*np.sqrt(alpha73) #in cm
rmserr73=dsmcerr_cm*np.sqrt(alpha73) #in cm
rms73_sol=rms73/(6.96e10)
rmserr73_sol=rmserr73/(6.96e10)

mms73=((rms73**2)/(6.67e-8))*10.**(float(kg73))
mmserr73=((2.*rms73)/(6.67e-8))*rmserr73*10.**(float(kg73))
mms73_sol=mms73/(1.99e33)
mmserr73_sol=mmserr73/(1.99e33)

#%%
#BB
chii_bb358,model_bb_name358,scale_bb358,flux_bb358,bandbb358,funcb358,indb358=wparfuncs.fits_bb(spec358,wavetable358)
bbtemp358=np.float(model_bb_name358[8:14])
#fbb358=twoobs_func.bb(bbtemp358,wkur358)
#fbb358=flux_bb358[1]


#fbb358=sin.interp1d(x=flux_bb358[0],y=flux_bb358[1],kind='linear')
wkuradjust358=wkur358[(wkur358>=np.min(flux_bb358[0])) & (wkur358<=np.max(flux_bb358[0]))]
fkuradjust358=fkur358[(wkur358>=np.min(flux_bb358[0])) & (wkur358<=np.max(flux_bb358[0]))]

fbb358k = blackbody.blackbody_lambda(wkuradjust358,temperature=bbtemp358)
fbb358 = fbb358k.value

beta358=scale_bb358 #=(r/d)^2
rbb358=dsmc_cm*np.sqrt(beta358) #in cm
rbb358_sol=rbb358/(6.96e10)
rbb358_e=rbb358/(6.37e8)

rbberr358=dsmcerr_cm*np.sqrt(beta358) #in cm
rbberr358_sol=rbberr358/(6.96e10)
rbberr358_e=rbberr358/(6.37e8)


#((0.3*(6.96e10))/dsmc_cm)**2

#%%
print("N73")
chii_bb73,model_bb_name73,scale_bb73,flux_bb73,bandbb73,funcb73,indb73=wparfuncs.fits_bb(spec73,wavetable73)
bbtemp73=model_bb_name73[8:14]
#fbb73=twoobs_func.bb(bbtemp73,wkur73)
#fbb73=flux_bb73[1]
#beta73=np.polyfit()

#fbb73=sin.interp1d(x=flux_bb73[0],y=flux_bb73[1],kind='linear')
wkuradjust73=wkur73[(wkur73>=np.min(flux_bb73[0])) & (wkur73<=np.max(flux_bb73[0]))]
fkuradjust73=fkur73[(wkur73>=np.min(flux_bb73[0])) & (wkur73<=np.max(flux_bb73[0]))]
#ifbb73=funcb73(wkuradjust73)

fbb73k = blackbody.blackbody_lambda(wkuradjust73,temperature=bbtemp73)
fbb73 = fbb73k.value

beta73=scale_bb73 #=(r/d)^2
rbb73=dsmc_cm*np.sqrt(beta73) #in cm
rbb73_sol=rbb73/(6.96e10)
rbb73_e=rbb73/(6.37e8)

rbberr73=dsmcerr_cm*np.sqrt(beta73) #in cm
rbberr73_sol=rbberr73/(6.96e10)
rbberr73_e=rbberr73/(6.37e8)

#beta=np.polyfit(wavelength_bb,beta73*flux_bb,2)

#%%
#BOTH
scaleflux_both358=(scale_kur358*fkuradjust358)+(scale_bb358*fbb358)
scaleflux_both73=(scale_kur73*fkuradjust73)+(scale_bb73*fbb73)

#scaleflux_both358=(scale_kur358*flux_kur358[1])+(scale_bb358*flux_bb358[1])
#scaleflux_both73=(scale_kur73*flux_kur73[1])+(scale_bb73*flux_bb73[1])
#%%

"""
wcommon358=np.array([])
scaleflux_both358=np.array([])

l=0
for l in np.arange(len(infokur358[0])):
    k=0
    if infokur358[0][l]<=0: continue
    for k in np.arange(len(infobb358[0])):
        if infobb358[0][k]<=0: continue
        if infokur358[0][l]==infobb358[0][k]:
            wcommon358=np.append(wcommon358,infokur358[0][l])
            scaleflux_both358=np.append(scaleflux_both358,(scale_kur358*infokur358[1][l])+(scale_bb358*infobb358[1][k]))
"""

"""
wcommon73=np.array([])
scaleflux_both73=np.array([])

d=0
for d in np.arange(len(infokur73[0])):
    j=0
    if infokur73[0][d]<=0: continue
    for j in np.arange(len(infobb73[0])):
        if infobb73[0][j]<=0: continue
        if infokur73[0][d]==infobb73[0][j]:
            wcommon73=np.append(wcommon73,infokur73[0][d])
            scaleflux_both73=np.append(scaleflux_both73,(scale_kur73*infokur73[1][d])+(scale_bb73*infobb73[1][j]))
"""

#chii358,model_kur_name358,model_bb_name358,scale358,flux_kur358,flux_bb358,flux_both358=twoobs_func.fits_both(spec358)
#chii73,model_kur_name73,model_bb_name73,scale73,flux_kur73,flux_bb73,flux_both73=twoobs_func.fits_both(spec73)

"""
temp358=twoobs_func.fits_temp(spec358)
if temp358<110000:
    chii358,model_kur_name358,model_bb_name358,scale358,flux_kur358,flux_bb358,flux_both358=twoobs_func.fits_both_specific(spec358,start=100000,stop=110000,step=1000)
    scale_kur358,scale_bb358=scale358
else:
    chii358,model_kur_name358,model_bb_name358,scale358,flux_kur358,flux_bb358,flux_both358=twoobs_func.fits_both_specific(spec358,temp358-5000,temp358+5000,1000)
    scale_kur358,scale_bb358=scale358

temp73=twoobs_func.fits_temp(spec73)
if temp73<110000:
    chii73,model_kur_name73,model_bb_name73,scale73,flux_kur73,flux_bb73,flux_both73=twoobs_func.fits_both_specific(spec73,100000,110000,1000)
    scale_kur73,scale_bb73=scale73
else:    
    chii73,model_kur_name73,model_bb_name73,scale73,flux_kur73,flux_bb73,flux_both73=twoobs_func.fits_both_specific(spec73,temp73-5000,temp73+5000,1000)
    scale_kur73,scale_bb73=scale73
"""
#%%
fig, axs = plt.subplots(figsize=(12,8))
plt.rc('font',family='Times New Roman')
axs.plot(np.log10(wkuradjust358),np.log10(scale_kur358*fkuradjust358),color='blue',label='Scaled Model Kurucz',zorder=1)
axs.plot(np.log10(wkuradjust358),np.log10(scale_bb358*fbb358),color='red',label='Scaled Model Blackbody',zorder=2)
#axs.scatter(np.log10(wkuradjust358),np.log10(scale_bb358*ifbb358),color='red',marker='*',label='Scaled Model Blackbody',zorder=2)
axs.plot(np.log10(wkuradjust358),np.log10(scaleflux_both358),color='purple',label='Scaled Model Sum',zorder=3)
axs.scatter(np.log10(spec358[0][:]), np.log10(spec358[1][:]),color='orange',label='Observed Flux',zorder=4)
e=0
err358=np.zeros((2,len(spec358[2])))
err358[0]=spec358[1]-spec358[2]
err358[1]=spec358[1]+spec358[2]
for e in np.arange(len(spec358[2])):
    axs.plot([np.log10(spec358[0][e]),np.log10(spec358[0][e])], [np.log10(err358[0][e]),np.log10(err358[1][e])],color='orange',marker='_',zorder=5,markersize=10)
axs.set_title("LIN 358")
axs.set_xlabel('log(wavelength) [A]')  # This is the a the ordinate, or vertical axis.
axs.set_ylim(bottom=-18.5,top=-13.6)#bscissa, or horizontal axis.
axs.set_ylabel('log(flux) [erg/cm2/s/A]')  # This is
axs.set_xlim(left=3.0,right=5.5)
axs.text(x=3.05,y=-13.9,s="GIANT: T=3800 K, R="+str(round(rms358_sol,1))+r'$\pm$'+str(round(rmserr358_sol,1))+r'$R_\odot$'+", M="+str(round(mms358_sol,1))+r'$\pm$'+str(round(mmserr358_sol,1))+r'$M_\odot$',fontsize=16)
axs.text(x=3.05,y=-14.1,s="BB: T="+str(bbtemp358)[0:6]+r'$\pm$'+"5000 K, R="+str(round(rbb358_sol,2))+r'$\pm$'+str(round(rbberr358_sol,2))+r"$R_\odot$"+"="+str(round(rbb358_e,1))+r'$\pm$'+str(round(rbberr358_e,1))+r"$R_\oplus$",fontsize=16)
i=0
#for i in range(5):
#    axs.text(x=np.log10(spec358[0][i]), y=np.log10(spec358[1][i]),s=str(np.round(spec358[0][i],0)),fontsize=12,zorder=6)
axs.legend()
#fig.savefig("lin358models05262.png",bbox_inches='tight')
plt.show()
#%%
fig, axs = plt.subplots(figsize=(12,8))
plt.rc('font',family='Times New Roman')
axs.plot(np.log10(wkuradjust73),np.log10(scale_kur73*fkuradjust73),color='blue',label='Scaled Model Kurucz',zorder=1)
#axs.plot(np.log10(wkuradjust73),np.log10(scale_bb73*ifbb73),color='green')
axs.plot(np.log10(wkuradjust73),np.log10(scale_bb73*fbb73),color='red',label='Scaled Model Blackbody',zorder=2)
axs.plot(np.log10(wkuradjust73),np.log10(scaleflux_both73),color='purple',label='Scaled Model Sum',zorder=3)
axs.scatter(np.log10(spec73[0][:]), np.log10(spec73[1][:]), color='orange',label='Observed Flux',zorder=4)
s=0
err=np.zeros((2,len(spec73[2])))
err[0]=spec73[1]-spec73[2]
err[1]=spec73[1]+spec73[2]
#axs.errorbar(np.log10(spec73[0]), np.log10(spec73[1]),yerr=err,ecolor='orange',capsize=5)
for s in np.arange(len(spec73[2])):
    axs.plot([np.log10(spec73[0][s]),np.log10(spec73[0][s])], [np.log10(err[0][s]),np.log10(err[1][s])],color='orange',marker='_',zorder=5,markersize=10)
axs.set_title("SMC N73")
axs.set_xlabel('log(wavelength) [A]')  # This is the abscissa, or horizontal axis.
axs.set_ylabel('log(flux) [erg/cm2/s/A]')  # This is the ordinate, or vertical axis.
axs.set_ylim(bottom=-18.5,top=-13.6)
axs.set_xlim(left=3.0,right=5.5)
axs.text(x=3.05,y=-13.9,s="GIANT: T=3800 K, R="+str(round(rms73_sol,1))+r'$\pm$'+str(round(rmserr73_sol,1))+r"$R_\odot$"+", M="+str(round(mms73_sol,1))+r'$\pm$'+str(round(mmserr73_sol,1))+"$M_\odot$",fontsize=16)
axs.text(x=3.05,y=-14.1,s="BB: T="+str(bbtemp73)+r'$\pm$'+"1000 K, R="+str(round(rbb73_sol,2))+r'$\pm$'+str(round(rbberr73_sol,2))+r"$R_\odot$"+"="+str(round(rbb73_e,1))+r'$\pm$'+str(round(rbberr73_e,1))+r"$R_\oplus$",fontsize=16)
axs.legend()
#fig.savefig("n73models05262.png",bbox_inches='tight')
plt.show()

#%%
fig, ax8 = plt.subplots(1, 1, figsize=(15.2, 10)) #creates plot

cmap = mpl.cm.get_cmap('brg') #choose what color map you want, matplotlib.org/3.2.1/tutorials/colors/colormaps.html

G = 2942.0 #solar radii cubed per solar mass day squared

j=0
mms=[mms358_sol,mms73_sol]
mwd=[0.8,1.0]
for j in range(len(mms)):
    mstar=mms[j]
    msini=mwd[j]
    q_exp = mstar/msini #expected mass ratio, based on the derived mass of the star and mass of the companion
    #msini is WD
    
    period_array = np.arange(100.0,1500.0,2.0) #days, range of periods to plot over, might want to make it span a larger range
    
    a = (period_array**2.0 * (G*(mstar + msini))/(4.0*np.pi**2.0))**(1.0/3.0) #calc semimajor axis
    
    r_RL = (0.49*q_exp**(2.0/3.0))/((0.6*q_exp**(2.0/3.0)) + (np.log(1 + q_exp**(1.0/3.0)))) * a #calc Roche lobe radius at a given period for the expected mass ratio
    
    mstar_array = np.arange(msini,msini*22.0,msini*2.0) #now generate value to do the same calculation for a range of mass ratios
    
    if msini>0.8:
        ax8.fill_between(r_RL,period_array,alpha=0.5,cmap='hsv',clim=(0.05,1.0))
    else:
        ax8.scatter(r_RL, period_array,s=0.4,
                                 c=[1.0/q_exp]*len(r_RL), cmap='brg',
                                 vmin=0.05, vmax=1.0,
                                 norm=mpl.colors.LogNorm())
    ax8.text(230, 460, r'$q_{exp}'+r' = {:.1f}$'.format(1.0/q_exp), rotation=36,
                     ha='left', fontsize=14)
    
    for i,m in enumerate(mstar_array):
            
        q = m/msini #mass ratio
        
        a = (period_array**2.0 * (G*(m + msini))/(4.0*np.pi**2.0))**(1.0/3.0) #semimajor axis
    
        r_RL = (0.49*q**(2.0/3.0))/((0.6*q**(2.0/3.0)) + (np.log(1 + q**(1.0/3.0)))) * a #Roche lobe radius
        
        if j==0:
            if i in [0,1,3,5,7,10]: #only have it plot every few lines, can change these values
                print(q)   
                im = ax8.scatter(r_RL, period_array,
                                 s=0.1,
                                 c=[1.0/q]*len(r_RL), cmap='brg',
                                 vmin=0.05, vmax=1.0,
                                 norm=mpl.colors.LogNorm())

        if j==1:
            if i in [0,1,3,5,7,10]: #only have it plot every few lines, can change these values
                print(q)   
                im = ax8.scatter(r_RL, period_array,
                                 s=0.1,
                                 c=[1.0/q]*len(r_RL), cmap='hsv',
                                 vmin=0.05, vmax=1.0,
                                 norm=mpl.colors.LogNorm())

            

        if i==0: #label the lines that are plotted with their mass ratio
            print(q)
            ax8.text(94, 500, r'$q = {:.1f}$'.format(1.0/q), rotation=66,
                     ha='left', fontsize=14)
            print(r_RL[560])
     
        if i==5:
            print(q)
            ax8.text(265, 450, r'$q = {:.1f}$'.format(1.0/q), rotation=36,
                     ha='left', fontsize=14)
            
        if i==10:
            print(q)
            ax8.text(315, 375, r'$q = {:.2f}$'.format(1.0/q), rotation=27,
                     ha='left', fontsize=14)
    
    rstar=rms358_sol
    rstar_err=rmserr358_sol
    per358=760.0
    
    rstar73=rms73_sol
    
    ax8.errorbar(rstar, per358, xerr=rstar_err, yerr=117., ms=5, marker='o', color='dodgerblue') #add the observed variables to the plot
    ax8.text(rstar, per358+40, r'LIN 358', ha='center', fontsize=16, color='dodgerblue')
    ax8.text(rstar73,440, r'SMC N73',ha='center',fontsize=16,color='dodgerblue')
            
    ax8.set_xlim(50.0,400.0)
    ax8.set_ylim(100.0,1500.0)
    
    ax8.set_xlabel('Radius of RG for RLOF [$R_\odot$]', fontsize=18)
    ax8.set_ylabel('Period [days]', fontsize=18)
    
    if j==0:
        cbar = fig.colorbar(im, ax=ax8, ticks=[0.1, 0.3, 1.0])
        cbar.ax.set_yticklabels(['0.1', '0.3', '1.0'])
        cbar.set_label(r'$q (LIN 358)$', fontsize=18)

    if j==1:
        cbar = fig.colorbar(im, ax=ax8, ticks=[0.1, 0.3, 1.0])
        cbar.ax.set_yticklabels(['0.1', '0.3', '1.0'])
        cbar.set_label(r'$q (SMC N73)$', fontsize=18)
        
    
    plt.tight_layout()
    
    fig.savefig('RL0902.png', dpi=300)