import numpy as np
import pandas as pd
import os
from astropy.io import ascii
import scipy.interpolate as sin

def spec_vot(table,waves):
    uwave=np.unique(waves)
    spec=np.zeros((3,len(uwave)))
    filters=np.array([])
  
    s=0
    for s in np.arange(len(uwave)):
        spec[0][s]=uwave[s]
        spec[1][s]=np.mean(table.at[uwave[s],'flux'])
        spec[2][s]=np.sqrt((np.sum(table.at[uwave[s],'ferror']))**2)
        fils=table.at[uwave[s],'filter']
        filters=np.append(filters,fils)
    
    ufilters=list(dict.fromkeys(filters))
    dat={'wavelength':uwave}
    wavetable=pd.DataFrame(data=dat,index=ufilters,columns=['wavelength'])

    return spec,wavetable

def spec_model(wmod,fmod,wspec):
    spec=np.zeros((2,len(wspec)))
    mod=sin.interp1d(x=wmod,y=fmod,kind='next')
    spec[0]=wspec
    spec[1]=mod(wspec)
    
    return spec,mod
#def fits(start,stop,step,name,Teff_i,logg_const,start_ms,stop_ms,step_ms,dist):

def fits_kur(spec,wavetable):
    chi2comp=2e62
    met=np.array(["m05","m10"])
    logg=np.array(["0.0","0.5"])
    temp=3500 #goes to 6500 in incr of 250
    #took out 1.0 log(g), is there a way to keep it/skip?
    while temp<4600:
        a=0
        for a in np.arange(len(met)):
            b=0
            for b in np.arange(len(logg)):
                filename="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                if os.path.exists(filename)==False: continue
                e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                header1 = e.readline()
                header2 = e.readline()
                header3 = e.readline()
                header4 = e.readline()
                header5 = e.readline()
                header6 = e.readline()
                header7 = e.readline()
                header8 = e.readline()
                header9 = e.readline()
                npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                wavelength_kur=np.zeros(npts_ms)
                flux_kur=np.zeros(npts_ms)
#                band_kur=([])
        

                m=0
                for line in e:
                    line = line.strip()
                    columns = line.split()
                    if len(columns)<3:
                #            band_kur=np.append(band_kur,columns[0])
                        wavelength_kur[m]=0.0
                        flux_kur[m]=0.0
                    else:
                #            band_kur=np.append(band_kur,columns[0])
                        wavelength_kur[m]=np.float(columns[1])
                        flux_kur[m]=np.float(columns[2])
                    m+=1
                
                #    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                #    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                
                minfo=np.zeros((2,len(flux_kur)))
                minfo[0]=wavelength_kur
                minfo[1]=flux_kur
                 
                bands=['2MASS:Ks','WISE:W1','WISE:W2','WISE:W3']
                wlist=list(spec[0])
                
                model_ms,func_ms=spec_model(wavelength_kur,flux_kur,spec[0])
                
                blow=wavetable.at['2MASS:J','wavelength']
                indlower=wlist.index(blow)
                
                wpos=wavelength_kur[flux_kur>0]
                fpos=flux_kur[flux_kur>0]
                mpos=np.zeros((2,len(wpos)))
                mpos[0]=wpos
                mpos[1]=fpos
                
                for b in bands:
                    bw=wavetable.at[b,'wavelength']
                    ind=wlist.index(bw)
                    
                    scale=spec[1][ind]/model_ms[1][ind]
                    smodel_ms=model_ms[1]*scale
                    chi2total=0.
                    for c in np.arange(indlower,len(spec[0])):
                        chi2=((smodel_ms[c]-spec[1][c])/spec[2][c])**2
                        chi2total+=chi2
                #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        model_best=filename
                        scale_best=scale
                        flux_best=model_ms
                        info=b
                        func_best=func_ms

                """
                m=0
                for line in e:
                    line = line.strip()
                    columns = line.split()
                    if len(columns)<3:
                        band_kur=np.append(band_kur,columns[0])
                        wavelength_kur[m]=0.0
                        flux_kur[m]=0.0
                    else:
                        band_kur=np.append(band_kur,columns[0])
                        wavelength_kur[m]=np.float(columns[1])
                        flux_kur[m]=np.float(columns[2])
                    m+=1
                
                dat={'wavelength':wavelength_kur,'flux':flux_kur}
                table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                
                minfo=np.zeros((2,len(flux_kur)))
                minfo[0]=wavelength_kur
                minfo[1]=flux_kur
                
                model_ms=spec_model(table)

                wpos=wavelength_kur[flux_kur>0]
                fpos=flux_kur[flux_kur>0]
                mpos=np.zeros((2,len(wpos)))
                mpos[0]=wpos
                mpos[1]=fpos
                
                scale=spec[1][4]/model_ms[1][4]
                smodel_ms=model_ms[1]*scale
                chi2total=0.
                for c in np.arange(2,7):
                    chi2=((smodel_ms[c]-spec[1][c])/spec[2][c])**2
                    chi2total+=chi2
        #                print(chi2total)
                if chi2total<chi2comp:
                    chi2comp=chi2total
                    model_best=filename
                    scale_best=scale
                    flux_best=mpos
                    info=minfo
                """
        temp+=250
    return chi2comp,model_best,scale_best,flux_best,info,func_best 

def fits_kur_apogee(spec,temp,logg,met,wavetable):
    chi2comp=2e62
    #3800
    #took out 0.0 log(g), is there a way to keep it/skip?
    filename="Kurucz_f"+met+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg+"0000.phot.dat"
#    if os.path.exists(filename)==False: continue
    e=open("Kurucz_f"+met+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg+"0000.phot.dat")
    header1 = e.readline()
    header2 = e.readline()
    header3 = e.readline()
    header4 = e.readline()
    header5 = e.readline()
    header6 = e.readline()
    header7 = e.readline()
    header8 = e.readline()
    header9 = e.readline()
    npts_ms = sum(1 for line in open("Kurucz_f"+met+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg+"0000.phot.dat"))-9
    wavelength_kur=np.zeros(npts_ms)
    flux_kur=np.zeros(npts_ms)
#    band_kur=([])

    m=0
    for line in e:
        line = line.strip()
        columns = line.split()
        if len(columns)<3:
#            band_kur=np.append(band_kur,columns[0])
            wavelength_kur[m]=0.0
            flux_kur[m]=0.0
        else:
#            band_kur=np.append(band_kur,columns[0])
            wavelength_kur[m]=np.float(columns[1])
            flux_kur[m]=np.float(columns[2])
        m+=1
    
#    dat={'wavelength':wavelength_kur,'flux':flux_kur}
#    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
    
    minfo=np.zeros((2,len(flux_kur)))
    minfo[0]=wavelength_kur
    minfo[1]=flux_kur
 
    bands=['2MASS:Ks','WISE:W1','WISE:W2','WISE:W3']
    wlist=list(spec[0])
    
    model_ms,func_ms=spec_model(wavelength_kur,flux_kur,spec[0])

    blow=wavetable.at['2MASS:J','wavelength']
    indlower=wlist.index(blow)

    wpos=wavelength_kur[flux_kur>0]
    fpos=flux_kur[flux_kur>0]
    mpos=np.zeros((2,len(wpos)))
    mpos[0]=wpos
    mpos[1]=fpos
    
    for b in bands:
        bw=wavetable.at[b,'wavelength']
        ind=wlist.index(bw)
        
        scale=spec[1][ind]/model_ms[1][ind]
        smodel_ms=model_ms[1]*scale
        chi2total=0.
        for c in np.arange(indlower,len(spec[0])):
            chi2=((smodel_ms[c]-spec[1][c])/spec[2][c])**2
            chi2total+=chi2
    #                print(chi2total)
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=filename
            scale_best=scale
            flux_best=model_ms
            info=b
            func_best=func_ms

    return chi2comp,model_best,scale_best,flux_best,info,func_best        

def fits_bb(spec,wavetable): 
    chi2comp=2e62
    temp=100000
    dsmc=60.6*1000 #pc
    dsmc_cm=dsmc*(3.086e18) #cm
    slim=((0.4*(6.96e10))/dsmc_cm)**2
    while temp<200050:
        filename="bbody_bb"+str(temp)+".phot.dat"
        if os.path.exists(filename)==False: continue
        g = open("bbody_bb"+str(temp)+".phot.dat")
        g.readline()
        g.readline()
        g.readline()
        g.readline()
        g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp)+".phot.dat"))-5
#        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
#                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
#                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
#        dat={'wavelength':wavelength_bb,'flux':flux_bb}
#        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
#        minfo=np.zeros((2,len(flux_bb)))
#        minfo[0]=wavelength_bb
#        minfo[1]=flux_bb
        
        bands=['GALEX:FUV','GALEX:NUV']
        wlist=list(spec[0])
        skip=['XMM-OT:U']
                            
        wpos=wavelength_bb[flux_bb>0]
        fpos=flux_bb[flux_bb>0]
        
        sortwave=np.sort(wpos)
        
        dat={'flux':fpos}
        table=pd.DataFrame(data=dat,index=wpos,columns=['flux'])
        
        sortflux=np.zeros(len(sortwave))
        
        u=0
        for u in np.arange(len(sortwave)):
            sortflux[u]=np.mean(table.at[sortwave[u],'flux'])
                
        x=sortwave
        y=sortflux
        
        z=0
        z1=1
        
        stop=3000
        for z in np.arange(stop):
            dif=abs(y[z]-y[z1])
            if dif>(2.4*(y[z]/10)):
                y[z1]=y[z]
            z1+=1
            if z1==stop: break

        mpos=np.zeros((2,len(x)))
        mpos[0]=x
        mpos[1]=y

#        mpos=np.zeros((2,len(wpos)))
#        mpos[0]=wpos
#       mpos[1]=fpos

#        model_bb=spec_model(wavelength_bb,flux_bb,spec[0])
        model_bb,func_bb=spec_model(x,y,spec[0])

        bup=wavetable.at['Johnson:B','wavelength']
#        indup=wlist.index(bup)
        indup=2
#        print("INDUP",indup)        
#        indup=[4,5,6,7,8,9,10]
        
        for b in bands:
            bw=wavetable.at[b,'wavelength']
            ind=wlist.index(bw)
            scale=spec[1][ind]/model_bb[1][ind]
            if scale>slim:continue
            smodel_bb=model_bb[1]*scale
#            print(smodel_bb[1])
            chi2total=0.
            for c in np.arange(0,len(spec[0])):
                if spec[0][c]>3600: break
                if abs(spec[0][c]-3443.0)<1.0: 
#                    print("got em!")
                    continue
#                if spec[1][c]<=0.0: continue
                chi2=((smodel_bb[c]-spec[1][c])/(spec[2][c]))**2
#                print("CHI2",chi2)
                chi2total+=chi2
#                print(c)
#            print("TOTAL",chi2total)
            if chi2total<chi2comp:
                chi2comp=chi2total
                model_best=filename
                scale_best=scale
                flux_best=mpos
                info=b
                func_best=func_bb
                ind_best=indup
                                 
        temp+=1000

    return chi2comp,model_best,scale_best,flux_best,info,func_best,ind_best       


def kurucz_getter(file):
    e = open(file+".dat.txt",'r')
    header1 = e.readline()
    header2 = e.readline()
    header3 = e.readline()
    header4 = e.readline()
    header5 = e.readline()
    header6 = e.readline()
    header7 = e.readline()
    header8 = e.readline()
    header9 = e.readline()
    npts_ms = sum(1 for line in open(file+".dat.txt"))-9
    wavelength_kur=np.zeros(npts_ms)
    flux_kur=np.zeros(npts_ms)
    
    m=0
    for line in e:
        line = line.strip()
        columns = line.split()
        wavelength_kur[m]=np.float(columns[0])
        flux_kur[m]=np.float(columns[1])
        m+=1
    return wavelength_kur,flux_kur

def blackbody(filename,spec):
    chi2comp=2e63
    g = open(filename)
    header1_2 = g.readline()
    header2_2 = g.readline()
    header3_2 = g.readline()
    header4_2 = g.readline()
    header5_2 = g.readline()

    npts_wd = sum(1 for line in open(filename))-5
    band_bb=([])
    flux_bb=np.zeros(npts_wd)
    wavelength_bb=np.zeros(npts_wd)
    
    q=0
    for line in g:
        line = line.strip()
        columns = line.split()
        if len(columns)<3:
            band_bb=np.append(band_bb,columns[0])
            wavelength_bb[q]=0.0
            flux_bb[q]=0.0
        else:
            band_bb=np.append(band_bb,columns[0])
            wavelength_bb[q]=columns[1]
            flux_bb[q]=columns[2]
        q+=1
        
    dat={'wavelength':wavelength_bb,'flux':flux_bb}
    table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
        
    minfo=np.zeros((2,len(flux_bb)))
    minfo[0]=wavelength_bb
    minfo[1]=flux_bb
    
    model_bb=spec_model(table)
    
    wpos=wavelength_bb[flux_bb>0]
    fpos=flux_bb[flux_bb>0]
    mpos=np.zeros((2,len(wpos)))
    mpos[0]=wpos
    mpos[1]=fpos
    
    for b in np.arange(0,1):
        scale=spec[1][b]/model_bb[1][b]
        smodel_bb=model_bb[1]*scale
        chi2total=0.
        for c in np.arange(0,4):
            chi2=(abs(smodel_bb[c]-spec[1][c]))**2
            chi2total+=chi2
        if chi2total<chi2comp:
            scale_best=scale

#    smodel_bb=model_bb[1]*scale

    return mpos,scale_best


def bb(temp,wavelength):
    j=(5.6704e-5)*(temp**4)
    i=0
    fbb=np.zeros(len(wavelength))
    for i in np.arange(len(wavelength)):
        fbb[i]=j/wavelength[i]
    return fbb

def fits_both(spec,met,temp,logg):
    chi2comp=2e62
    temp_bb=100000
    while temp_bb<200050:
#        print("BB:",str(temp_bb))
        filename_bb="bbody_bb"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table)

        wpos=wavelength_bb[flux_bb>0]
        fpos=flux_bb[flux_bb>0]
        mpos=np.zeros((2,len(wpos)))
        mpos[0]=wpos
        mpos[1]=fpos

            
        filename_kur="Kurucz_f"+met+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg+"0000.phot.dat"
        if os.path.exists(filename_kur)==False: continue
        e=open("Kurucz_f"+met+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg+"0000.phot.dat")
        header1 = e.readline()
        header2 = e.readline()
        header3 = e.readline()
        header4 = e.readline()
        header5 = e.readline()
        header6 = e.readline()
        header7 = e.readline()
        header8 = e.readline()
        header9 = e.readline()
        npts_ms = sum(1 for line in open("Kurucz_f"+met+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg+"0000.phot.dat"))-9
        wavelength_kur=np.zeros(npts_ms)
        flux_kur=np.zeros(npts_ms)
        band_kur=([])

        m=0
        for line in e:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_kur=np.append(band_kur,columns[0])
                wavelength_kur[m]=0.0
                flux_kur[m]=0.0
            else:
                band_kur=np.append(band_kur,columns[0])
                wavelength_kur[m]=np.float(columns[1])
                flux_kur[m]=np.float(columns[2])
            m+=1
        
        dat={'wavelength':wavelength_kur,'flux':flux_kur}
        table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
        
        model_kur=spec_model(table)
        
        model_both=model_kur+model_bb
        
        for i in np.arange(0,len(model_both)):
            scale=spec[1][4]/model_both[1][4]
            smodel=model_both[1]*scale
            chi2total=0.
            for c in np.arange(0,len(smodel)):
                chi2=(abs(smodel[c]-spec[1][c]))**2
                chi2total+=chi2
    #                print(chi2total)
            if chi2total<chi2comp:
                chi2comp=chi2total
                model_best_bb=filename_bb
                model_best_kur=filename_kur
                scale_best=scale
                flux_best_kur=model_kur
                flux_best_bb=mpos
                flux_best_both=model_both
                band=band_kur[i]
        temp_bb+=1000
        
    return chi2comp,model_best_kur,model_best_bb,scale_best,flux_best_kur,flux_best_bb,flux_best_both,band


def fits_temp(spec):
    chi2comp=2e62
    temp_bb=100000
    while temp_bb<200050:
        print("BB",str(temp_bb))
        filename_bb="bbody_bb"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table)
            
        met=np.array(["p00","p02","p05","m05","m10","m15","m20","m25"])
        logg=np.array(["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"])
        temp=3500 #goes to 6500 in incr of 250
        while temp<4600:
            a=0
            for a in np.arange(len(met)):
                b=0
                for b in np.arange(len(logg)):
                    filename_kur="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                    if os.path.exists(filename_kur)==False: continue
                    e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                    header1 = e.readline()
                    header2 = e.readline()
                    header3 = e.readline()
                    header4 = e.readline()
                    header5 = e.readline()
                    header6 = e.readline()
                    header7 = e.readline()
                    header8 = e.readline()
                    header9 = e.readline()
                    npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                    wavelength_kur=np.zeros(npts_ms)
                    flux_kur=np.zeros(npts_ms)
                    band_kur=([])
            
                    m=0
                    for line in e:
                        line = line.strip()
                        columns = line.split()
                        if len(columns)<3:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=0.0
                            flux_kur[m]=0.0
                        else:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=np.float(columns[1])
                            flux_kur[m]=np.float(columns[2])
                        m+=1
                    
                    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                    
                    model_kur=spec_model(table)
                    
                    model_both=model_kur+model_bb
                    scale=spec[1][4]/model_both[1][4]
                    smodel=model_both[1]*scale
                    chi2total=0.
                    for c in np.arange(0,len(smodel)):
                        chi2=(abs(smodel[c]-spec[1][c]))**2
                        chi2total+=chi2
            #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        temp_mid=temp_bb
            print("KUR: ",str(temp))
            print(chi2comp)
            temp+=250    
        temp_bb+=10000
    return temp_mid

def fits_both_specific(spec,start,stop,step):
    chi2comp=2e62
    temp_bb=start
    while temp_bb<stop:
        print("BB: ",str(temp_bb))
        filename_bb="bbody_bb"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table_bb=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table_bb)
        
        met=np.array(["p00","p02","p05","m05","m10","m15","m20","m25"])
        logg=np.array(["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"])
        temp=3500 #goes to 6500 in incr of 250
        while temp<4600:
            a=0
            for a in np.arange(len(met)):
                b=0
                for b in np.arange(len(logg)):
                    filename_kur="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                    if os.path.exists(filename_kur)==False: continue
                    e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                    header1 = e.readline()
                    header2 = e.readline()
                    header3 = e.readline()
                    header4 = e.readline()
                    header5 = e.readline()
                    header6 = e.readline()
                    header7 = e.readline()
                    header8 = e.readline()
                    header9 = e.readline()
                    npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                    wavelength_kur=np.zeros(npts_ms)
                    flux_kur=np.zeros(npts_ms)
                    band_kur=([])
            
                    m=0
                    for line in e:
                        line = line.strip()
                        columns = line.split()
                        if len(columns)<3:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=0.0
                            flux_kur[m]=0.0
                        else:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=np.float(columns[1])
                            flux_kur[m]=np.float(columns[2])
                        m+=1
                    
                    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                    
                    model_kur=spec_model(table)
                    
                    model_both=model_kur+model_bb
                    scale=spec[1][4]/model_both[1][4]
                    smodel=model_both[1]*scale
                    chi2total=0.
                    for c in np.arange(0,len(smodel)):
                        chi2=(abs(smodel[c]-spec[1][c]))**2
                        chi2total+=chi2
            #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        model_best_bb=filename_bb
                        model_best_kur=filename_kur
                        scale_best=scale
                        flux_best_kur=model_kur
                        flux_best_bb=model_bb
                        flux_best_both=model_both
            print("KUR: ",str(temp))
            print(chi2comp)
            temp+=250    
        temp_bb+=step

    return chi2comp,model_best_kur,model_best_bb,scale_best,flux_best_kur,flux_best_bb,flux_best_both

def fits_bb8(spec): 
    chi2comp=2e62
    temp=100000
    while temp<200050:
        filename="bbody_bb"+str(temp)+".phot.dat"
        if os.path.exists(filename)==False: continue
        g = open("bbody_bb"+str(temp)+".phot.dat")
        
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp)+".phot.dat"))-5
        band_bb=([])
#        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
#                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
#                flux_bb[q]=columns[2]
            q+=1
                
        tbb=np.float(filename[8:14])
        flux_bb=bb(tbb,wavelength_bb)
        
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
#        minfo=np.zeros((2,len(flux_bb)))
#        minfo[0]=wavelength_bb
#        minfo[1]=flux_bb
        
        model_bb=spec_model(table)
        
        wpos=wavelength_bb[flux_bb>0]
        fpos=flux_bb[flux_bb>0]
        mpos=np.zeros((2,len(wpos)))
        mpos[0]=wpos
        mpos[1]=fpos
            
        scale=spec[1][0]/model_bb[1][0]
        smodel_bb=model_bb[1]*scale
        chi2total=0.
        for c in np.arange(0,4):
            chi2=(abs(smodel_bb[c]-spec[1][c]))**2
            chi2total+=chi2
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=filename
            scale_best=scale
            flux_best=mpos
#            info=minfo                 
        temp+=1000

    return chi2comp,model_best,scale_best,flux_best        

#fits_both comment out

"""
    temp_bb=5000
    while temp_bb<9999:
        print("BB: ",str(temp_bb))
        filename_bb="bbody_bb0"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb0"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb0"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table_bb=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table_bb)
        
        met=np.array(["p00","p02","p05","m05","m10","m15","m20","m25"])
        logg=np.array(["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"])
        temp=3500 #goes to 6500 in incr of 250
        while temp<6600:
            a=0
            for a in np.arange(len(met)):
                b=0
                for b in np.arange(len(logg)):
                    filename_kur="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                    if os.path.exists(filename_kur)==False: continue
                    e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                    header1 = e.readline()
                    header2 = e.readline()
                    header3 = e.readline()
                    header4 = e.readline()
                    header5 = e.readline()
                    header6 = e.readline()
                    header7 = e.readline()
                    header8 = e.readline()
                    header9 = e.readline()
                    npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                    wavelength_kur=np.zeros(npts_ms)
                    flux_kur=np.zeros(npts_ms)
                    band_kur=([])
            
                    m=0
                    for line in e:
                        line = line.strip()
                        columns = line.split()
                        if len(columns)<3:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=0.0
                            flux_kur[m]=0.0
                        else:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=np.float(columns[1])
                            flux_kur[m]=np.float(columns[2])
                        m+=1
                    
                    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                    
                    model_kur=spec_model(table)
                    
                    model_both=model_kur+model_bb
                    scale=spec[1][4]/model_both[1][4]
                    smodel=model_both[1]*scale
                    chi2total=0.
                    for c in np.arange(0,len(smodel)):
                        chi2=(abs(smodel[c]-spec[1][c]))**2
                        chi2total+=chi2
            #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        model_best_bb=filename_bb
                        model_best_kur=filename_kur
                        scale_best=scale
                        flux_best_kur=model_kur
                        flux_best_bb=model_bb
                        flux_best_both=model_both
            print("KUR: ",str(temp))
            print(chi2comp)
            temp+=250    
        temp_bb+=50
    print("5000 cleared")
    temp_bb=10000
    while temp_bb<20020:
        filename_bb="bbody_bb"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table)
            
        met=np.array(["p00","p02","p05","m05","m10","m15","m20","m25"])
        logg=np.array(["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"])
        temp=3500 #goes to 6500 in incr of 250
        while temp<4600:
            a=0
            for a in np.arange(len(met)):
                b=0
                for b in np.arange(len(logg)):
                    filename_kur="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                    if os.path.exists(filename_kur)==False: continue
                    e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                    header1 = e.readline()
                    header2 = e.readline()
                    header3 = e.readline()
                    header4 = e.readline()
                    header5 = e.readline()
                    header6 = e.readline()
                    header7 = e.readline()
                    header8 = e.readline()
                    header9 = e.readline()
                    npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                    wavelength_kur=np.zeros(npts_ms)
                    flux_kur=np.zeros(npts_ms)
                    band_kur=([])
            
                    m=0
                    for line in e:
                        line = line.strip()
                        columns = line.split()
                        if len(columns)<3:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=0.0
                            flux_kur[m]=0.0
                        else:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=np.float(columns[1])
                            flux_kur[m]=np.float(columns[2])
                        m+=1
                    
                    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                    
                    model_kur=spec_model(table)
                    
                    model_both=model_kur+model_bb
                    scale=spec[1][4]/model_both[1][4]
                    smodel=model_both[1]*scale
                    chi2total=0.
                    for c in np.arange(0,len(smodel)):
                        chi2=(abs(smodel[c]-spec[1][c]))**2
                        chi2total+=chi2
            #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        model_best_bb=filename_bb
                        model_best_kur=filename_kur
                        scale_best=scale
                        flux_best_kur=model_kur
                        flux_best_bb=model_bb
                        flux_best_both=model_both
            print(chi2comp)
            temp+=250    
        temp_bb+=50
    print("10000 cleared")

    temp_bb=20000
"""


###############################################################################
"""
def fits_wd_tophat(name):
    chi2comp=2e62
    
    d = open(name+".dat",'r')
    header = d.readline()
    npts = sum(1 for line in open(name+".dat"))-1
    
    band=np.array([])
    wavelength=np.zeros(npts)
    flux=np.zeros(npts)
    ferror=np.zeros(npts)
    
    l=0
    for line in d:
        line = line.strip()
        columns = line.split()
        band=np.append(band,columns[0])
        wavelength[l]=np.float(columns[1])
        flux[l]=np.float(columns[4])
        ferror[l]=np.float(columns[5])
        l = l+1
    
    dat={'wavelength':wavelength,'flux':flux}
    table=pd.DataFrame(data=dat,index=band,columns=['wavelength','flux'])
    
    spec=np.zeros((10,2))
    spec[0][1]=table.at['GALEX/GALEX.FUV','flux']
    spec[1][1]=table.at['GALEX/GALEX.NUV','flux']
    spec[2][1]=table.at['GAIA/GAIA2.G','flux']
    spec[3][1]=table.at['2MASS/2MASS.J','flux']
    spec[4][1]=table.at['2MASS/2MASS.H','flux']
    spec[5][1]=table.at['2MASS/2MASS.Ks','flux']
    spec[6][1]=table.at['WISE/WISE.W1','flux']
    spec[7][1]=table.at['WISE/WISE.W2','flux']
    spec[8][1]=table.at['WISE/WISE.W3','flux']
    spec[9][1]=table.at['WISE/WISE.W4','flux']

    spec[0][0]=table.at['GALEX/GALEX.FUV','wavelength']
    spec[1][0]=table.at['GALEX/GALEX.NUV','wavelength']
    spec[2][0]=table.at['GAIA/GAIA2.G','wavelength']
    spec[3][0]=table.at['2MASS/2MASS.J','wavelength']
    spec[4][0]=table.at['2MASS/2MASS.H','wavelength']
    spec[5][0]=table.at['2MASS/2MASS.Ks','wavelength']
    spec[6][0]=table.at['WISE/WISE.W1','wavelength']
    spec[7][0]=table.at['WISE/WISE.W2','wavelength']
    spec[8][0]=table.at['WISE/WISE.W3','wavelength']
    spec[9][0]=table.at['WISE/WISE.W4','wavelength']

    modwd,modelname=flux_models("models_flux_filters_KOESTER.dat")
    
    filters=np.array(['GALEX_FUV', 'GALEX_NUV', 'GAIA2_G', 'SDSS_i', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4'])

    jasmin=0
    for jasmin in range(len(modelname)):
        g = open(modelname[jasmin],'r')

        q=0
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
        header6_2 = g.readline()
    
        npts_wd = sum(1 for line in open(modelname[jasmin]))-6
        flux_koe=np.zeros(npts_wd)
        wavelength_koe=np.zeros(npts_wd)
        for line in g:
            line = line.strip()
            columns = line.split()
            wavelength_koe[q]=columns[0]
            flux_koe[q]=columns[1]
            q+=1

        model_wd=tophat(wavelength_koe,flux_koe)
        nans=np.isnan(model_wd)
        model_wd[nans]=0
    
        scale=spec[0][1]/model_wd[0]
        smodel_wd=model_wd*scale
        chi2total=0.
        for c in np.arange(0,4):
            chi2=(abs(smodel_wd[c]-spec[c][1]))**2
            chi2total+=chi2
#               print(chi2total)
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=modelname[jasmin]
            scale_best=scale
            flux_best=smodel_wd 

    return chi2comp,model_best,scale_best,flux_best             

def fits_ms_tophat(name):
    chi2comp=2e62
    
    d = open(name+".dat",'r')
    header = d.readline()
    npts = sum(1 for line in open(name+".dat"))-1
    
    band=np.array([])
    wavelength=np.zeros(npts)
    flux=np.zeros(npts)
    ferror=np.zeros(npts)
    
    l=0
    for line in d:
        line = line.strip()
        columns = line.split()
        band=np.append(band,columns[0])
        wavelength[l]=np.float(columns[1])
        flux[l]=np.float(columns[4])
        ferror[l]=np.float(columns[5])
        l = l+1
    
    dat={'wavelength':wavelength,'flux':flux}
    table=pd.DataFrame(data=dat,index=band,columns=['wavelength','flux'])
    
    spec=np.zeros((10,2))
    spec[0][1]=table.at['GALEX/GALEX.FUV','flux']
    spec[1][1]=table.at['GALEX/GALEX.NUV','flux']
    spec[2][1]=table.at['GAIA/GAIA2.G','flux']
    spec[3][1]=table.at['2MASS/2MASS.J','flux']
    spec[4][1]=table.at['2MASS/2MASS.H','flux']
    spec[5][1]=table.at['2MASS/2MASS.Ks','flux']
    spec[6][1]=table.at['WISE/WISE.W1','flux']
    spec[7][1]=table.at['WISE/WISE.W2','flux']
    spec[8][1]=table.at['WISE/WISE.W3','flux']
    spec[9][1]=table.at['WISE/WISE.W4','flux']

    spec[0][0]=table.at['GALEX/GALEX.FUV','wavelength']
    spec[1][0]=table.at['GALEX/GALEX.NUV','wavelength']
    spec[2][0]=table.at['GAIA/GAIA2.G','wavelength']
    spec[3][0]=table.at['2MASS/2MASS.J','wavelength']
    spec[4][0]=table.at['2MASS/2MASS.H','wavelength']
    spec[5][0]=table.at['2MASS/2MASS.Ks','wavelength']
    spec[6][0]=table.at['WISE/WISE.W1','wavelength']
    spec[7][0]=table.at['WISE/WISE.W2','wavelength']
    spec[8][0]=table.at['WISE/WISE.W3','wavelength']
    spec[9][0]=table.at['WISE/WISE.W4','wavelength']

    modms,modelname=flux_models("models_flux_filters_KURUCZ.dat")
    
#    filters=np.array(['GALEX_FUV', 'GALEX_NUV', 'GAIA2_G', 'SDSS_i', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4'])
    
    morgan=0
    for morgan in range(len(modelname)):
        e = open(modelname[morgan]+"..logg=4.00000.dat.txt",'r')
        header1 = e.readline()
        header2 = e.readline()
        header3 = e.readline()
        header4 = e.readline()
        header5 = e.readline()
        header6 = e.readline()
        header7 = e.readline()
        header8 = e.readline()
        header9 = e.readline()
        npts_ms = sum(1 for line in open(modelname[morgan]+"..logg=4.00000.dat.txt"))-9
        wavelength_kur=np.zeros(npts_ms)
        flux_kur=np.zeros(npts_ms)
        
        m=0
        for line in e:
            line = line.strip()
            columns = line.split()
            wavelength_kur[m]=np.float(columns[0])
            flux_kur[m]=np.float(columns[1])
            m+=1
            
        model_ms=tophat(wavelength_kur,flux_kur)
        nans=np.isnan(model_ms)
        model_ms[nans]=0

        scale=spec[4][1]/model_ms[4]
        smodel_ms=model_ms*scale
        chi2total=0.
        for c in np.arange(2,7):
            chi2=(abs(smodel_ms[c]-spec[c][1]))**2
            chi2total+=chi2
#                print(chi2total)
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=modelname[morgan]
            scale_best=scale
            flux_best=smodel_ms 

    return chi2comp,model_best,scale_best,flux_best        


#    n=0
#    for n in np.arange(len(modelname)):
#        model_ms=np.zeros(10,dtype=float)
#        f=0
#        for f in np.arange(len(model_ms)):
#            model_ms[f]=np.mean(modms.at[modelname[n],filters[f]])


def fits_both_tophat(name):
    chi2comp=2e62
    
    d = open(name+".dat",'r')
    header = d.readline()
    npts = sum(1 for line in open(name+".dat"))-1
    
    band=np.array([])
    wavelength=np.zeros(npts)
    flux=np.zeros(npts)
    ferror=np.zeros(npts)
    
    l=0
    for line in d:
        line = line.strip()
        columns = line.split()
        band=np.append(band,columns[0])
        wavelength[l]=np.float(columns[1])
        flux[l]=np.float(columns[4])
        ferror[l]=np.float(columns[5])
        l = l+1
    
    dat={'wavelength':wavelength,'flux':flux}
    table=pd.DataFrame(data=dat,index=band,columns=['wavelength','flux'])
    
    spec=np.zeros((10,2))
    spec[0][1]=table.at['GALEX/GALEX.FUV','flux']
    spec[1][1]=table.at['GALEX/GALEX.NUV','flux']
    spec[2][1]=table.at['GAIA/GAIA2.G','flux']
    spec[3][1]=table.at['2MASS/2MASS.J','flux']
    spec[4][1]=table.at['2MASS/2MASS.H','flux']
    spec[5][1]=table.at['2MASS/2MASS.Ks','flux']
    spec[6][1]=table.at['WISE/WISE.W1','flux']
    spec[7][1]=table.at['WISE/WISE.W2','flux']
    spec[8][1]=table.at['WISE/WISE.W3','flux']
    spec[9][1]=table.at['WISE/WISE.W4','flux']

    spec[0][0]=table.at['GALEX/GALEX.FUV','wavelength']
    spec[1][0]=table.at['GALEX/GALEX.NUV','wavelength']
    spec[2][0]=table.at['GAIA/GAIA2.G','wavelength']
    spec[3][0]=table.at['2MASS/2MASS.J','wavelength']
    spec[4][0]=table.at['2MASS/2MASS.H','wavelength']
    spec[5][0]=table.at['2MASS/2MASS.Ks','wavelength']
    spec[6][0]=table.at['WISE/WISE.W1','wavelength']
    spec[7][0]=table.at['WISE/WISE.W2','wavelength']
    spec[8][0]=table.at['WISE/WISE.W3','wavelength']
    spec[9][0]=table.at['WISE/WISE.W4','wavelength']

    modms,modelname_ms=flux_models("models_flux_filters_KURUCZ.dat")
    
#    filters=np.array(['GALEX_FUV', 'GALEX_NUV', 'GAIA2_G', 'SDSS_i', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4'])
    
    modwd,modelname_wd=flux_models("models_flux_filters_KOESTER.dat")
    
    j=0
    for j in range(len(modelname_ms)):
        e = open(modelname_ms[j]+"..logg=4.00000.dat.txt",'r')
        header1 = e.readline()
        header2 = e.readline()
        header3 = e.readline()
        header4 = e.readline()
        header5 = e.readline()
        header6 = e.readline()
        header7 = e.readline()
        header8 = e.readline()
        header9 = e.readline()
        npts_ms = sum(1 for line in open(modelname_ms[j]+"..logg=4.00000.dat.txt"))-9
        wavelength_kur=np.zeros(npts_ms)
        flux_kur=np.zeros(npts_ms)
        
        m=0
        for line in e:
            line = line.strip()
            columns = line.split()
            wavelength_kur[m]=np.float(columns[0])
            flux_kur[m]=np.float(columns[1])
            m+=1
            
        model_ms=tophat(wavelength_kur,flux_kur)
        nans=np.isnan(model_ms)
        model_ms[nans]=0

        scale_ms=spec[4][1]/model_ms[4]
        smodel_ms=model_ms*scale_ms

        g = open(modelname_wd[j],'r')

        q=0
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
        header6_2 = g.readline()
    
        npts_wd = sum(1 for line in open(modelname_wd[j]))-6
        flux_koe=np.zeros(npts_wd)
        wavelength_koe=np.zeros(npts_wd)
        for line in g:
            line = line.strip()
            columns = line.split()
            wavelength_koe[q]=columns[0]
            flux_koe[q]=columns[1]
            q+=1

        model_wd=tophat(wavelength_koe,flux_koe)
        nans=np.isnan(model_wd)
        model_wd[nans]=0
    
        scale_wd=spec[0][1]/model_wd[0]
        smodel_wd=model_wd*scale_wd
        smodel=smodel_wd+smodel_ms
        chi2total=0.
        for c in np.arange(len(smodel)):
            chi2=(abs(smodel[c]-spec[c][1]))**2
            chi2total+=chi2
#               print(chi2total)
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best_ms=modelname_ms[j]
            model_best_wd=modelname_ms[j]
            scale_best_ms=scale_ms
            scale_best_wd=scale_wd
            flux_best=smodel

    return chi2comp,model_best_ms,model_best_wd,scale_best_ms,scale_best_wd,flux_best             


##################################

def fits_ms(name):
    chi2comp=2e62
    
    d = open(name+".dat",'r')
    header = d.readline()
    npts = sum(1 for line in open(name+".dat"))-1
    
    band=np.array([])
    wavelength=np.zeros(npts)
    flux=np.zeros(npts)
    ferror=np.zeros(npts)
    
    l=0
    for line in d:
        line = line.strip()
        columns = line.split()
        band=np.append(band,columns[0])
        wavelength[l]=np.float(columns[1])
        flux[l]=np.float(columns[4])
        ferror[l]=np.float(columns[5])
        l = l+1
    
    dat={'wavelength':wavelength,'flux':flux}
    table=pd.DataFrame(data=dat,index=band,columns=['wavelength','flux'])
    
    spec=np.zeros((10,2))
    spec[0][1]=table.at['GALEX/GALEX.FUV','flux']
    spec[1][1]=table.at['GALEX/GALEX.NUV','flux']
    spec[2][1]=table.at['GAIA/GAIA2.G','flux']
    spec[3][1]=table.at['2MASS/2MASS.J','flux']
    spec[4][1]=table.at['2MASS/2MASS.H','flux']
    spec[5][1]=table.at['2MASS/2MASS.Ks','flux']
    spec[6][1]=table.at['WISE/WISE.W1','flux']
    spec[7][1]=table.at['WISE/WISE.W2','flux']
    spec[8][1]=table.at['WISE/WISE.W3','flux']
    spec[9][1]=table.at['WISE/WISE.W4','flux']

    spec[0][0]=table.at['GALEX/GALEX.FUV','wavelength']
    spec[1][0]=table.at['GALEX/GALEX.NUV','wavelength']
    spec[2][0]=table.at['GAIA/GAIA2.G','wavelength']
    spec[3][0]=table.at['2MASS/2MASS.J','wavelength']
    spec[4][0]=table.at['2MASS/2MASS.H','wavelength']
    spec[5][0]=table.at['2MASS/2MASS.Ks','wavelength']
    spec[6][0]=table.at['WISE/WISE.W1','wavelength']
    spec[7][0]=table.at['WISE/WISE.W2','wavelength']
    spec[8][0]=table.at['WISE/WISE.W3','wavelength']
    spec[9][0]=table.at['WISE/WISE.W4','wavelength']

    modms,modelname=flux_models("models_flux_filters_KURUCZ.dat")
    
    filters=np.array(['GALEX_FUV', 'GALEX_NUV', 'GAIA2_G', 'SDSS_i', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4'])
    
    n=0
    for n in np.arange(len(modelname)):
        model_ms=np.zeros(10,dtype=float)
        f=0
        for f in np.arange(len(model_ms)):
            model_ms[f]=np.mean(modms.at[modelname[n],filters[f]])
        scale=spec[4][1]/model_ms[4]
        smodel_ms=model_ms*scale
        chi2total=0.
        for c in np.arange(2,7):
            chi2=(abs(smodel_ms[c]-spec[c][1]))**2
            chi2total+=chi2
#                print(chi2total)
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=modelname[n]
            scale_best=scale
            flux_best=smodel_ms 

    return chi2comp,model_best,scale_best,flux_best        
    

def fits_wd(name):
    chi2comp=2e62
    
    d = open(name+".dat",'r')
    header = d.readline()
    npts = sum(1 for line in open(name+".dat"))-1
    
    band=np.array([])
    wavelength=np.zeros(npts)
    flux=np.zeros(npts)
    ferror=np.zeros(npts)
    
    l=0
    for line in d:
        line = line.strip()
        columns = line.split()
        band=np.append(band,columns[0])
        wavelength[l]=np.float(columns[1])
        flux[l]=np.float(columns[4])
        ferror[l]=np.float(columns[5])
        l = l+1
    
    dat={'wavelength':wavelength,'flux':flux}
    table=pd.DataFrame(data=dat,index=band,columns=['wavelength','flux'])
    
    spec=np.zeros((10,2))
    spec[0][1]=table.at['GALEX/GALEX.FUV','flux']
    spec[1][1]=table.at['GALEX/GALEX.NUV','flux']
    spec[2][1]=table.at['GAIA/GAIA2.G','flux']
    spec[3][1]=table.at['2MASS/2MASS.J','flux']
    spec[4][1]=table.at['2MASS/2MASS.H','flux']
    spec[5][1]=table.at['2MASS/2MASS.Ks','flux']
    spec[6][1]=table.at['WISE/WISE.W1','flux']
    spec[7][1]=table.at['WISE/WISE.W2','flux']
    spec[8][1]=table.at['WISE/WISE.W3','flux']
    spec[9][1]=table.at['WISE/WISE.W4','flux']

    spec[0][0]=table.at['GALEX/GALEX.FUV','wavelength']
    spec[1][0]=table.at['GALEX/GALEX.NUV','wavelength']
    spec[2][0]=table.at['GAIA/GAIA2.G','wavelength']
    spec[3][0]=table.at['2MASS/2MASS.J','wavelength']
    spec[4][0]=table.at['2MASS/2MASS.H','wavelength']
    spec[5][0]=table.at['2MASS/2MASS.Ks','wavelength']
    spec[6][0]=table.at['WISE/WISE.W1','wavelength']
    spec[7][0]=table.at['WISE/WISE.W2','wavelength']
    spec[8][0]=table.at['WISE/WISE.W3','wavelength']
    spec[9][0]=table.at['WISE/WISE.W4','wavelength']

    modwd,modelname=flux_models("models_flux_filters_KOESTER.dat")
    
    filters=np.array(['GALEX_FUV', 'GALEX_NUV', 'GAIA2_G', 'SDSS_i', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4'])
    
    n=0
    for n in np.arange(len(modelname)):
        model_wd=np.zeros(10,dtype=float)
        f=0
        for f in np.arange(len(model_wd)):
            model_wd[f]=np.mean(modwd.at[modelname[n],filters[f]])
        scale=spec[1][1]/model_wd[1]
        smodel_wd=model_wd*scale
        chi2total=0.
        for c in np.arange(0,2):
            chi2=(abs(smodel_wd[c]-spec[c][1]))**2
            chi2total+=chi2
#                print(chi2total)
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=modelname[n]
            scale_best=scale
            flux_best=smodel_wd 

    return chi2comp,model_best,scale_best,flux_best   

def flux_models(filename):
#    f=open(filename,'r')
#    header1=f.readline()
#    header2=f.readline()

    table=ascii.read(filename,data_start=2,delimiter=' ')
#    print(table)
#    table['col2'][0]
    npts=len(table['col1'])

#    floatVar = float(table['col2'][0])
#    print(floatVar)
  
    modelname=np.array([])
    GALEX_FUV=np.zeros((npts),dtype=float)
    GALEX_NUV=np.zeros((npts),dtype=float)
    SDSS_u=np.zeros((npts),dtype=float)
    APASS_B=np.zeros((npts),dtype=float)
    SDSS_g=np.zeros((npts),dtype=float)
    GAIA2_Gbp=np.zeros((npts),dtype=float)
    APASS_V=np.zeros((npts),dtype=float)
    SDSS_r=np.zeros((npts),dtype=float)
    GAIA2_G=np.zeros((npts),dtype=float)
    SDSS_i=np.zeros((npts),dtype=float)
    GAIA2_Grp=np.zeros((npts),dtype=float)
    twoMASS_J=np.zeros((npts),dtype=float)
    twoMASS_H=np.zeros((npts),dtype=float)
    twoMASS_Ks=np.zeros((npts),dtype=float)
    WISE_W1=np.zeros((npts),dtype=float)
    WISE_W2=np.zeros((npts),dtype=float)
    WISE_W3=np.zeros((npts),dtype=float)
    WISE_W4 = np.zeros(npts,dtype=float)

    g=0
#    for line in f:
    for g in np.arange(len(table['col1'])):
#        line = f.readline()
#        columns = line.split(' ')
#        modelname=np.append(modelname,columns[0])
#        print(columns[1])
        modelname=np.append(modelname,table['col1'][g])
        GALEX_FUV[g]=float(table['col2'][g])
        GALEX_NUV[g]=float(table['col3'][g])
        SDSS_u[g]=float(table['col4'][g])
        APASS_B[g]=float(table['col5'][g])
        SDSS_g[g]=float(table['col6'][g])
        GAIA2_Gbp[g]=float(table['col7'][g])
        APASS_V[g]=float(table['col8'][g])
        SDSS_r[g]=float(table['col9'][g])
        GAIA2_G[g]=float(table['col10'][g])
        SDSS_i[g]=float(table['col11'][g])
        GAIA2_Grp[g]=float(table['col12'][g])
        twoMASS_J[g]=float(table['col13'][g])
        twoMASS_H[g]=float(table['col14'][g])
        twoMASS_Ks[g]=float(table['col15'][g])
        WISE_W1[g]=float(table['col16'][g])
        WISE_W2[g]=float(table['col17'][g])
        WISE_W3[g]=float(table['col18'][g])
        WISE_W4[g]=float(table['col19'][g])
        g+=1


    dat={'MODEL':modelname,'GALEX_FUV':GALEX_FUV, 'GALEX_NUV':GALEX_NUV, 'SDSS_u':SDSS_u, 'APASS_B':APASS_B, 'SDSS_g':SDSS_g, 'GAIA2_Gbp':GAIA2_Gbp, 'APASS_V':APASS_V, 'SDSS_r':SDSS_r, 'GAIA2_G':GAIA2_G, 'SDSS_i':SDSS_i, 'GAIA2_Grp':GAIA2_Grp, '2MASS_J':twoMASS_J, '2MASS_H':twoMASS_H, '2MASS_Ks':twoMASS_Ks, 'WISE_W1':WISE_W1, 'WISE_W2':WISE_W2, 'WISE_W3':WISE_W3, 'WISE_W4':WISE_W4}
    mod=pd.DataFrame(data=dat,columns=['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'APASS_B', 'SDSS_g', 'GAIA2_Gbp', 'APASS_V', 'SDSS_r', 'GAIA2_G', 'SDSS_i', 'GAIA2_Grp', '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4'],index=modelname)
    return mod,modelname
"""  
"""
#5000-10000
    while temp_bb<9999:
        print("BB: ",str(temp_bb))
        filename_bb="bbody_bb0"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb0"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb0"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table_bb=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table_bb)
        
        met=np.array(["p00","p02","p05","m05","m10","m15","m20","m25"])
        logg=np.array(["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"])
        temp=3500 #goes to 6500 in incr of 250
        while temp<4600:
            a=0
            for a in np.arange(len(met)):
                b=0
                for b in np.arange(len(logg)):
                    filename_kur="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                    if os.path.exists(filename_kur)==False: continue
                    e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                    header1 = e.readline()
                    header2 = e.readline()
                    header3 = e.readline()
                    header4 = e.readline()
                    header5 = e.readline()
                    header6 = e.readline()
                    header7 = e.readline()
                    header8 = e.readline()
                    header9 = e.readline()
                    npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                    wavelength_kur=np.zeros(npts_ms)
                    flux_kur=np.zeros(npts_ms)
                    band_kur=([])
            
                    m=0
                    for line in e:
                        line = line.strip()
                        columns = line.split()
                        if len(columns)<3:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=0.0
                            flux_kur[m]=0.0
                        else:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=np.float(columns[1])
                            flux_kur[m]=np.float(columns[2])
                        m+=1
                    
                    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                    
                    model_kur=spec_model(table)
                    
                    model_both=model_kur+model_bb
                    scale=spec[1][4]/model_both[1][4]
                    smodel=model_both[1]*scale
                    chi2total=0.
                    for c in np.arange(0,len(smodel)):
                        chi2=(abs(smodel[c]-spec[1][c]))**2
                        chi2total+=chi2
            #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        temp_mid=temp_bb 
            print("KUR: ",str(temp))
            print(chi2comp)
            temp+=250    
        temp_bb+=1000
    print("5000 cleared")

#10000-20000
    temp_bb=10000
    while temp_bb<20020:
        print("BB:",str(temp_bb))
        filename_bb="bbody_bb"+str(temp_bb)+".phot.dat"
        if os.path.exists(filename_bb)==False: continue
        g = open("bbody_bb"+str(temp_bb)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp_bb)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table)
            
        met=np.array(["p00","p02","p05","m05","m10","m15","m20","m25"])
        logg=np.array(["0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"])
        temp=3500 #goes to 6500 in incr of 250
        while temp<6600:
            a=0
            for a in np.arange(len(met)):
                b=0
                for b in np.arange(len(logg)):
                    filename_kur="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                    if os.path.exists(filename_kur)==False: continue
                    e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                    header1 = e.readline()
                    header2 = e.readline()
                    header3 = e.readline()
                    header4 = e.readline()
                    header5 = e.readline()
                    header6 = e.readline()
                    header7 = e.readline()
                    header8 = e.readline()
                    header9 = e.readline()
                    npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                    wavelength_kur=np.zeros(npts_ms)
                    flux_kur=np.zeros(npts_ms)
                    band_kur=([])
            
                    m=0
                    for line in e:
                        line = line.strip()
                        columns = line.split()
                        if len(columns)<3:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=0.0
                            flux_kur[m]=0.0
                        else:
                            band_kur=np.append(band_kur,columns[0])
                            wavelength_kur[m]=np.float(columns[1])
                            flux_kur[m]=np.float(columns[2])
                        m+=1
                    
                    dat={'wavelength':wavelength_kur,'flux':flux_kur}
                    table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                    
                    model_kur=spec_model(table)
                    
                    model_both=model_kur+model_bb
                    scale=spec[1][4]/model_both[1][4]
                    smodel=model_both[1]*scale
                    chi2total=0.
                    for c in np.arange(0,len(smodel)):
                        chi2=(abs(smodel[c]-spec[1][c]))**2
                        chi2total+=chi2
            #                print(chi2total)
                    if chi2total<chi2comp:
                        chi2comp=chi2total
                        temp_mid=temp_bb
            print("KUR: ",str(temp))
            print(chi2comp)
            temp+=250    
        temp_bb+=1000
    print("10000 cleared")

    #20000-200000 K
    temp_bb=20000
"""
    
#rest of temps
"""    
    temp=5000
    while temp<9999:
        filename="bbody_bb0"+str(temp)+".phot.dat"
        if os.path.exists(filename)==False: continue
        g = open("bbody_bb0"+str(temp)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb0"+str(temp)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table_bb=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table_bb)
            
        scale=spec[1][1]/model_bb[1][1]
        smodel_bb=model_bb[1]*scale
        chi2total=0.
        for c in np.arange(0,2):
            chi2=(abs(smodel_bb[c]-spec[1][c]))**2
            chi2total+=chi2
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=filename
            scale_best=scale
            flux_best=model_bb         
        temp+=50
    temp=10000
    while temp<20020:
        filename="bbody_bb"+str(temp)+".phot.dat"
        if os.path.exists(filename)==False: continue
        g = open("bbody_bb"+str(temp)+".phot.dat")
        header1_2 = g.readline()
        header2_2 = g.readline()
        header3_2 = g.readline()
        header4_2 = g.readline()
        header5_2 = g.readline()
    
        npts_wd = sum(1 for line in open("bbody_bb"+str(temp)+".phot.dat"))-5
        band_bb=([])
        flux_bb=np.zeros(npts_wd)
        wavelength_bb=np.zeros(npts_wd)
        
        q=0
        for line in g:
            line = line.strip()
            columns = line.split()
            if len(columns)<3:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=0.0
                flux_bb[q]=0.0
            else:
                band_bb=np.append(band_bb,columns[0])
                wavelength_bb[q]=columns[1]
                flux_bb[q]=columns[2]
            q+=1
            
        dat={'wavelength':wavelength_bb,'flux':flux_bb}
        table=pd.DataFrame(data=dat,index=band_bb,columns=['wavelength','flux'])
            
        model_bb=spec_model(table)
            
        scale=spec[1][1]/model_bb[1][1]
        smodel_bb=model_bb[1]*scale
        chi2total=0.
        for c in np.arange(0,2):
            chi2=(abs(smodel_bb[c]-spec[1][c]))**2
            chi2total+=chi2
        if chi2total<chi2comp:
            chi2comp=chi2total
            model_best=filename
            scale_best=scale
            flux_best=model_bb         
        temp+=50
    temp=20000

def fits_kur(spec):
    chi2comp=2e62
    met=np.array(["m05","m10"])
    logg=np.array(["0.5","1.0"])
    temp=3500 #goes to 6500 in incr of 250
    #took out 0.0 log(g), is there a way to keep it/skip?
    while temp<4600:
        a=0
        for a in np.arange(len(met)):
            b=0
            for b in np.arange(len(logg)):
                filename="Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"
                if os.path.exists(filename)==False: continue
                e=open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat")
                header1 = e.readline()
                header2 = e.readline()
                header3 = e.readline()
                header4 = e.readline()
                header5 = e.readline()
                header6 = e.readline()
                header7 = e.readline()
                header8 = e.readline()
                header9 = e.readline()
                npts_ms = sum(1 for line in open("Kurucz_f"+met[a]+"k2odfnew.pck.teff="+str(temp)+"..logg="+logg[b]+"0000.phot.dat"))-9
                wavelength_kur=np.zeros(npts_ms)
                flux_kur=np.zeros(npts_ms)
                band_kur=([])
        
                m=0
                for line in e:
                    line = line.strip()
                    columns = line.split()
                    if len(columns)<3:
                        band_kur=np.append(band_kur,columns[0])
                        wavelength_kur[m]=0.0
                        flux_kur[m]=0.0
                    else:
                        band_kur=np.append(band_kur,columns[0])
                        wavelength_kur[m]=np.float(columns[1])
                        flux_kur[m]=np.float(columns[2])
                    m+=1
                
                dat={'wavelength':wavelength_kur,'flux':flux_kur}
                table=pd.DataFrame(data=dat,index=band_kur,columns=['wavelength','flux'])
                
                minfo=np.zeros((2,len(flux_kur)))
                minfo[0]=wavelength_kur
                minfo[1]=flux_kur
                
                model_ms=spec_model(table)

                wpos=wavelength_kur[flux_kur>0]
                fpos=flux_kur[flux_kur>0]
                mpos=np.zeros((2,len(wpos)))
                mpos[0]=wpos
                mpos[1]=fpos
                
                scale=spec[1][4]/model_ms[1][4]
                smodel_ms=model_ms[1]*scale
                chi2total=0.
                for c in np.arange(2,7):
                    chi2=((smodel_ms[c]-spec[1][c])/spec[2][c])**2
                    chi2total+=chi2
        #                print(chi2total)
                if chi2total<chi2comp:
                    chi2comp=chi2total
                    model_best=filename
                    scale_best=scale
                    flux_best=mpos
                    info=minfo

        temp+=250
    return chi2comp,model_best,scale_best,flux_best  
"""
    