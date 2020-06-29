from sat_graphers import sat
import numpy as np
import pandas as pd
class ratio_utils(sat):
    
    def __init__(self,*args,**kwargs):
        super(ratio_utils,self).__init__(*args,**kwargs)
        
        def CM_to_FEH(CM):
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)            
        jhc=np.array([0.74,0.70,0.77,10])
        hkc=np.array([0.44,0.48,0.60,8])
        jkc=np.array([1.33,1.30,1.42,1.50])
        jk4m=np.array([0.99,0.988,0.96,0.86])
        trgb=np.array([18,17.8,17.8,18])
        
        cuts=np.array([jhc,hkc,jkc,jk4m,trgb])
        cut_frame=[]
        for i in range(len(cuts)):

            for j in range(len(cuts[i])):
                if self.infilenames[j]==self.galaxy:
                    cut_frame.append(cuts[i][j])
                    break
        cut_frame=np.array(cut_frame)
        
        self.jhc=cut_frame[0]
        self.hkc=cut_frame[1]
        self.jkc=cut_frame[2]
        self.jk4m=cut_frame[3]
        self.trgb=cut_frame[4]
        self.CM_to_FEH=CM_to_FEH

        data=self.data.dropna()
        hk=data.hmag-data.kmag
        jh=data.jmag-data.hmag
        jk=data.jmag-data.kmag
        mdata=data.copy()
        cdata=data.copy()


        for i in data.index:
            if data.kmag[i] > self.trgb or jk[i] < self.jk4m:
                mdata.loc[i]=np.nan
                cdata.loc[i]=np.nan
            elif (hk[i] > self.hkc and jh[i] > self.jhc)  or jk[i] > self.jkc:
                mdata.loc[i]=np.nan
                
            else:
                cdata.loc[i]=np.nan
        
        self.cdata=cdata.dropna()
        self.mdata=mdata.dropna()
    
    def define_marginals(self,jkmult=1,hkmult=1,jhmult=1):
        
        for j in range(len(self.infilenames)):
            
            if self.galaxy==self.infilenames[j]:
                k=j
                break
            
        
        jkdata=self.data.copy()
        hkdata=self.data.copy()

        jherrs=[0.1,0.15,0.13,0]
        hkerrs=[0.1,0.1,0.1,0]
        jkerrs=[0.1,0.1,0.1,0.1]
        
        hk=self.data.hmag-self.data.kmag
        jh=self.data.jmag-self.data.hmag
        jk=self.data.jmag-self.data.kmag
        
        
        
        
        jkmmdata=self.data.copy()
        jkMdata=self.data.copy()
        jkcmdata=self.data.copy()
        jkCdata=self.data.copy()
        
        hkmmdata=self.data.copy()
        hkMdata=self.data.copy()
        hkcmdata=self.data.copy()
        hkCdata=self.data.copy()
        
        
        
        for i in jkdata.index:
            if jkdata.kmag[i] > self.trgb or jk[i] < self.jk4m:
                jkMdata.loc[i]=np.nan
                jkmmdata.loc[i]=np.nan
                jkcmdata.loc[i]=np.nan
                jkCdata.loc[i]=np.nan
                jkdata.loc[i]=np.nan
                
                hkMdata.loc[i]=np.nan
                hkmmdata.loc[i]=np.nan
                hkcmdata.loc[i]=np.nan
                hkCdata.loc[i]=np.nan
                hkdata.loc[i]=np.nan

                

        
        jkMdata=jkMdata.dropna()
        jkmmdata=jkmmdata.dropna()
        jkcmdata=jkcmdata.dropna()
        jkCdata=jkCdata.dropna()
        jkdata=jkdata.dropna()
        
        hkMdata=hkMdata.dropna()
        hkmmdata=hkmmdata.dropna()
        hkcmdata=hkcmdata.dropna()
        hkCdata=hkCdata.dropna()
        hkdata=hkdata.dropna()
        for i in jkdata.index:
            if jk[i]>self.jkc-jkerrs[k]:
                jkMdata.loc[i]=np.nan
        
        for i in jkdata.index:
            if jk[i]<self.jkc - jkerrs[k] or jk[i] > self.jkc:
                jkmmdata.loc[i]=np.nan

        
        for i in jkdata.index:
            if jk[i] < self.jkc or jk[i] > self.jkc + jkerrs[k]:
                jkcmdata.loc[i]=np.nan

        
        for i in jkdata.index:
            if jk[i] < self.jkc + jkerrs[k]:
                jkCdata.loc[i]=np.nan

        
        
        for i in hkdata.index:
            
            if hk[i] > self.hkc-hkerrs[k] and jh[i] > self.jhc-jherrs[k]:
                hkMdata.loc[i]=np.nan

        
        for i in hkdata.index:
            
            if (hk[i] > self.hkc and jh[i] > self.jhc) or (hk[i] < self.hkc-hkerrs[k] or jh[i] < self.jhc-jherrs[k]):
                hkmmdata.loc[i]=np.nan
                

        
        for i in hkdata.index:
            
            if (hk[i] > self.hkc+ hkerrs[k] and jh[i] > self.jhc + jherrs[k]) or (hk[i]<self.hkc or jh[i] < self.jhc):
                hkcmdata.loc[i]=np.nan

        
        for i in hkdata.index:
            
            if hk[i] < self.hkc+hkerrs[k] or jh[i] < self.jhc+jherrs[k]:
                hkCdata.loc[i]=np.nan

        
        self.jkMdata=jkMdata
        self.jkmmdata=jkmmdata
        self.jkcmdata=jkcmdata
        self.jkCdata=jkCdata
        
        self.hkMdata=hkMdata
        self.hkmmdata=hkmmdata
        self.hkcmdata=hkcmdata
        self.hkCdata=hkCdata
        
    def CM_marginals(self):
        
        for j in range(3):
            
            
        
            jkMdata=self.jkMdata
            jkmmdata=self.jkmmdata
            jkcmdata=self.jkcmdata
            jkCdata=self.jkCdata
            
            hkMdata=self.hkMdata
            hkmmdata=self.hkmmdata
            hkcmdata=self.hkcmdata
            hkCdata=self.hkCdata
            
    
            final_cdata=[]
            final_mdata=[]
        
    
    #1.5/0.7 for hk, 1/0.5 for jk
    
    #mm/
        
            for i in self.hkMdata.index:
                #if hkMdata.kmag[i]==jkMdata.kmag[i] or hkMdata.kmag[i]==jkmmdata.kmag[i]:
                
                if np.isnan(self.hkMdata.kmag[i])==False:
                
                    final_mdata.append(hkMdata.loc[i])
                
            for i in self.hkCdata.index:
                
                if np.isnan(self.hkCdata.kmag[i])==False:
                
                    final_cdata.append(hkCdata.loc[i])
                
            for i in self.hkmmdata.index:
                if hkmmdata.kmag[i]==jkMdata.kmag[i] or hkmmdata.kmag[i]==jkcmdata.kmag[i] or hkmmdata.kmag[i]==jkmmdata.kmag[i]:
                    
                    if j!=2 or hkmmdata.kmag[i]==jkMdata.kmag[i]:
                    
                        final_mdata.append(hkmmdata.loc[i])
                        
                    else:
                        
                        final_cdata.append(hkmmdata.loc[i])
                    
                elif hkmmdata.kmag[i]==jkCdata.kmag[i]:
                    
                        final_cdata.append(hkmmdata.loc[i])

                    
            for i in self.hkcmdata.index:
                if hkcmdata.kmag[i]==jkCdata.kmag[i] or hkcmdata.kmag[i]==jkmmdata.kmag[i] or hkcmdata.kmag[i]==jkcmdata.kmag[i]:
                    
                    if j!=1 or hkcmdata.kmag[i]==jkCdata.kmag[i]:
                    
                        final_cdata.append(hkcmdata.loc[i])
                        
                    else:
                        
                        final_mdata.append(hkcmdata.loc[i])
                    
                elif hkcmdata.kmag[i]==jkMdata.kmag[i]:
                    
                        final_mdata.append(hkcmdata.loc[i])
                        


            
            final_cdata=pd.concat(final_cdata)
            final_mdata=pd.concat(final_mdata)

            
            final_cdata=final_cdata.dropna()
            final_mdata=final_mdata.dropna()
            
            cpois=np.sqrt(len(final_cdata.kmag))
            mpois=np.sqrt(len(final_mdata.kmag))
            
            mnum=len(final_cdata.kmag)
            cnum=len(final_mdata.kmag)
            
            cm=np.divide(len(final_cdata.kmag),len(final_mdata.kmag))
            
            if j==0:
                cm_final=cm
                FEH= self.CM_to_FEH(cm)
                cmpois= cm * np.sqrt((mpois/mnum)**2 + (cpois/cnum)**2)
                
            elif j==1:
                Mcm=cm
                
            elif j==2:
                Ccm=cm
                
        upsystunc= np.abs(cm_final-Mcm)
        downsystunc=np.abs(cm_final-Ccm)
        
        finupcm=np.sqrt(cmpois**2 + upsystunc**2)
        findowncm=np.sqrt(cmpois**2 + downsystunc**2)
        
        print(finupcm)
        print(findowncm)
        
        fehup=self.CM_to_FEH(cm_final + finupcm)
        fehdown=self.CM_to_FEH(cm_final-findowncm)
            
        print('Average C/M ratio = ' + str(cm_final) +', average [Fe/H] = ' + str(FEH) + '+/- ' + str(np.abs(fehup-FEH)) + '/' + str(np.abs(fehdown-FEH)) )
                

                
                
                
        
                
                



            
        
        
    
    def total_CM_FEH(self):
        
        cm=np.divide(len(self.cdata),len(self.mdata))
        FEH= self.CM_to_FEH(cm)
        
        print('Average C/M ratio = ' + str(cm) +', average [Fe/H] = ' + str(FEH))
        
        