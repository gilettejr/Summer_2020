from data_read import data_read
from astropy.io import ascii
import numpy as np

import matplotlib.pyplot as plt
class isochrones(data_read):
    
    def plot_isos(self,isofile,graph='kj_cmd',overlay=True):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        def merr(dist,derr):
            
            m=5*np.log10(dist/10)-5*np.log10((dist-derr)/10)
            
            return m
        def mag(a,b):
            c=np.sqrt(a**2+b**2)
            return c
        
        self.merr=merr
        self.mag=mag
        
        dists=[(680000,30000),(630000,30000),(820000,30000),(785000,25000),(0,0)]
        
        for i in range(len(self.galaxies)):
            
            if self.galaxy==self.galaxies[i]:
                
                distance=dists[i][0]
                disterr=dists[i][1]
                
                break
        
        names=['z','MH','age','Mini','int_IMF','mass','logL','logTe','logg','label','McoreTP','C_O','period0','period1','pmode','Mloss','tau1m','x','y','xc','xn','xo','cexcess','z_','bolmag','zmag','ymag','jmag','hmag','kmag']
        
        iso=ascii.read(isofile,names=names)
        
        iso=iso.to_pandas()
        

        
        iso.jmag=apparent(iso.jmag,distance)
        iso.hmag=apparent(iso.hmag,distance)
        iso.kmag=apparent(iso.kmag,distance)
        
        print('Is this the long bit?')
        
        for i in iso.index:
            
            #label=8 is the AGB phase in padova isochrones
            
            if iso.label[i]!=7 and iso.label[i]!=8:

                iso.loc[i]=np.nan
        
        #same method for seperating out the different isochrones in the set
        #as in self.isoplot
        
        iso=iso.dropna()
        iso=iso.reset_index(drop=True)
        print(iso)
        
        if graph=='kj_cmd':
            
            xdata=iso.jmag-iso.kmag
            ydata=iso.kmag
            
        elif graph=='cc':
            
            xdata=iso.jmag-iso.hmag
            ydata=iso.hmag-iso.kmag
            
        else: 
            
            print('Invalid graph input in isochrone method')
        
        indices=[]
        
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 

        
        plt.plot(xdata[:indices[0]],ydata[:indices[0]],label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
        
        for i in range(len(indices)):
            if i==(len(indices)-1):
                plt.plot(xdata[indices[i]:],ydata[indices[i]:],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                break
            else:
                plt.plot(xdata[indices[i]:indices[i+1]],ydata[indices[i]:indices[i+1]],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        
        if overlay==False:
            
            if graph=='kj_cmd':
            
                plt.gca().invert_yaxis()
                plt.ylabel('$K_0$')
                plt.xlabel('$J_0$-$K_0$')

                
            elif graph=='cc':
                
                plt.ylabel('$H_0$-$K_0$')
                plt.xlabel('$J_0$-$H_0$')
        
        plt.legend()
        
    def plot_gaisos(self,isofile,overlay=False):
        
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        def merr(dist,derr):
            
            m=5*np.log10(dist/10)-5*np.log10((dist-derr)/10)
            
            return m
        def mag(a,b):
            c=np.sqrt(a**2+b**2)
            return c
        
        self.merr=merr
        self.mag=mag
        
        dists=[(680000,30000),(630000,30000),(820000,30000),(785000,25000),(0,0)]
        
        for i in range(len(self.galaxies)):
            
            if self.galaxy==self.galaxies[i]:
                
                distance=dists[i][0]
                disterr=dists[i][1]
                
                break
        
        names=['z','MH','age','Mini','int_IMF','mass','logL','logTe','logg','label','McoreTP','C_O','period0','period1','pmode','Mloss','tau1m','x','y','xc','xn','xo','cexcess','z_','bolmag','Gmag','bpprmag','bpftmag','rpmag']
        
        iso=ascii.read(isofile,names=names)
        
        iso=iso.to_pandas()
        

        
        iso.Gmag=apparent(iso.Gmag,distance)
        
        for i in iso.index:
            
            #label=8 is the AGB phase in padova isochrones
            
            if iso.label[i]!=7 and iso.label[i]!=8:

                iso.loc[i]=np.nan
        
        #same method for seperating out the different isochrones in the set
        #as in self.isoplot
        
        iso=iso.dropna()
        iso=iso.reset_index(drop=True)


            
        xdata=iso.Gmag-iso.rpmag
        ydata=iso.Gmag
        
        indices=[]
        
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 

        
        plt.plot(xdata[:indices[0]],ydata[:indices[0]],label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
        
        for i in range(len(indices)):
            if i==(len(indices)-1):
                plt.plot(xdata[indices[i]:],ydata[indices[i]:],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                break
            else:
                plt.plot(xdata[indices[i]:indices[i+1]],ydata[indices[i]:indices[i+1]],label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        
        if overlay==False:
                       
            plt.gca().invert_yaxis()
            plt.ylabel('$G_0$')
            plt.xlabel('$G_0')

        plt.legend()
        
    def plot_simlum(self,isofile,linestyle='solid',overlay=False,pop='multiple'):
                
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        def merr(dist,derr):
            
            m=5*np.log10(dist/10)-5*np.log10((dist-derr)/10)
            
            return m
        def mag(a,b):
            c=np.sqrt(a**2+b**2)
            return c
        
        self.merr=merr
        self.mag=mag
        
        dists=[(680000,30000),(630000,30000),(820000,30000),(785000,25000),(0,0)]
        
        for i in range(len(self.galaxies)):
            
            if self.galaxy==self.galaxies[i]:
                
                distance=dists[i][0]
                disterr=dists[i][1]
                
                break
        
        names=['age','z','magbin','mbolmag','zmag','ymag','jmag','hmag','kmag']
        
        iso=ascii.read(isofile,names=names)
        
        iso=iso.to_pandas()
        
        print(iso.kmag)
        
        iso.magbin=apparent(iso.magbin,distance)
        
        print(iso.magbin)
            #label=8 is the AGB phase in padova isochrones
            

        #same method for seperating out the different isochrones in the set
        #as in self.isoplot
        

        
        xdata=iso.magbin
        ydata=iso.kmag
        
        indices=[]
        
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 
        if pop!='single':
        
        
            plt.plot(xdata[:indices[0]],ydata[:indices[0]],linestyle=linestyle,label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
            print(xdata[:indices[0]])
            print(ydata[:indices[0]])
            for i in range(len(indices)):
                if i==(len(indices)-1):
                    plt.plot(xdata[indices[i]:],ydata[indices[i]:],linestyle=linestyle,label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                    break
                else:
                    plt.plot(xdata[indices[i]:indices[i+1]],ydata[indices[i]:indices[i+1]],linestyle=linestyle,label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        else:
            
            plt.plot(xdata,ydata,linestyle=linestyle,marker='o',markersize='1',label='log(t) = ' + str(iso.age[0]) + ', Z = ' + str(iso.z[0]))
        
        if overlay==False:
            

                plt.ylabel('N_stars')
                plt.xlabel('$K_0$')
                
                plt.yscale('log')
                
        
        plt.legend()            
        
    def plot_simpop(self,isofile,graph='kj_cmd',overlay=False,pop='single'):
                
        def apparent(M,dguess):
            m = M + 5*np.log10(dguess/10)
            return m
        
        def merr(dist,derr):
            
            m=5*np.log10(dist/10)-5*np.log10((dist-derr)/10)
            
            return m
        def mag(a,b):
            c=np.sqrt(a**2+b**2)
            return c
        
        self.merr=merr
        self.mag=mag
        
        dists=[(680000,30000),(630000,30000),(820000,30000),(785000,25000),(0,0)]
        
        for i in range(len(self.galaxies)):
            
            if self.galaxy==self.galaxies[i]:
                
                distance=dists[i][0]
                disterr=dists[i][1]
                
                break
        
        names=['z','age','Mini','mass','logL','logTe','logg','label','McoreTP','C_O','period0','period1','pmode','Mloss','tau1m','x','y','xc','xn','xo','cexcess','z_','bolmag','zmag','ymag','jmag','hmag','kmag']
        
        iso=ascii.read(isofile,names=names)
        
        iso=iso.to_pandas()
        

        
        iso.jmag=apparent(iso.jmag,distance)
        iso.hmag=apparent(iso.hmag,distance)
        iso.kmag=apparent(iso.kmag,distance)
        
            
            #label=8 is the AGB phase in padova isochrones
            

        #same method for seperating out the different isochrones in the set
        #as in self.isoplot
        
        iso=iso.dropna()
        iso=iso.reset_index(drop=True)

        
        if graph=='kj_cmd':
            
            xdata=iso.jmag-iso.kmag
            ydata=iso.kmag
            
        elif graph=='cc':
            
            xdata=iso.jmag-iso.hmag
            ydata=iso.hmag-iso.kmag
            
        else: 
            
            print('Invalid graph input in isochrone method')
        
        indices=[]
        
        for i in range(1,len(iso.age)):
            if iso.age[i]!=iso.age[i-1] or iso.z[i]!=iso.z[i-1]:
                indices.append(i)
 
        if pop!='single':
        
        
            plt.plot(xdata[:indices[0]],ydata[:indices[0]],linestyle='none',label='log(t) = ' + str(iso.age[0]) + ', Z = ' +str(iso.z[0]))
            
            for i in range(len(indices)):
                if i==(len(indices)-1):
                    plt.plot(xdata[indices[i]:],ydata[indices[i]:],linestyle='none',label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
                    break
                else:
                    plt.plot(xdata[indices[i]:indices[i+1]],ydata[indices[i]:indices[i+1]],linestyle='none',label='log(t) = ' + str(iso.age[indices[i]]) + ', Z = ' +str(iso.z[indices[i]]))
        
        else:
            
            plt.plot(xdata,ydata,linestyle='none',marker='o',markersize='1',label='log(t) = ' + str(iso.age[0]) + ', Z = ' + str(iso.z[0]))
        
        if overlay==False:
            
            if graph=='kj_cmd':
            
                plt.gca().invert_yaxis()
                plt.ylabel('$K_0$')
                plt.xlabel('$J_0$-$K_0$')

                
            elif graph=='cc':
                
                plt.ylabel('$H_0$-$K_0$')
                plt.xlabel('$J_0$-$H_0$')
        
        plt.legend()                
                
        