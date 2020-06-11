#basic Python imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#imports for manipulating astronomical data
from astropy.coordinates import SkyCoord
from astropy.io import ascii

#extinction correction package
from dustmaps.sfd import SFDQuery

#curve fitting imports
from scipy.stats import gaussian_kde,norm
from scipy.interpolate import make_interp_spline, BSpline

from sklearn import neighbors
from scipy.signal import savgol_filter

#class for processing starting from the raw WFCAM datasets
class data_load:
    #reads in and performs initial cuts on data from chosen galaxy
    #change optional arguments to false to skip initial cuts
    #path_to_file argument used to specify where WFCAM data is stored
    def __init__(self,galaxy, CLS=True, mag=True, ext=True, path_to_file='initial_data/'):


        
        #convert hhmmss ra and ddmmss dec into decimal values
        
        def make_deg_coord(frame):
            
            #takes relevant columns as arguements, converts into single ra and single dec column
            
            def coordtransform(rah,ram,ras,dd,dm,ds):
            
                #transformations carried out
                
                #ra
                
                long = (rah + np.divide(ram,60) + np.divide(ras,3600))*15
                
                #dec
                
                lat = dd + np.divide(dm,60) + np.divide(ds,3600)
                
                #array containing transformed co-ordinates returned
        
                return(np.array([long,lat]))
                
            #function called    
            
            coords=coordtransform(frame.rah,frame.ram,frame.ras,frame.dd,frame.dm,frame.ds)
                
            #columns added to DataFrame
                
            frame['RA']=coords[0]
            
            frame['DEC']=coords[1]
            
        def make_tan_coord(frame):
            
                def create_tangent_coords(tangentra,tangentdec):
        
        #ra and dec attributes converted to variables in radians
        
                    ra = np.radians(self.ra)
                    dec = np.radians(self.dec)
                    
                    #tangent co-ordinates also converted to radians
                    
                    tanra = np.radians(tangentra)
                    tandec = np.radians(tangentdec)
                    
                    #conversion for xi carried out
                    
                    xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
                    
                    #conversion for eta carried out
                    
                    eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
                    
                    #co-ordinates converted to degrees and set as attributes
                    
                    self.xi = xi * (180/np.pi)
                    self.eta = eta * (180/np.pi)
                    
                    return (np.array([self.xi,self.eta]))
                
                if frame.galaxy=='ngc147':
                    tra=8.300500
                    tdec=48.50850
                elif frame.galaxy=='ngc185':
                    tra=9.7415417
                    tdec=48.3373778
                elif frame.galaxy=='ngc205':
                    tra=10.09189356
                    tdec=41.68541564
                elif frame.galaxy=='m32':
                    tra=10.6742708
                    tdec=40.8651694
                    
                coords=create_tangent_coords(tra,tdec)
                
                frame['xi']=coords[0]
                frame['eta']=coords[1]
                    
        #function to cut DataFrame based on cls index
        
        def CLS_cut(frame):
            
            #removes any non-stellar sources from data
            for i in range(len(frame.jcis)):
                if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0) or (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0): 
                    frame.loc[i]=np.nan
            
            #percentage decrease of catalogue length calculated and printed
            
            k=len(frame.jcis)
            l=len(frame.jcis.dropna())            
            decrease=100-(l/k)*100
            
            print('CLS cut reduced data by ' + str(int(decrease)) + '%')
            
            #for some reason if I put frame=frame.dropna() here it doesn't do anything, so I have to wipe the NaNs outside the function
            
            
        #cuts DataFrame with large magnitude error
        
        def mag_err_cut(frame,magerr=0.2):
            
            #have to start being careful about referencing index column of dataframe now since dropna() will have been used
            
            for i in frame.index:
                if frame.kerr[i] > magerr or frame.herr[i] > magerr or frame.kerr[i] > magerr:
                
                    frame.loc[i]=np.nan
            
            #percentage decrease of catalogue length calculated and printed
            
            k=len(frame.jcis)
            l=len(frame.jcis.dropna())
            
            decrease=100-(l/k)*100
            
            print('Magerr cut reduced data by ' + str(int(decrease)) + '%')
            
        #applies extinction correction to data in DataFrame format    
        
        def ext_corr(frame):
            

            #attributes set as variables for ease
            
            ra = frame.RA
            dec = frame.DEC
            
            #astropy skycoord called, units and frame set
            
            coords=SkyCoord(ra,dec,unit='deg',frame='icrs')
            
            #sfd map chosen
            
            sfd=SFDQuery()
            
            #E(B-V) for each star loaded into array from map
            
            sfdred=sfd(coords)
               
            #corrections from table 6 of Schlafly and Finkbeiner(2011) made
            
            jext = sfdred * 0.709
            hext = sfdred * 0.449
            kext = sfdred * 0.302
            
            #extinction corrections carried out on class attributes
            
            frame.jmag=frame.jmag - jext
            frame.hmag=frame.hmag - hext
            frame.kmag=frame.kmag - kext
            
            print('UKIRT extinction corrections done')
            



        #with initial functions defined, we now need to convert the ASCII WFCAM data into a dataframe
        
        #column headers of WFCAM data defined
        
        columnames=['rah','ram','ras','dd','dm','ds','unknown','xj','yj','jmag','jerr','jcis','xh','yh','hmag','herr','hcis','xk','yk','kmag','kerr','kcis']
        
        #galaxy list for reading in the correct file
        
        galaxies=['ngc147','ngc185','ngc205','m32']
        
        #names of datafiles
        
        self.galaxies=galaxies
        self.galaxy=galaxy
        
        infilenames=['lot_n147.unique','lot_n185.unique','N205_new_trimmed.unique','M32_new.asc']
        
        
        
        #associates galaxy input with appropriate file
        
        for i in range(len(infilenames)):
                
            if galaxy==galaxies[i]:
                    
                file=infilenames[i]
                break
        
        #ascii data read in to astropy table
        
        frame=ascii.read(path_to_file + file)
        
        #converted to dataframe
        
        frame=frame.to_pandas()
        
        #columns assigned
        
        print(frame)
        
        if self.galaxy=='ngc205' or self.galaxy=='m32':
            
            frame=frame.drop(['col26'],axis=1)
            frame=frame.drop(['col25'],axis=1)
            frame=frame.drop(['col24'],axis=1)
            frame=frame.drop(['col23'],axis=1)
        
        
        frame.columns=columnames

        
        #functions called to create decimal coordinate columns and flag column
        make_deg_coord(frame)
        make_flag_col(frame)
        
        #optional skip when initiating class for the cuts/extinction corrections
        
        if CLS==True:
            
            
            #cls cut carried out, NaN values purged
            CLS_cut(frame)
            frame=frame.dropna()
        if mag==True:
            
            #mag cut carried out, NaN values purged
            mag_err_cut(frame)
            
            frame=frame.dropna()
            
        if ext==True:
            
            #extinction corrections done
            ext_corr(frame)
        
        #data set as class attribute
        
        self.data=frame
        
    def forecut(self):
        
        foredata=self.data.copy()
        data=self.data
        
        forecuts=[0.98,0.98,0.965,0.98]
        
        galaxies=self.galaxies
        
        for i in range(len(forecuts)):
            
            if self.galaxy==galaxies[i]:
                cut=forecuts[i]
                break
        
        jk=data.jmag-data.kmag
        
        for i in data.index:
            
            if jk[i] < cut:
                
                data.loc[i]=np.nan
                
            else:
                foredata.loc[i]=np.nan
        k=len(data.copy())
        data=data.dropna()
        foredata=foredata.dropna()
        decrease=100-(len(data)/k) * 100
        
        print('Foreground cut reduced data by ' + str(int(decrease)) + '%, ' + str(len(data)) + ' sources retained')
        
        self.foredata=foredata
        
    def trgbcut(self):
        
        data=self.data
        rgbdata=self.data.copy()
        trgbcuts=[18.13,17.84,17.94,17.9]
        
        galaxies=self.galaxies
        
        for i in range(len(galaxies)):
            
            if self.galaxy==galaxies[i]:
                
                cut=trgbcuts[i]
                break
                
        
        for i in data.index:
            
            if data.kmag[i] > cut:
                
                data.loc[i]=np.nan
                
            else:
                
                rgbdata.loc[i]=np.nan
                
                
                
        k=len(data.copy())
        data=data.dropna()
        rgbdata=rgbdata.dropna()
        
        decrease=100-(len(data)/k) * 100
        
        print('Foreground cut reduced data by ' + str(int(decrease)) + '%, ' + str(len(data)) + ' sources retained')
        
        self.rgbdata=rgbdata
    
    def plot_kj_cmd(self,marker='o',markersize=1,color='blue'):
        
        data=self.data
        


        plt.figure()        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.kmag,data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
    
    def plot_lum(self):
        
        data=self.data
        
        lum_func=data.kmag
        
        plt.figure()
        

        sns.kdeplot(lum_func.dropna(),color='white',gridsize=100)
        plt.xlabel('$K_0$')
        line = plt.gca().get_lines()[0]
        xd = line.get_xdata()
        yd = line.get_ydata()
        
        diffy = np.diff(yd) / np.diff(xd)
        diffx = (np.array(xd)[:-1] + np.array(xd)[1:]) / 2
        
        #plt.plot(diffx,diffy,label='diff_lum_func')
        
        twodiffy = np.diff(diffy) / np.diff(diffx)
        twodiffx = (np.array(diffx)[:-1] + np.array(diffx)[1:]) / 2        
        
        #plt.plot(twodiffx,twodiffy)
        
        # 300 represents number of points to make between T.min and T.max
        xnew = np.linspace(twodiffx.min(), twodiffx.max(), 300) 
        
        spl = make_interp_spline(twodiffx, twodiffy, k=3)  # type: BSpline
        ynew = spl(xnew)
        
        plt.plot(xnew, ynew,label='diff2_lum_func')
        plt.legend()
        plt.show()
    
    def trgbfind(self, magname = 'Kmag', dmagname = 'eKmag', niter = 1000, kernel = 'epanechnikov'):
    # Set the data to find the TRGB
        mk = self.data.kmag.dropna()
        dmk = self.data.kerr.dropna()
            
        #Initialise stuff
        niter = niter           # Number of itterations 
        rtol = 1e-5             # Relative tolerance of the result
        kernel = 'epanechnikov' # Parabolic kernel for the KDE
        
        mx = np.linspace(max(mk)*1.2, min(mk)*0.8, 1000)
        trgbloc = np.zeros(niter)
        
        #----------------------------------------
        #Generate NITER realisations of the Kernel Density Estimation
        for i in range(niter):
            msamp = np.random.normal(mk, dmk)     # Add Noise to data -> diff. each loop -> more reliable TRGB
            
            # Find an ideal binwidth for the luminosity function  
            # PS: Monte Carlo already smooths the distribution, so reduce the ideal binwidth a bit.
            
            bandwidth_factor = 0.25
            bandwidth = bandwidth_factor*(np.std(msamp)*(len(msamp)**(-0.2)))
                
            #----------------------------------------
            # Implement the Kernel density estimation using a KD Tree for efficient queries 
            
            kde = neighbors.KernelDensity(bandwidth = bandwidth, rtol = rtol, kernel = kernel)  # Inialise
            kde.fit(msamp[:, np.newaxis])        # Fit the Kernel Density model on the data. 
            #kde.score_samples returns ln(pdf)   # Evaluate the density model on data - probablility density function
            pdf = np.exp(kde.score_samples(mx[:, np.newaxis]))  #MX is x-axis range which the PDF is computed/plotted.
            
            #----------------------------------------
            # Set the Edge Detection part using a savgol_filter
            
            smooth_window = 31
            poly_degree = 3
            dpdf = savgol_filter(pdf, smooth_window, poly_degree, deriv = 1)
            trgbloc[i] = mx[np.argmin(dpdf)]     # Most negative value corresponds to highest rate of decrease 
            
        trgbloc_mean = np.mean(trgbloc)  # Find the TRGB
        trgbloc_sd = np.std(trgbloc)     # Find the Error in the TRGB estimate
        
        return trgbloc_mean, trgbloc_sd
            
            
        #method to save data as binary parquet file, to be used by data_read class    
    def save_to_parquet(self,fileloc):
            
        self.data.to_parquet(fileloc)
        
        
        
        
        
        
        
        