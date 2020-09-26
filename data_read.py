from data_load import data_loader

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#libraries for finding trgb
from sklearn import neighbors
from scipy.signal import savgol_filter

#libraries for seperating C and M with non horizontal/perpendicular lines




#class inherits all methods from data_loads, but takes data from defined binary file rather than loading from raw data

class data_reader(data_loader):
    #load in data from appropriate file depending on stage argument and galaxy
    def __init__(self,stage='cm',galaxy='ngc147'):
        
        def data_corners(data):
            
            corner1=(np.max(data.xi),np.max(data.eta))
            corner2=(np.max(data.xi),np.min(data.eta))
            corner3=(np.min(data.xi),np.min(data.eta))
            corner4=(np.min(data.xi),np.max(data.eta))
            
            return np.array([corner1,corner2,corner3,corner4])
        
        self.data_corners=data_corners
        
        #function to cut data to the left/right of non vertical line
        #pretty much redundant
        
        def eq_to_tan(ra,dec,tangentra,tangentdec):
            

        
            ra = np.radians(ra)
            dec = np.radians(dec)
            
            #tangent co-ordinates also converted to radians
            
            tanra = np.radians(tangentra)
            tandec = np.radians(tangentdec)
            
            #conversion for xi carried out
            
            xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
            
            #conversion for eta carried out
            
            eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
            
            #co-ordinates converted to degrees and set as attributes

            return xi,eta
        
        
        self.eq_to_tan=eq_to_tan
        
        def rotate_coords(xi,eta,theta):
            
            xi_rad=np.radians(xi.copy())
            eta_rad=np.radians(eta.copy())
            
            alpha=eta_rad*np.sin(theta)+xi_rad * np.cos(theta)
            beta=eta_rad * np.cos(theta)+ xi_rad * -np.sin(theta)
            
            alpha=np.degrees(alpha)
            beta=np.degrees(beta)
            
            return np.array([alpha,beta])
        
        self.rotate_coords=rotate_coords
        
        #set of if/elif/else statements to choose correct file
        
        if stage=='agb':
            path_to_folder='processed_data/agb_data/'
            
        elif stage=='agb_crossed':
            path_to_folder='processed_data/agb_crossed_data/'
            
        elif stage=='cls_crossed':
            path_to_folder='processed_data/cls_crossed_data/'
            
        elif stage=='cls_cut':
            path_to_folder='processed_data/cls_cut_data/'
        
        elif stage=='fore_cut':
            path_to_folder='processed_data/fore_cut_data/'
            
        elif stage=='cm':
            mpath_to_folder='processed_data/m_agb_data/'
            cpath_to_folder='processed_data/c_agb_data/'
            path_to_folder='processed_data/agb_crossed_data/'
        else:
            path_to_folder=stage
            
            
        #slightly different method needed for c and m stars
        #since these are in two separate files
        if stage=='cm':
            mframe=pd.read_parquet(mpath_to_folder + galaxy)
            cframe=pd.read_parquet(cpath_to_folder + galaxy)
            frame=pd.read_parquet(path_to_folder + galaxy)
            
            self.mdata=mframe
            self.cdata=cframe
            self.data=frame
            
        #otherwise, whole dataset set to class atribute
        else:
            
            frame=pd.read_parquet(path_to_folder + galaxy)
            
            self.data=frame


            
        self.galaxy=galaxy
        
        #galaxies list created for future matching
        
        galaxies=['ngc147','ngc185','ngc205','m32','m31','and1','and2','and3','and6','and7','and10','and14','and15','and16','and17','and18','and19','and20']
        
        self.galaxies=galaxies
        
        #class function defined for C/m tp [Fe/H] conversion
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
            
        self.CM_to_FEH=CM_to_FEH
        
        def CM_to_FEH_uncs(CM,deltaCM):
            
            deltalogCM=deltaCM/(CM*np.log(10))
            
            logCM=np.log10(CM)
            
            rhs=0.47 * logCM
            
            deltarhs= rhs * np.sqrt((deltalogCM/logCM)**2 + (0.1/0.47)**2)
            

            
            deltaFEH=np.sqrt(0.06**2+deltarhs**2)
            
            return deltaFEH
            
        self.CM_to_FEH_uncs=CM_to_FEH_uncs
    #produces trgb location and associated error
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
            
            bandwidth_factor = 0.1
            bandwidth = bandwidth_factor*(np.std(msamp)*(len(msamp)**(-0.2)))
                
            #----------------------------------------
            # Implement the Kernel density estimation using a KD Tree for efficient queries 
            
            kde = neighbors.KernelDensity(bandwidth = bandwidth, rtol = rtol, kernel = kernel)  # Inialise
            kde.fit(msamp[:, np.newaxis])        # Fit the Kernel Density model on the data. 
            #kde.score_samples #returns ln(pdf)   # Evaluate the density model on data - probablility density function
            pdf = np.exp(kde.score_samples(mx[:, np.newaxis]))  #MX is x-axis range which the PDF is computed/plotted.
            
            plt.plot(mx,pdf)
            
            #----------------------------------------
            # Set the Edge Detection part using a savgol_filter
            
            smooth_window = 31
            poly_degree = 3
            dpdf = savgol_filter(pdf, smooth_window, poly_degree, deriv = 1)
            trgbloc[i] = mx[np.argmin(dpdf)]     # Most negative value corresponds to highest rate of decrease 
            
        trgbloc_mean = np.mean(trgbloc)  # Find the TRGB
        trgbloc_sd = np.std(trgbloc)     # Find the Error in the TRGB estimate
        
        return [trgbloc_mean, trgbloc_sd]
    
    def hkcmfind(self, magname = 'Kmag', dmagname = 'eKmag', niter = 1000, kernel = 'epanechnikov'):
    # Set the data to find the TRGB
        mk = self.data.hmag.dropna()-self.data.kmag.dropna()
        dmk = np.sqrt(self.data.kerr.dropna()**2 + self.data.herr.dropna()**2)
        
        mk=-mk
        
            
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
            #kde.score_samples #returns ln(pdf)   # Evaluate the density model on data - probablility density function
            pdf = np.exp(kde.score_samples(mx[:, np.newaxis]))  #MX is x-axis range which the PDF is computed/plotted.
            
            plt.plot(mx,pdf)
            
            #----------------------------------------
            # Set the Edge Detection part using a savgol_filter
            
            smooth_window = 31
            poly_degree = 3
            dpdf = savgol_filter(pdf, smooth_window, poly_degree, deriv = 1)
            trgbloc[i] = mx[np.argmin(dpdf)]     # Most negative value corresponds to highest rate of decrease 
            
        trgbloc_mean = np.mean(trgbloc)  # Find the TRGB
        trgbloc_sd = np.std(trgbloc)     # Find the Error in the TRGB estimate
        
        return [trgbloc_mean, trgbloc_sd]
    
    def jhcmfind(self, magname = 'Kmag', dmagname = 'eKmag', niter = 1000, kernel = 'epanechnikov'):
    # Set the data to find the TRGB
        mk = self.data.jmag.dropna()-self.data.hmag.dropna()
        dmk = np.sqrt(self.data.jerr.dropna()**2 + self.data.herr.dropna()**2)
        
        mk=-mk
        
            
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
            #kde.score_samples #returns ln(pdf)   # Evaluate the density model on data - probablility density function
            pdf = np.exp(kde.score_samples(mx[:, np.newaxis]))  #MX is x-axis range which the PDF is computed/plotted.
            
            plt.plot(mx,pdf)
            
            #----------------------------------------
            # Set the Edge Detection part using a savgol_filter
            
            smooth_window = 31
            poly_degree = 3
            dpdf = savgol_filter(pdf, smooth_window, poly_degree, deriv = 1)
            trgbloc[i] = mx[np.argmin(dpdf)]     # Most negative value corresponds to highest rate of decrease 
            
        trgbloc_mean = np.mean(trgbloc)  # Find the TRGB
        trgbloc_sd = np.std(trgbloc)     # Find the Error in the TRGB estimate
        
        return [-trgbloc_mean, trgbloc_sd]
    
    def forefind(self, magname = 'Kmag', dmagname = 'eKmag', niter = 1000, kernel = 'epanechnikov',cut=18, up=False):
    # Set the data to find the TRGB
        if self.galaxy=='ngc205' or self.galaxy=='m32':
            
            cut=17
        
        for i in self.data.index:
            
            if self.data.kmag[i] < cut:
                
                self.data.loc[i]=np.nan
    
        mk = self.data.jmag.dropna()-self.data.kmag.dropna()
        dmk = np.sqrt(self.data.jerr.dropna()**2 + self.data.kerr.dropna()**2)
        
        if self.galaxy!='ngc205' and self.galaxy!='m32':
        

            mk=-mk
        
            
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
            #kde.score_samples #returns ln(pdf)   # Evaluate the density model on data - probablility density function
            pdf = np.exp(kde.score_samples(mx[:, np.newaxis]))  #MX is x-axis range which the PDF is computed/plotted.
            
            plt.plot(mx,pdf)
            
            #----------------------------------------
            # Set the Edge Detection part using a savgol_filter
            
            smooth_window = 31
            poly_degree = 3
            dpdf = savgol_filter(pdf, smooth_window, poly_degree, deriv = 1)
            trgbloc[i] = mx[np.argmin(dpdf)]     # Most negative value corresponds to highest rate of decrease 
            
        trgbloc_mean = -np.mean(trgbloc)  # Find the TRGB
        trgbloc_sd = np.std(trgbloc)     # Find the Error in the TRGB estimate
        
        return [trgbloc_mean, trgbloc_sd]
        
    def find_FEH(self):
        
        #function performs conversion between CM and [Fe/H] from Cioni(2009)
        #define c/m ratio
        CM=len(self.cdata.dropna())/len(self.mdata.dropna())
        
        #carry out conversion
        FEH=self.CM_to_FEH(CM)
        #print out result
        print('C/M = ' + str(CM) + ' [Fe/H] = ' + str(FEH))
        

        
        
        
