#basic Python imports
import numpy as np
import matplotlib.pyplot as plt


#libraries for finding trgb
from sklearn import neighbors
from scipy.signal import savgol_filter

#libraries for seperating C and M with non horizontal/perpendicular lines

from data_read import data_read
class edge_detectors(data_read):
    
    
    
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
        
        return [-trgbloc_mean, trgbloc_sd]
    
    def forefind(self, magname = 'Kmag', dmagname = 'eKmag', niter = 1000, kernel = 'epanechnikov',cut=17.5, up=False):
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
            
        trgbloc_mean = -np.mean(trgbloc)  # Find the TRGB
        trgbloc_sd = np.std(trgbloc)     # Find the Error in the TRGB estimate
        
        return [trgbloc_mean, trgbloc_sd]

    
    