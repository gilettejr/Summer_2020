#basic Python imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#from graphing_utils import graphing_utils
from selection_utils import selection_utils
from initial_processes import initial_processes


#imports for manipulating astronomical data
from astropy.coordinates import SkyCoord
from astropy.io import ascii


#extinction correction package


#curve fitting imports
from scipy.stats import gaussian_kde,norm
from scipy.interpolate import make_interp_spline, BSpline

#libraries for finding trgb
from sklearn import neighbors
from scipy.signal import savgol_filter

#libraries for seperating C and M with non horizontal/perpendicular lines
from shapely.geometry import Point, box, MultiPoint
from shapely.geometry.polygon import Polygon,LinearRing

#libraries for defining ellipses
import shapely.affinity
from descartes import PolygonPatch
#class for processing starting from the raw WFCAM datasets

from matplotlib.patches import Ellipse,Rectangle






# =============================================================================
# __init__ carries out cls, magerr cuts and corrects for reddening, converts
# data into dataframe format. Adds tangent coordinate and decimal ra/dec
#coordinates to dataframe, referenced by xi/etc, RA/DEC respectively

#methods: forecut,trgbcut,trgbfind,plot_kj_cmd,plot_spatial,save_as_parquet
# =============================================================================
class data_loader:
    #reads in and performs initial cuts on data from chosen galaxy
    #change optional arguments to false to skip initial cuts
    #path_to_file argument used to specify where WFCAM data is stored
    def __init__(self,galaxy, CLS=True,cls_bands='norm', mag=True, ext=True, path_to_file='initial_data/'):
        
        
        excolumnames=['rah','ram','ras','dd','dm','ds','unknown','xj','yj','jmag','jerr','jcis','xh','yh','hmag','herr','hcis','xk','yk','kmag','kerr','kcis','u1','u2','u3','u4']
        columnames=['rah','ram','ras','dd','dm','ds','unknown','xj','yj','jmag','jerr','jcis','xh','yh','hmag','herr','hcis','xk','yk','kmag','kerr','kcis']
        
        #galaxy list for reading in the correct file
        
        galaxies=['ngc147','ngc185','ngc205','m32','m31','and1','and2','and3','and6','and7','and10','and14','and15','and16','and17','and18','and19','and20']
        
        #names of datafiles
        
        #galaxy list and galaxy name set as attributes
        #used for matching correct parameters for appropriate
        #galaxy when running class methods
        self.galaxies=galaxies
        self.galaxy=galaxy
        
        
        
        
        #associate galaxy input with appropriate file
        
        for i in range(len(galaxies)):
                
            if galaxy==galaxies[i]:
                    
                file=galaxies[i]
                break
        
        #ascii data read in to astropy table
        #m31 treated differently here because of occasional extra columns
        #end result is the same
        if galaxy=='m31':
        
            galaxy_data=pd.read_csv(path_to_file + file,sep="\s+",names=excolumnames)
            
        else:
            
            galaxy_data=ascii.read(path_to_file + file)
            
            #dataframe created
            
            galaxy_data=galaxy_data.to_pandas()


    
            
            #columns assigned
            
            #extra columns dropped from ngc205 and m32
            
            if self.galaxy=='ngc205' or self.galaxy=='m32':
                
                galaxy_data=galaxy_data.drop(['col26'],axis=1)
                galaxy_data=galaxy_data.drop(['col25'],axis=1)
                galaxy_data=galaxy_data.drop(['col24'],axis=1)
                galaxy_data=galaxy_data.drop(['col23'],axis=1)
            
            #column names assigned
            
            galaxy_data.columns=columnames
            
        self.galaxy_data=galaxy_data

        startup_processor=initial_processes(self)
        
        #decimal and tangent coordinates constructed
        startup_processor.make_deg_coord()
        startup_processor.make_tan_coord()

        
        #optional skip when initiating class for the cuts/extinction corrections
        #all True by default, so mostly will be run
        
        if CLS==True:
            
            
            #cls cut carried out, NaN values purged
            startup_processor.CLS_cut(cls_bands=cls_bands)
            galaxy_data=galaxy_data.dropna()
            print(str(len(galaxy_data)) + ' sources retained after CLS cut')
        if mag==True:
            
            #mag cut carried out, NaN values purged
            startup_processor.mag_err_cut()
            
            galaxy_data=galaxy_data.dropna()
            print(str(len(galaxy_data)) + ' sources retained after mag cut')
        if ext==True:
            
            #extinction corrections done
            startup_processor.ext_corr()
        
        #data set as class attribute
        
        self.data=galaxy_data
    #define line with two points, make a cut on left and right side
    #point arguments given in tuples
    #xdata, ydata lists of colour/magnitude data

        
        
    
    #remove bluer foreground data above a specified blue limit
    def do_forecut(self):
    
    #create copy to hold cut data for completeness
    
        foredata=self.data.copy()
        
        #set class attribute for frame holding data to variable for ease
        
        data=self.data
        
        #set galaxy name list to variable for matching between cut and galaxy
        
        galaxies=self.galaxies
        
        #cuts for each galaxy placed in list. Defined from inspection of
        #j-k CMD
        
        forecuts=[0.992,0.964,1.00,0.92,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,]
        #loop through galaxies to match galaxy with foreground cut
            
        
        for i in range(len(forecuts)):
            
            if self.galaxy==galaxies[i]:
                cut=forecuts[i]
                break
        #j-k colour defined for data
        
        jk=data.jmag-data.kmag
        
        #remove foreground sources in data, retain only foreground sources
        #in foredata
        
        for i in data.index:
            
            if jk[i] < cut:
                
                data.loc[i]=np.nan
                
            else:
                foredata.loc[i]=np.nan
        
        #find % decrease in sources after foreground cut
        #wipe NaN values from DataFrames     
        
        k=len(data.copy())
        data=data.dropna()
        foredata=foredata.dropna()
        decrease=((k-len(data))/k) * 100
        
        print('Foreground cut reduced data by ' + str(int(decrease)) + '%, ' + str(len(data)) + ' sources retained')
        
        #set foreground dataframe as attribute so it can be accessed
        
        self.foredata=foredata
        self.data=data

    #method to cut all data below a certain defined magnitude

    def do_trgbcut(self):
        
        #similar process as in forecut
        #set attributes to variables, for holding data with rgb stars removed
        #and also for holding just rgb stars
        
        data=self.data
        rgbdata=self.data.copy()
        
        #cuts defined from running trgbtip on foreground removed data

        trgbcuts=[18.137,17.862,17.930,17.8,17.8,18.27,18.05,18.43,18.41,18.27,18.27,18.77,18.38,17.92,18.52,19.65,18.75,18.66]
        
        #galaxies attribute used to match galaxy to associated trgb cut
        
        galaxies=self.galaxies
        #by default, uses trgb cut coded into the trgbcuts list
        #optionally allows cut to be included as argument when method is run
        #in which case it will make that cut instead
      
        for i in range(len(galaxies)):
            
            if self.galaxy==galaxies[i]:
                
                cut=trgbcuts[i]
                break
        #make cut on main data attribute, convert rgbdata to frame holding only
        #cut data (rgb stars)        
        
        for i in data.index:
            
            if data.kmag[i] > cut:
                
                data.loc[i]=np.nan
        
            else:
                
                rgbdata.loc[i]=np.nan
                
        #% decrease in sources from cuts printed, NaN values purged
        k=len(data.copy())

        data=data.dropna()

        rgbdata=rgbdata.dropna()
        
        decrease=((k-len(data))/k) * 100
        
        print('RGB cut reduced data by ' + str(int(decrease)) + '%, ' + str(len(data)) + ' sources retained')
        
        #set cut data as attribute for completeness
        
        self.rgbdata=rgbdata
        self.data=data
    
    #separate AGB stars into C and M, using areas defined by boxes in 2D colour space
    #perpendicular to axes
    
    

                
    #separate AGB stars into C-AGB and M-AGB using 2 colour cuts
    #takes cuts from lists inside method, or optional arguments when
    #method is called
    def do_CM_cut(self):
        
        #variables set
        
        data=self.data
        mdata=data.copy()
        cdata=data.copy()
        
        #define cuts in each colour space for each galaxy
        
        #hkcuts=[0.44,0.44,0.60,0.57]
        #jhcuts=[0.82,0.82,0.77,0.93]
        
        hkcuts=[0.337,0.323,0.407,0.477,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        jhcuts=[0.883,0.857,0.930,0.913,0.91,0,0,0,0,0,0,0,0,0,0,0,0,0]
        
        #match galaxy to cut, or make cut from argument if used
        

        
        for i in range(len(self.galaxies)):
            
            if self.galaxies[i]==self.galaxy:
                
                hkcut=hkcuts[i]
                jhcut=jhcuts[i]
                
                break

            
            
        
        #set colour arrays
                
    
        hk=data.hmag-data.kmag
        jh=data.jmag-data.hmag
        
        #loop through and separate data into M and C stars based on colour cuts
        
        for i in data.index:
            
            if hk[i] > hkcut and jh[i] > jhcut:
                
                mdata.loc[i]=np.nan
                
            else:
                
                cdata.loc[i]=np.nan
                

        
        #wipe NaN values
        
        #mdata=mdata.dropna()
        #cdata=cdata.dropna()
        
        #set as class attributes
        self.mdata=mdata
        self.cdata=cdata
    
    #method to remove sources crossmatched with gaia, with well measured
    #parallaxes or proper motions
    
    def gaia_remove(self,path):
        
        #data variable set
        
        data=self.data
        
        
        #csv of crossmatched data produce by topcat read in to DataFrame
        
        cross=pd.read_csv(path)
        
        #remove data from foreground removal routine if Gaia noise measurement
        #is too high
        
        for i in cross.index:
            
            if cross.astrometric_excess_noise_sig[i] > 2:
                
                cross.loc[i]=np.nan
                
        #produce list of indices of sources in original datasets with well measured
        #parallaxes or proper motions
        
        for i in cross.index:
            
            
            if np.abs(cross.pmra[i]/cross.pmra_error[i]) < 0.33 and np.abs(cross.pmdec[i]/cross.pmdec_error[i]) < 0.33:
                
                cross.loc[i]=np.nan
                
        
        #wipe NaN values from index list
        
        cross=cross.dropna()
        
        #remove sources with matching indices to cross list
        #to remove foreground data
        
        
        for i in cross.orig_index:
            
            data.loc[i]=np.nan

        #NaN values wiped
        
        self.data=data.dropna()

        #print(self.data)
    

        
        #plt.errorbar(background.m31_bin_locs,background.bin_densities,yerr=background.bin_uncs,capsize=2,marker=marker,color=color)
        
        

    
    def read_background(self):
        
        self.background=pd.read_parquet('205_background')
        
        #plot graph in new coordinates
            
        #method to save data as binary parquet file, to be used by data_read class    
    def save_to_parquet(self,fileloc):
            
        self.data.to_parquet(fileloc)
        
        #save m and agb dataframes to separate binary files
    def cm_save_to_parquet(self,mfileloc,cfileloc):
        
        self.mdata.to_parquet(mfileloc)
        
        
        self.cdata.to_parquet(cfileloc)
        
        #save data to _csv, extra index column added for bookkeeping
        #when taken in by topcat for crossmatching with Gaia DR2
    def save_to_csv(self,fileloc):
        
        data=self.data
        
        orig_index=data.index
        
        data['orig_index']=orig_index
        
        data.to_csv(fileloc)
        
        

    #28500    
    
        
    def overplot_ellipse(self,ellipticity,a,PA,marker='o',markersize=1,color='black'):
        
        n147ra=[]
        
        data=self.data
        
        b=a*(1-ellipticity)
        ell=Ellipse(xy=[0,0],height=a*2,width=b*2,angle=360-PA,facecolor='none',edgecolor='red',linestyle='--',linewidth=2)
        
        plt.rc('axes',labelsize=20)
        fig,ax=plt.subplots()
        ax.plot(data.xi,data.eta,linestyle='none',marker=marker,markersize=markersize,color=color,zorder=1)
        ax.add_artist(ell)
        ax.invert_xaxis()
        ax.set_ylabel(r'$\eta$')
        ax.set_xlabel(r'$\xi$')
        
    

        
        
        
        
        
        
        
        
                
        
        
        
        
        
        