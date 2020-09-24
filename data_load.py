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
from astropy.modeling import fitting
from astropy.modeling.models import Sersic1D

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
        
        def data_corners(data):
            
            corner1=(np.max(data.xi),np.max(data.eta))
            corner2=(np.max(data.xi),np.min(data.eta))
            corner3=(np.min(data.xi),np.min(data.eta))
            corner4=(np.min(data.xi),np.max(data.eta))
            
            return np.array([corner1,corner2,corner3,corner4])
        self.data_corners=data_corners
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 - 0.47*np.log10(CM)
        
            return(FEH)
            
        self.CM_to_FEH=CM_to_FEH
        
        def CM_to_FEH_unc(CM,deltaCM):
            
            deltalogCM=deltaCM/(CM*np.log(10))
            
            logCM=np.log10(CM)
            
            rhs=0.47 * logCM
            
            deltarhs= rhs * np.sqrt((deltalogCM/logCM)**2 + (0.1/0.47)**2)
            
            deltaFEH=np.sqrt(0.06**2+deltarhs**2)
            
            return deltaFEH
            
        self.CM_to_FEH_unc=CM_to_FEH_unc

        #convert hhmmss ra and ddmmss dec into decimal values
        

        #creates standard coordinates from ra and dec inputs. dE and dSph coordinates are hard coded here    


        #with initial functions defined, we now need to convert the ASCII WFCAM data into a dataframe
        
        #column headers of WFCAM data defined
        
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
    
    
    #REDUNDANT
    #could be useful for some more awkward cutting techniques

    
    #same as above method, but with h-k on x axis
    
    
    #graphing method for plotting spatial distributions
    #can plot m, c, all agb stars, or m and c on separate subplots

    
    #plot contour map of stars

        
        
        
            
    #basically redundant method for finding trgb
    #far less elegant and robust than trgbfind
    

    

    
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
        
        
    #produce C/M ratio and [Fe/H] values from mdata and cdata        
    
    def FEH_find(self):
        
        #function performs conversion between CM and [Fe/H] from Cioni(2009)
        #define c/m ratio
        CM=len(self.cdata.dropna())/len(self.mdata.dropna())
        
        #carry out conversion
        FEH=self.CM_to_FEH(CM)
        #print out result
        print('C/M = ' + str(CM) + ' [Fe/H] = ' + str(FEH))
        
        
        #method to construct elliptical slices and find [Fe/H] in each slices
        #takes semimajor axis of smallest ellipse and outside slice as inputs
    def make_slices(self,stars='agb',a_width=0.02,outer_rad=0.3):
        
        #read in AGB data together, and individual C and M catalogues
        if stars=='agb':
            data=self.data
            
        elif stars=='c':
            
            data=self.cdata
            
        elif stars=='m':
            
            data=self.mdata
            
        self.a_width=a_width
        self.outer_rad=outer_rad
        self.stars=stars
        
        #list of eccentricities of galaxies
        eccentricities=[0.46,0.22,0.43,0.14]
        #list of rotation angles of galaxies
        rotations=[34.2,45.9,169.2,157.9]
        
        for i in range(len(self.galaxies)):
            if self.galaxy==self.galaxies[i]:
                
                eccentricity=eccentricities[i]
                rotation=rotations[i]
        #redundant, for estimating eccentricity

        #initialise selection_utils instance for constructing ellipses
        check=selection_utils()
        #list for holding data slices
        slices=[]
        areas=[]
        ellipse_shapes=[]
        #fill slices with elliptical selections of decreasing size, from outer_rad to a_width
        for i in range(int((outer_rad*1000)/(a_width*1000))):
        
            ellipse=check.select_ellipse(data,afl=outer_rad-(a_width * i),eccentricity=eccentricity,clockrot=rotation)
            slices.append(ellipse[0])
            areas.append(ellipse[1])
            ellipse_shapes.append(ellipse[2])
        
        self.ellipticity=eccentricity
            
            
            
        #
        #remove overlapping data between ellipses, leaving only elliptical slices in the list
        #
        


        
        
        for i in range(len(slices)-1):
            
            areas[i]=areas[i]-areas[i+1]
            intersection=ellipse_shapes[i].intersection(ellipse_shapes[i+1])
            ellipse_shapes[i]=ellipse_shapes[i].difference(intersection)
            #looks ugly but far faster than using more nested loops
            
            #produce list of common points between outer and inner ellipse at each element in slices list
            locs=np.where(slices[i].orig_index==slices[i+1].orig_index)
            #np.where produces extra array for some reason
            #remove that
            locs=locs[0]
            #match index in array to original pandas index
            indices=[]
            for j in locs:
                indices.append(slices[i].index[j])
                
            #wipe all common values in outer slice
            for k in indices:
                slices[i].loc[k]=np.nan
            
            

        for i in slices:
        
            i=i.dropna()

        self.areas=areas
        self.slices=slices
        self.slice_shapes=ellipse_shapes

            
        
            
            
        #now match slices to m and c datasets to determine c/m ratio


        
            
            
        
    def FEH_slices(self,m_background,c_background,crowding_num):
        
        #read in AGB data together, and individual C and M catalogues
        cdata=self.cdata
        mdata=self.mdata
        slices=self.slices
        outer_rad=self.outer_rad
        a_width=self.a_width
        
        
        
        #list of eccentricities of galaxies



        mslices=[]
        cslices=[]
        
        #need to create copies to prevent dataframes being
        #continuously altered
        
        for i in slices:
            
            mslices.append(i.copy())
            cslices.append(i.copy())
        
        #find locations of m/c data in slices
        
        for i in range(len(slices)):
            
            clocs=np.where(mslices[i].orig_index==mdata.orig_index)
            mlocs=np.where(cslices[i].orig_index==cdata.orig_index)
            
            #match array index to pandas index
            
            mindices=[]
            cindices=[]
            
            for j in mlocs:
                
                mindices.append(cslices[i].index[j])
            
            for j in clocs:
                
                cindices.append(mslices[i].index[j])
                
            #wipe C stars from m data slices
                
            for k in mindices:
                
                mslices[i].loc[k]=np.nan
            
            #wipe M stars from c data slices

            for k in cindices:
                
                cslices[i].loc[k]=np.nan


            mslices[i]=mslices[i].dropna()
            cslices[i]=cslices[i].dropna()
            
        #construct new lists to find c/m ratio
            
        mnum=[]
        cnum=[]            
        
        #find number of c and m stars in each slices
        
        for i in range(len(slices)):
            mnum.append(len(mslices[i]))
            cnum.append(len(cslices[i]))
            
        mnum=np.array(mnum)
        cnum=np.array(cnum)
        
        #create arrray of c/m ratios in each slices
        
        cm=cnum/mnum
        
        c_unc=np.sqrt(cnum)
        m_unc=np.sqrt(mnum)
        
        back_num=c_background[2]/m_background[2]
        
        back_cm=c_background[0]/m_background[0]
        
        slice_all=cnum+mnum
        
        all_stars=np.sum(slice_all)
        
        weights=(slice_all/all_stars)
        


        
        back_m_uncs=np.sqrt(m_background[2])
        back_c_uncs=np.sqrt(c_background[2])
        
        back_cm_unc= back_cm * np.sqrt((back_m_uncs/m_unc)**2 + (back_c_uncs/c_unc)**2)
        
        
        cm_slice_unc= cm * np.sqrt((c_unc/cnum)**2 + (m_unc/mnum)**2)
        
        cm_slice_unc= np.sqrt(cm_slice_unc**2 + back_cm_unc**2)
        
        cm=cm-(back_cm * (back_num/all_stars))
        
        avg_cm=np.average(cm[:(len(cm)-crowding_num)],weights=weights[:(len(cm)-crowding_num)])
        
        avg_cm_unc=np.sqrt(np.sum(cm_slice_unc**2))/len(cm_slice_unc)
        
        print(avg_cm)
        
        print(avg_cm_unc)

        

        
        #convert to [Fe/H]
        
        FEH = self.CM_to_FEH(cm)
        
        
        FEH_unc=self.CM_to_FEH_uncs(cm,cm_slice_unc)
        

        
        avgFEH=self.CM_to_FEH(avg_cm)
        
        avgFEH_unc=self.CM_to_FEH_uncs(avg_cm,avg_cm_unc)
        
        print(avgFEH)
        
        print(avgFEH_unc)
        
        #plot radial distribution
        
        xdata=np.linspace(outer_rad-a_width/2,0+a_width/2,num=(outer_rad*1000)/(a_width*1000))
        #plt.errorbar(xdata,FEH,yerr=FEH_unc,capsize=2,linestyle='none',marker='o',markersize='3',color='black')
        m,b=np.polyfit(xdata,FEH,1)
        #plt.plot(xdata, m*xdata + b,color='red')
        
        return xdata, FEH, FEH_unc, avgFEH,avgFEH_unc
        
        #plt.xlabel('a/deg')
        #plt.ylabel('[Fe/H]')


        
             
    def slice_count_profile(self,background_deg,crowding_num=0):
        
        print(background_deg)
        
        slices=self.slices
        areas=self.areas
        outer_rad=self.outer_rad
        a_width=self.a_width
        ellipticity=self.ellipticity
        slice_nums=[]
        
        for i in slices:
            slice_nums.append(len(i.dropna()))
            
        slice_nums=np.array(slice_nums)

        areas=np.array(areas)
        
        slice_count_uncs=np.sqrt(slice_nums)



        slice_star_densities=slice_nums/areas
        slice_star_densities_uncs=(np.sqrt((slice_count_uncs/areas)**2 + background_deg[1]**2))

        slice_star_densities=slice_star_densities-background_deg[0]
        
        slice_star_densities=slice_star_densities/3600
        slice_star_densities_uncs=slice_star_densities_uncs/3600
        
        
        
        xdata=np.linspace(outer_rad-a_width/2,0+a_width/2,num=(outer_rad*1000)/(a_width*1000))
        xdata=xdata*60
        #convert from a to r_eff
        ydata=slice_star_densities
        yerr=slice_star_densities_uncs

        
        #plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',marker='o',markersize=5,capsize=2,color='black',label='AGB data')
        #plt.xlabel('Semi-major axis/arcmins')
        #plt.ylabel('Nstars/arcmins$^2$')
        
        #n147 params
        
        xdata_orig=xdata
        ydata_orig=ydata
        yerr_orig=yerr
            
        xdata=xdata[:len(xdata)-(crowding_num)]
        ydata=ydata[:len(ydata)-(crowding_num)]
        yerr=yerr[:len(yerr)-(crowding_num)]
        
        
        
        if self.galaxy=='ngc147':
        
            model=Sersic1D(r_eff=5.2,n=1.8)
            
        elif self.galaxy=='ngc185':
            
            model=Sersic1D(r_eff=1.7,n=1.00)
            
            
        
        fit=fitting.LevMarLSQFitter()
        
        s=fit(model,xdata,ydata,weights=1/yerr)
        def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
            '''
            fit (array) values for the fit
            x,y,yerr (arrays) data
            N total number of points
            n_free number of parameters we are fitting
            '''
            return 1.0/(N-n_free)*sum(((fit - y)/yerr)**2)
    
        print(calc_reduced_chi_square(s(xdata),xdata,ydata,yerr=yerr,N=len(ydata),n_free=3))
        
        xdata=np.linspace(min(xdata_orig),max(xdata_orig),1000)
        
        #plt.plot(xdata,s(xdata),label='Sersic Fit with r$_{eff}$ = ' + str(round(s.r_eff.value,3)) + ', n = ' + str(round(s.n.value,3)))
        
        #plt.legend()
        print(s)
        self.r_eff=s.r_eff.value
        self.n=s.n.value
        
        return xdata_orig,ydata_orig,yerr_orig,xdata,s(xdata)
        
        #plt.yscale('log')
        #r_eff=a*(1-n/3)
        
        #print('r_eff = ' + str(r_eff)+ ' arcmins')
        
        
        
    def close_slice_count_profile(self):
        
        print('Make sure youve got a background dataframe for subtraction!')
        #remake alpha coordinates
        slice_shapes=self.slice_shapes
        background=self.background
        rotate_coords=self.rotate_coords
        stars=self.stars
        if stars=='agb':
            
            data=self.data
        elif stars=='m':
            
            
            data=self.mdata
            
        elif stars=='c':
            
            data=self.cdata
            
        areas=self.areas
        slices=self.slices
        bin_shapes=self.bin_shapes
        
        outer_rad=self.outer_rad
        a_width=self.a_width
        
        m31ra=10.68470833
        m31dec=41.26875
        
        if self.galaxy=='ngc205':
            tra=10.09189356
            tdec=41.68541564
        elif self.galaxy=='m32':
            tra=10.6742708
            tdec=40.8651694
        
        m31xi,m31eta=self.eq_to_tan(m31ra,m31dec,tra,tdec)
        

        
        #find angle of rotation
        theta=np.arctan((m31eta)/(m31xi))
        
        xi=data.xi.copy()
        eta=data.eta.copy()
        theta_deg=np.degrees(theta)
        #rotate
        
        
        
        for i in slices:
            i['alpha']=rotate_coords(xi,eta,theta)[0]
            i['beta']=rotate_coords(xi,eta,theta)[1]
        slice_shapes_rot=[]
        for i in slice_shapes:
        
            slice_shapes_rot.append(shapely.affinity.rotate(i,-theta_deg))
        
        
        slice_nums=[]
        
        for i in slices:
            
            slice_nums.append(len(i.dropna()))
            
        slice_nums=np.array(slice_nums)
        areas=np.array(areas)
        
        slice_densities_deg=slice_nums/areas
        slice_densities_unc=np.sqrt(slice_nums)
        slice_densities_min=slice_densities_deg/3600
        slice_densities_unc=(slice_densities_unc/areas)/3600
        #background columns: bin_densities, bin_uncs, bin_locs
        #areas in degrees for ease, background in mins
        new_density=[]
        new_uncs=[]
        for i in range(len(slice_shapes_rot)):
            slice_bin_densities=[]
            slice_bin_densities_uncs=[]
            slice_bin_areas=[]
            for j in background.index:
                

                
                if slice_shapes_rot[i].intersects(bin_shapes[j])==True:
                    
                    overlap_shape=slice_shapes_rot[i].intersection(bin_shapes[j])
                    overlap_area=overlap_shape.area
                    
                    
                    #subtract background
                    overlap_density=slice_densities_min[i]-background.density_fit[j]
                    overlap_density_unc=np.sqrt(slice_densities_unc[i]**2 + background.bin_uncs[j]**2)
                    
                    slice_bin_densities.append(overlap_density)
                    slice_bin_densities_uncs.append(overlap_density_unc)
                    slice_bin_areas.append(overlap_area)
                    
            slice_bin_areas=np.array(slice_bin_areas)
            slice_bin_densities=np.array(slice_bin_densities)
            
            weightings=slice_bin_areas/slice_shapes_rot[i].area
            
            avg_density=np.average(slice_bin_densities,weights=weightings)
            density_err=(np.sqrt(np.sum(np.square(slice_bin_densities_uncs))))/len(slice_bin_densities)
            new_density.append(avg_density)
            new_uncs.append(density_err)

            
        xdata=np.linspace(outer_rad-a_width/2,0+a_width/2,num=(outer_rad*1000)/(a_width*1000))
        xdata=xdata*60
            
        ydata=new_density
        yerr=new_uncs
            
        fig=plt.figure()
        

        
        density_distribution=pd.DataFrame({'a':xdata,'density':ydata,'density_err':yerr,'slice_area_deg':areas,'slice_nums':slice_nums})
        

            
        outfilename=self.galaxy+'_radial_density_' + stars
        
        density_distribution.to_parquet(outfilename)
            
        plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',capsize=2,marker='o',color='black')
            
    def close_FEH_slices(self):
        
        if self.galaxy=='ngc205':
            
            cut=5
            
        elif self.galaxy=='m32':
            
            cut=3
        
        c_slice_dataframe=pd.read_parquet(self.galaxy + '_radial_density_c')
        
        m_slice_dataframe=pd.read_parquet(self.galaxy +'_radial_density_m')
        
        full_slice_nums=c_slice_dataframe.slice_nums + m_slice_dataframe.slice_nums
        
        total_star_num=np.sum(full_slice_nums)

        slice_a=c_slice_dataframe.a
        
        c_slice_dataframe
        
        cm_slices=c_slice_dataframe.density/m_slice_dataframe.density
        print(m_slice_dataframe.density)
        print(c_slice_dataframe.density)



        cm_slices_uncs=cm_slices*np.sqrt((c_slice_dataframe.density_err/c_slice_dataframe.density)**2 + (m_slice_dataframe.density_err/m_slice_dataframe.density)**2)

        FEH_slices=self.CM_to_FEH(cm_slices)
        #print(cm_slices)
        #print(cm_slices_uncs)
        FEH_slice_uncs=self.CM_to_FEH_uncs(cm_slices,cm_slices_uncs)

        
        plt.errorbar(slice_a,FEH_slices,marker='o',linestyle='none',color='black',capsize=2,yerr=FEH_slice_uncs)
        plt.ylabel('[Fe/H]')
        plt.xlabel('Semi-major axis (arcmins)')
        
        weights=full_slice_nums[cut:]/total_star_num

        
        
        
        

        
        new_cm_slices=[]
        new_cm_slices_uncs=[]
        new_weights=[]
        for i in range(len(cm_slices)):
            
            if m_slice_dataframe.density[i]<0 or c_slice_dataframe.density[i]<0:
                
                continue
            
            else:
                
                new_cm_slices.append(cm_slices[i])
                new_cm_slices_uncs.append(cm_slices_uncs[i])
                new_weights.append(weights[i])

        avg_CM=np.average(new_cm_slices,weights=new_weights)        
        avg_CM_unc=np.sqrt(np.sum(np.square(new_cm_slices_uncs)))/len(new_cm_slices_uncs)
        
        FEH_avg=self.CM_to_FEH(avg_CM)
        
        FEH_avg_unc=self.CM_to_FEH_uncs(avg_CM,avg_CM_unc)
        print(avg_CM)
        print(avg_CM_unc)
        print(FEH_avg)
        print(FEH_avg_unc)
        
    
        
        
            
    def fit_close_profiles(self,crowding_num=1):
        stars='agb'
        distribution=pd.read_parquet(self.galaxy + '_radial_density_' + stars)
        
        xdata=distribution.a
        ydata=distribution.density
        yerr=distribution.density_err
        
        
        plt.errorbar(xdata,ydata,yerr=yerr,marker='o',linestyle='none',color='black')
        
        plt.ylabel('Density (N/arcmins$^2$)')
        plt.xlabel('Semi-major axis (arcmins)')
        
        xdata_orig=xdata
            
        xdata=xdata[:len(xdata)-(crowding_num)]
        ydata=ydata[:len(ydata)-(crowding_num)]
        yerr=yerr[:len(yerr)-(crowding_num)]
        
        
        
        if self.galaxy=='ngc205':
        
            model=Sersic1D(r_eff=2.0,n=1.5)
            
        elif self.galaxy=='m32':
            
            model=Sersic1D(r_eff=1.7,n=1.00)
            
            
        
        fit=fitting.LevMarLSQFitter()
        
        s=fit(model,xdata,ydata,weights=1/yerr)

        xdata=np.linspace(min(xdata_orig),max(xdata_orig),num=1000)                   
                    
                    
        plt.plot(xdata,s(xdata),label='Sersic Fit with a$_{eff}$ = ' + str(round(s.r_eff.value,3)) + ', n = ' + str(round(s.n.value,3)))
        plt.legend()
                
                
                
 
        
        #calculate density in each slice again
        
        #find overlapping area of each bin on each slice
        
        #make appropriate density reduction for each overlapping segment
        
        #average out slice density
        
        
    def slice_mag_profile(self):
        
        def mag_to_flux(m,zpflux):
            
            f=zpflux*10**(-m/2.5)
            
            return f
        
        def flux_to_mag(f,zpflux):
            
            m=2.5*np.log10(zpflux/f)
            
            return m
            


        
        
        jzpoint=1556.8
        hzpoint=1038.3
        kzpoint=644.1
        

        
        slices=self.slices
        areas=self.areas
        outer_rad=self.outer_rad
        a_width=self.a_width
        
        slice_fluxes=[]
        
        for i in slices:
            slice_fluxes.append(np.sum(mag_to_flux(i.jmag,jzpoint)))
            
        areas=np.array(areas)

        
        
        slice_fluxes=np.array(slice_fluxes)
            

        
        slice_fluxes_deg=slice_fluxes/areas
        slice_fluxes_sec=slice_fluxes_deg/(3600^2)
        
        slice_mags=flux_to_mag(slice_fluxes_sec,jzpoint)
        
        xdata=np.linspace(outer_rad-a_width/2,0+a_width/2,num=(outer_rad*1000)/(a_width*1000))
        ydata=slice_mags
        
        plt.plot(xdata,ydata,marker='o',color='black')
        plt.gca().invert_yaxis()
        plt.xlabel('Semi-major axis/degrees')
        plt.ylabel('Surface Brightness in Jmag/arcsec$^2$')
        

        #make ellipses
        
        #cut redundant points
        
        #calculate C/M, FEH
        
    def find_background_density_border(self,stars='agb',marker='o',markersize='1',color='black',borderwidth=0.07,show_figure=False):
        
        if stars=='agb':
            
            data=self.data
            
        elif stars=='m':
            
            data=self.mdata
            
        elif stars =='c':
            
            data=self.cdata
        
        border_data=data.copy()
        corner_locs=self.data_corners(self.data)
        
        for i in corner_locs:
            
            i=np.asarray(i)
            
        
        subtraction_identity=np.array([[-1,-1],[-1,1],[1,1],[1,-1]])
        
        spatial_subtraction_matrix=borderwidth * subtraction_identity
        
        inner_corner_locs=corner_locs + spatial_subtraction_matrix
            
            
        inner_corner_locs=map(tuple,inner_corner_locs)
        inner_corner_locs=tuple(inner_corner_locs)
        
        select_border=selection_utils()
        select_border.select_stars(border_data,inner_corner_locs[0],inner_corner_locs[2],unselect=True)
        
        inner_area=(inner_corner_locs[0][0]-inner_corner_locs[2][0]) * (inner_corner_locs[0][1]-inner_corner_locs[2][1])
        outer_area=(corner_locs[0][0]-corner_locs[2][0]) * (corner_locs[0][1]-corner_locs[2][1])
        

        
        
        
        border_area=outer_area-inner_area

        
        
        border_num=len(border_data.dropna())
        
        border_err=np.sqrt(border_num)
        
        
        if show_figure==True:
            
            borders=[Rectangle(xy=inner_corner_locs[2],width=inner_corner_locs[0][0]-inner_corner_locs[2][0],height=inner_corner_locs[0][1]-inner_corner_locs[2][1],fill=False)]
            
            
            plt.rc('axes',labelsize=20)
            fig,ax=plt.subplots()
            ax.plot(data.xi,data.eta,linestyle='none',marker=marker,markersize=markersize,color=color,zorder=1)
            for i in borders:
                
                ax.add_artist(i)
            ax.invert_xaxis()
            ax.set_ylabel(r'$\eta$')
            ax.set_xlabel(r'$\xi$')
            
         
            
        return np.array([border_num/border_area,border_err/border_area,border_num]) 
    
    def find_background_density_boxes(self,marker='o',markersize='1',color='black',boxsize=0.2,show_figure=False):
        
        data=self.data
        c1dat=data.copy()
        c2dat=data.copy()
        c3dat=data.copy()
        c4dat=data.copy()
        
        #set corners
        
        corner1=(np.max(data.xi),np.max(data.eta))
        corner2=(np.max(data.xi),np.min(data.eta))
        corner3=(np.min(data.xi),np.min(data.eta))
        corner4=(np.min(data.xi),np.max(data.eta))
        
        corner1a=(np.max(data.xi)-boxsize,np.max(data.eta)-boxsize)
        corner2a=(np.max(data.xi)-boxsize,np.min(data.eta)+boxsize)
        corner3a=(np.min(data.xi)+boxsize,np.min(data.eta)+boxsize)
        corner4a=(np.min(data.xi)+boxsize,np.max(data.eta)-boxsize)
        
        c1=selection_utils()
        corner1s=c1.select_stars(c1dat,corner1,corner1a)
        
        c2=selection_utils()
        corner2s=c2.select_stars(c2dat,corner2,corner2a)
        
        c3=selection_utils()
        corner3s=c3.select_stars(c3dat,corner3,corner3a)
        
        c4=selection_utils()
        corner4s=c4.select_stars(c4dat,corner4,corner4a)
        
        
        corner_dats=[c1dat.dropna(),c2dat.dropna(),c3dat.dropna(),c4dat.dropna()]
        corner_nums=[]
        
        for i in corner_dats:
            corner_nums.append(len(i))
        

        area=boxsize**2
        
        corner_num=np.sum(np.array(corner_nums))
        
        corner_p=corner_num/(area*4)
        

        
        
        if show_figure==True:
        
            squares=[Rectangle(xy=corner3,width=boxsize,height=boxsize,fill=False),Rectangle(xy=corner1a,width=boxsize,height=boxsize,fill=False),Rectangle(xy=corner2s[1],width=boxsize,height=boxsize,fill=False),Rectangle(xy=corner4s[3],width=boxsize,height=boxsize,fill=False)]
            
            plt.rc('axes',labelsize=20)
            fig,ax=plt.subplots()
            ax.plot(data.xi,data.eta,linestyle='none',marker=marker,markersize=markersize,color=color,zorder=1)
            for i in squares:
                
                ax.add_artist(i)
            ax.invert_xaxis()
            ax.set_ylabel(r'$\eta$')
            ax.set_xlabel(r'$\xi$')
            
        
        return corner_p
    
    def find_background_grad(self,ellipse_a=0.15,ellipticity=0.43,clockrot=169.2,marker='o',markersize='1',color='black',show_figure=False,binwidth=0.02):
        
        if self.galaxy=='m32':
            
            ellipticity=0.14
            clockrot=157.9
            ellipse_a=0.1
        
        def outer_boundary(data,x,y):
            
            corner1loc=np.where(data==np.min(data.alpha))[0][0]
            corner2loc=np.where(data==np.min(data.beta))[0][0]
            corner3loc=np.where(data==np.max(data.alpha))[0][0]
            corner4loc=np.where(data==np.max(data.beta))[0][0]
            
            corner1loc=data.index[data['alpha'] == np.min(data.alpha)].tolist()[0]
            corner2loc=data.index[data['beta'] == np.min(data.beta)].tolist()[0]
            corner3loc=data.index[data['alpha'] == np.max(data.alpha)].tolist()[0]
            corner4loc=data.index[data['beta'] == np.max(data.beta)].tolist()[0]
            
            
            
            corner1=(data.alpha[corner1loc],data.beta[corner1loc])
            corner2=(data.alpha[corner2loc],data.beta[corner2loc])
            corner3=(data.alpha[corner3loc],data.beta[corner3loc])
            corner4=(data.alpha[corner4loc],data.beta[corner4loc])
            
            boundary=LinearRing([corner1,corner2,corner3,corner4])
            
            return boundary
        
        stars=self.stars
        
        if stars=='agb':
            
            data=self.data
            
        elif stars=='m':
            
            data=self.mdata
            
        elif stars=='c':
            
            data=self.cdata
        

        
        #wipe galaxy
        
        e=selection_utils()
        
        dE_stars=e.select_ellipse(data,afl=ellipse_a,eccentricity=ellipticity,clockrot=clockrot,unselect=True)
        dE_ellipse=dE_stars[2]
        
        data=dE_stars[0]

        
        #convert m31 coords to local tangent coords

        
        m31ra=10.68470833
        m31dec=41.26875
        
        if self.galaxy=='ngc205':
            tra=10.09189356
            tdec=41.68541564
        elif self.galaxy=='m32':
            tra=10.6742708
            tdec=40.8651694
        
        m31xi,m31eta=self.eq_to_tan(m31ra,m31dec,tra,tdec)
        

        
        #find angle of rotation
        theta=np.arctan((m31eta)/(m31xi))
        
        print(theta)
        
        xi=data.xi.copy()
        eta=data.eta.copy()

        #rotate
        data['alpha']=self.rotate_coords(xi,eta,theta)[0]
        data['beta']=self.rotate_coords(xi,eta,theta)[1]
        
        m31xi=np.degrees(m31xi)
        m31eta=np.degrees(m31eta)
        
        m31_alpha=self.rotate_coords(m31xi,m31eta,theta)[0]
        m31_beta=self.rotate_coords(m31xi,m31eta,theta)[1]

        
        if show_figure==True:
        
            plt.plot(data.alpha,data.beta,marker=marker,markersize=markersize,linestyle='none',color=color)
            plt.gca().invert_xaxis()
            plt.gca().set_ylabel(r'$\beta$')
            plt.gca().set_xlabel(r'$\alpha$')
        
        corner1=np.array([np.min(data.xi),np.min(data.eta)])
        corner2=np.array([np.min(data.xi),np.max(data.eta)])
        corner3=np.array([np.max(data.xi),np.max(data.eta)])
        corner4=np.array([np.max(data.xi),np.min(data.eta)])
        
        corners=[corner1,corner2,corner3,corner4]
        corners_rot=[]
        for i in corners:
            
            corners_rot.append(self.rotate_coords(i[0],i[1],theta))
        
        
        


        border=LinearRing(corners_rot)
        s=Polygon(border)
        corners_x=np.array([corners_rot[0][0],corners_rot[1][0],corners_rot[2][0],corners_rot[3][0]])
        
        upper_bound=np.max(corners_x)-binwidth
        lower_bound=np.min(corners_x)
        
        print(upper_bound)
        print(lower_bound)
    
        
        #make bins
        


        #bin1
        
        bin_nums=[]
        bin_areas=[]
        bin_locs=[]
        #define masked ellipse
        mask=dE_ellipse
        theta_deg=np.degrees(theta)
        mask=shapely.affinity.rotate(dE_ellipse,-theta_deg)
        x,y=mask.exterior.xy

        bin_shapes=[]
        while upper_bound-lower_bound > binwidth:
            
        
            bin_slice_fill=box(upper_bound-binwidth,-1.0,upper_bound,1.0)

            bin_slice=bin_slice_fill.exterior
        
            #find intercepts
            bin_corners=border.intersection(bin_slice)
            
            swap0=bin_corners[0]
            swap1=bin_corners[1]
            swap2=bin_corners[3]
            swap3=bin_corners[2]
            
            corners=[swap0,swap1,swap2,swap3]
            bin_corners=MultiPoint([swap0,swap1,swap2,swap3])
            final_bin=Polygon(bin_corners)
            bin_shapes.append(final_bin)
            x,y=final_bin.exterior.xy
            if show_figure==True:
                
                plt.plot(x,y,color='orange')

            final_bin_area=final_bin.area
            
            box_data=data.copy()
            
            for i in box_data.index:
                
                if final_bin.contains(Point(box_data.alpha[i],box_data.beta[i]))==False:
                    
                    box_data.loc[i]=np.nan
                    
            bin_nums.append(len(box_data.dropna()))
            
            #subtract area overlapping with ellipse
            
            if bin_slice_fill.intersects(mask)==True:
                
                overlap_area=mask.intersection(bin_slice_fill).area

                final_bin_area=final_bin_area-overlap_area
            
            bin_locs.append(upper_bound-(binwidth/2))
            bin_areas.append(final_bin_area)
            upper_bound=upper_bound-binwidth
            
            if self.galaxy=='m32' and upper_bound-lower_bound < (binwidth * 2):
                
                break
            
        bin_nums=np.array(bin_nums)
        bin_locs=np.array(bin_locs)
        bin_areas=np.array(bin_areas)
        
        
        
        bin_uncs=np.sqrt(bin_nums)
        
        bin_densities=(bin_nums/bin_areas)/3600
        
        bin_uncs=(bin_uncs/bin_areas)/3600
        
        bin_locs=bin_locs * 60
        
        m31_alpha=m31_alpha*60

        m31_bin_locs=bin_locs-m31_alpha

        background=pd.DataFrame({'bin_densities':bin_densities,'bin_uncs':bin_uncs,'bin_locs':bin_locs,'m31_bin_locs':m31_bin_locs})
        
        self.background=background
        background.to_parquet('205_c_tests')
        
        self.bin_shapes=bin_shapes
        #fit with Sersic profile
        
        #plt.errorbar(background.m31_bin_locs,background.bin_densities,yerr=background.bin_uncs,capsize=2,marker=marker,color=color)
        

        
    def fit_close_background(self,marker='o',markersize='3',color='black'):
        
        #background=pd.read_parquet('205_c_tests')
    
        
        #stars='agb'
        
        background =self.background.copy()
        
        
        xdata=-background.m31_bin_locs
        
        ydata=background.bin_densities
        
        yerr=background.bin_uncs
        

        
        model=Sersic1D(r_eff=2,n=10.00)
            
            
        fit=fitting.LevMarLSQFitter()
        
        if self.galaxy=='ngc205' and self.stars=='c':
            

            model=Sersic1D(r_eff=10,n=20)
            s=fit(model,xdata[2:],ydata[2:])
            
            
        elif self.galaxy=='m32' and self.stars=='c':
            model=Sersic1D(r_eff=20,n=10)
            s=fit(model,xdata[3:],ydata[3:])
        elif self.galaxy=='m32':
            
            model=Sersic1D(r_eff=20,n=2)
            s=fit(model,xdata[3:],ydata[3:],weights=1/yerr[3:])
        
        else:
            
            s=fit(model,xdata[2:],ydata[2:],weights=1/yerr[2:])
        
        plt.errorbar(xdata,ydata,yerr=yerr,capsize=2,marker=marker,color=color)
        
        plt.plot(xdata,s(xdata))
        
        background['density_fit']=s(xdata)
        
        self.background=background
        
        print(background)
        

    
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
    
        
        
    def plot_cross(self,file='crossmatching/spitzer/m32_cross'):
        
        data=self.data
        
        crossdata=pd.read_csv(file)
        
        plotter=graphing_utils()
        plotter.plot_kj_cmd(data)
        plotter.plot_kj_cmd(crossdata,color='red',overlay=True)
        
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
        
    
    def overplot_ellipses(self,ellipticity,a_inner,a_outer,PA,marker='o',markersize=1,color='black'):
        
        data=self.data
        majors=np.linspace(a_inner,a_outer,num=int((a_outer*1000)/(a_inner*1000)))
        minors=majors*(1-ellipticity)
        
        ells=[]
        
        for i in range(len(majors)):
            ells.append(Ellipse(xy=[0,0],height=majors[i]*2,width=minors[i]*2,angle=360-PA,facecolor='none',edgecolor='red',linestyle='--',linewidth=1))
            
        plt.rc('axes',labelsize=20)
        fig,ax=plt.subplots()
        ax.plot(data.xi,data.eta,linestyle='none',marker=marker,markersize=markersize,color=color,zorder=1)
        for i in ells:
            ax.add_artist(i)
        ax.invert_xaxis()
        ax.set_ylabel(r'$\eta$')
        ax.set_xlabel(r'$\xi$')
        
        
        
        
        
        
        
        
                
        
        
        
        
        
        