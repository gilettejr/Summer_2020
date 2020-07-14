#basic Python imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from graphing_utils import graphing_utils
from selection_utils import selection_utils

#imports for manipulating astronomical data
from astropy.coordinates import SkyCoord
from astropy.io import ascii

#extinction correction package
from dustmaps.sfd import SFDQuery

#curve fitting imports
from scipy.stats import gaussian_kde,norm
from scipy.interpolate import make_interp_spline, BSpline

#libraries for finding trgb
from sklearn import neighbors
from scipy.signal import savgol_filter

#libraries for seperating C and M with non horizontal/perpendicular lines
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#libraries for defining ellipses
import shapely.affinity
from descartes import PolygonPatch
#class for processing starting from the raw WFCAM datasets
# =============================================================================
# __init__ carries out cls, magerr cuts and corrects for reddening, converts
# data into dataframe format. Adds tangent coordinate and decimal ra/dec
#coordinates to dataframe, referenced by xi/etc, RA/DEC respectively

#methods: forecut,trgbcut,trgbfind,plot_kj_cmd,plot_spatial,save_as_parquet
# =============================================================================
class data_load:
    #reads in and performs initial cuts on data from chosen galaxy
    #change optional arguments to false to skip initial cuts
    #path_to_file argument used to specify where WFCAM data is stored
    def __init__(self,galaxy, CLS=True,CLS_mags='all', mag=True, ext=True, path_to_file='initial_data/'):
        
        def CM_to_FEH(CM):
        
            FEH=-1.39 -0.47*np.log10(CM)
        
            return(FEH)
            
        self.CM_to_FEH=CM_to_FEH

        def linecut(frame,xdata,ydata,point1,point2):
        
            upper=frame.copy()
            lower=frame.copy()
            
            m=(point2[1]-point1[1])/(point2[0]-point1[0])
            
            c=point1[1]-point1[0] * m
            
            for i in frame.index:
                
                if ydata[i] > m * xdata[i] + c:
                    
                    lower.loc[i]=np.nan
                    
                else:
                    
                    upper.loc[i]=np.nan
            
            return([upper,lower])
        
        self.linecut=linecut
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
        #creates standard coordinates from ra and dec inputs. dE and dSph coordinates are hard coded here    
        def make_tan_coord(frame,galaxy):
                #intermediate function to do the heavy lifting
                def create_tangent_coords(frame,tangentra,tangentdec):
        
                    #ra and dec attributes converted to variables in radians
        
                    ra = np.radians(frame.RA)
                    dec = np.radians(frame.DEC)
                    
                    #tangent co-ordinates also converted to radians
                    
                    tanra = np.radians(tangentra)
                    tandec = np.radians(tangentdec)
                    
                    #conversion for xi carried out
                    
                    xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
                    
                    #conversion for eta carried out
                    
                    eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
                    
                    #co-ordinates converted to degrees and set as attributes
                    
                    xi = xi * (180/np.pi)
                    eta = eta * (180/np.pi)
                    
                    return (np.array([xi,eta]))
                #ra and dec positions in decimal degrees for each galaxy
                if galaxy=='ngc147':
                    tra=8.300500
                    tdec=48.50850
                elif galaxy=='ngc185':
                    tra=9.7415417
                    tdec=48.3373778
                elif galaxy=='ngc205':
                    tra=10.09189356
                    tdec=41.68541564
                elif galaxy=='m32':
                    tra=10.6742708
                    tdec=40.8651694
                
                elif galaxy=='m31':
                    tra=10.68470833
                    tdec=41.26875
                    
                elif galaxy=='and1':
                    tra=11.41583333
                    tdec=38.04111111
                    
                elif galaxy=='and2':
                    tra=19.12416667
                    tdec=33.41916667
                    
                elif galaxy=='and3':
                    tra=8.84458333
                    tdec=36.50472222
                    
                elif galaxy=='and6':
                    tra=357.94333333
                    tdec=24.58638889
                    
                elif galaxy=='and7':
                    tra=351.63208333
                    tdec=50.67583333
                    
                elif galaxy=='and10':
                    tra=16.64041667
                    tdec=44.80438889
                    
                    
                elif galaxy=='and14':
                    tra=12.89583333
                    tdec=29.69694444
                    
                elif galaxy=='and15':
                    tra=18.57791667
                    tdec=38.1175
                    
                elif galaxy=='and16':
                    tra=14.87416667
                    tdec=32.37666667
                    
                elif galaxy=='and17':
                    tra=9.27916667
                    tdec=44.32222222
                    
                elif galaxy=='and18':
                    tra=0.56041667
                    tdec=45.08888889
                    
                    
                elif galaxy=='and19':
                    tra=4.88375
                    tdec=35.04363889
                    
                elif galaxy=='and20':
                    tra=1.8779166700000003
                    tdec=35.13233333
                #function run dependent on galaxy
                coords=create_tangent_coords(frame,tra,tdec)
                #output produced and xi and eta columns added to dataframe
                frame['xi']=coords[0]
                frame['eta']=coords[1]
                    
        #function to cut DataFrame based on cls index
        
        def CLS_cut(frame,bands='norm'):
            #set of cuts used for normal data processing
            if bands=='norm':
                #m32 treated differently to prevent gradient issues in spatial distribution
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                        #standard cls cuts made for the rest
                        if (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0 and frame.jcis[i]!=-3.0) or (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0 and frame.kcis[i]!=-3.0) or (frame.hcis[i]==-8.0) or frame.jmag[i]-frame.hmag[i] > 17:
                            frame.loc[i]=np.nan
                
                else:
                

                #standard cls cuts made for the rest
                    for i in range(len(frame.jcis)):
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0) or (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
            #different combinations of cls cuts included below
            #can be used to test effect of different cls cuts

                
            elif bands=='jh' or bands=='hj':
                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0 and frame.jcis[i]!=-3.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0 and frame.hcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0) or (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0): 
                            frame.loc[i]=np.nan      
                            
            elif bands=='hk':
                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0 and frame.kcis[i]!=-3.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0 and frame.hcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
                            
            elif bands =='jk':
                

                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0 and frame.kcis[i]!=-3.0) or (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0 and frame.jcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0) or (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
                            
            elif bands == 'j':
                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0 and frame.jcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
                            
            elif bands == 'h':
                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0 and frame.hcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
                            
            elif bands == 'k':
                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0 and frame.kcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
                            
            elif bands=='all':
                
                if self.galaxy=='m32':
                    for i in range(len(frame.jcis)):
                    
                        if (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0 and frame.jcis[i]!=-3.0) or (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0 and frame.kcis[i]!=-3.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0 and frame.hcis[i]!=-3.0):
                            frame.loc[i]=np.nan
                
                else:
                
                #removes any non-stellar sources from data
                    for i in range(len(frame.jcis)):
                        if (frame.kcis[i] != -1.0 and frame.kcis[i]!=-2.0) or (frame.hcis[i] != -1.0 and frame.hcis[i]!=-2.0) or (frame.jcis[i] != -1.0 and frame.jcis[i]!=-2.0): 
                            frame.loc[i]=np.nan
                            
            else:
                #failure if invalid cls argument chosen upon class initialisation
                print('Invalid CLS argument chosen, no CLS cut made')
            #percentage decrease of catalogue length calculated and printed
            
            k=len(frame.jcis)
            l=len(frame.jcis.dropna())            
            decrease=100-(l/k)*100
            
            print('CLS cut reduced data by ' + str(int(decrease)) + '%')
            
            #for some reason if I put frame=frame.dropna() here it doesn't do anything, so I have to wipe the NaNs outside the function
            
            
        #function to cut data with large magnitude error
        #takes dataframe as input, only used internally in this class
        
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
        #again used in the class only
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
        
            frame=pd.read_csv(path_to_file + file,sep="\s+",names=excolumnames)
            
        else:
            
            frame=ascii.read(path_to_file + file)
            
            #dataframe created
            
            frame=frame.to_pandas()


    
            
            #columns assigned
            
            #extra columns dropped from ngc205 and m32
            
            if self.galaxy=='ngc205' or self.galaxy=='m32':
                
                frame=frame.drop(['col26'],axis=1)
                frame=frame.drop(['col25'],axis=1)
                frame=frame.drop(['col24'],axis=1)
                frame=frame.drop(['col23'],axis=1)
            
            #column names assigned
            
            frame.columns=columnames

        
        #decimal and tangent coordinates constructed
        make_deg_coord(frame)
        make_tan_coord(frame,self.galaxy)
        
        #optional skip when initiating class for the cuts/extinction corrections
        #all True by default, so mostly will be run
        
        if CLS==True:
            
            
            #cls cut carried out, NaN values purged
            CLS_cut(frame,bands=CLS_mags)
            frame=frame.dropna()
            print(str(len(frame)) + ' sources retained after CLS cut')
        if mag==True:
            
            #mag cut carried out, NaN values purged
            mag_err_cut(frame)
            
            frame=frame.dropna()
            print(str(len(frame)) + ' sources retained after mag cut')
        if ext==True:
            
            #extinction corrections done
            ext_corr(frame)
        
        #data set as class attribute
        
        self.data=frame
    #define line with two points, make a cut on left and right side
    #point arguments given in tuples
    #xdata, ydata lists of colour/magnitude data

        
        
    
    #remove bluer foreground data above a specified blue limit
    def forecut(self,cut=0):
        
        #create copy to hold cut data for completeness
        
        foredata=self.data.copy()
        
        #set class attribute for frame holding data to variable for ease
        
        data=self.data
        
        #set galaxy name list to variable for matching between cut and galaxy
        
        galaxies=self.galaxies
        
        #cuts for each galaxy placed in list. Defined from inspection of
        #j-k CMD
        
        forecuts=[0.992,0.964,1.00,0.92,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,1.002,]
        foresigs=[0.103,0.079,0.085,0.01]
        #loop through galaxies to match galaxy with foreground cut
        
        if cut==0:
            
        
            for i in range(len(forecuts)):
                
                if self.galaxy==galaxies[i]:
                    cut=forecuts[i]
                    break
            
        else:
            
            cut=cut
        
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
        
    def forecut_slant(self,slant_m='+'):
        
        #create copy to hold cut data for completeness
        
        foredata=self.data.copy()
        
        #set class attribute for frame holding data to variable for ease
        
        data=self.data
        jk=data.jmag-data.kmag
        k=data.kmag
        
        #set galaxy name list to variable for matching between cut and galaxy
        
        galaxies=self.galaxies
        
        #cuts for each galaxy placed in list. Defined from inspection of
        #j-k CMD
        
        oneforecuts=[(0.95,17.9),(0.71,18.35),(0,0),(0,0)]
        twoforecuts=[(0.99,17.0),(0,0),(0,0),(0,0)]
        
        #loop through galaxies to match galaxy with foreground cut
        
        for i in range(len(oneforecuts)):
            
            if self.galaxy==galaxies[i]:
                cut1=oneforecuts[i]
                cut2=twoforecuts[i]
                break
        
        result=self.linecut(data,jk,k,cut1,cut2)
        
        if slant_m=='+':
            foredata=result[1]
            data=result[0]
            
        elif slant_m=='-':
            foredata=result[0]
            data=result[1]
        else:
            print('Gradient argument incorrect in forecut_slant method')
            
        
        
        #set foreground dataframe as attribute so it can be accessed
        
        self.foredata=foredata.dropna()
        self.data=data.dropna()

    #method to cut all data below a certain defined magnitude

    def trgbcut(self,cut=0):
        
        #similar process as in forecut
        #set attributes to variables, for holding data with rgb stars removed
        #and also for holding just rgb stars
        
        data=self.data
        rgbdata=self.data.copy()
        
        #cuts defined from running trgbtip on foreground removed data
        

        
        trgbcuts=[18.137,17.862,17.930,17.8,17.8,18.27,18.05,18.43,18.41,18.27,18.27,18.77,18.38,17.92,18.52,19.65,18.75,18.66]
        trgbsigs=[0.141,0.110,0.066,0.67]
            
 

        
        #galaxies attribute used to match galaxy to associated trgb cut
        
        galaxies=self.galaxies
        #by default, uses trgb cut coded into the trgbcuts list
        #optionally allows cut to be included as argument when method is run
        #in which case it will make that cut instead
        if cut==0:        
            for i in range(len(galaxies)):
                
                if self.galaxy==galaxies[i]:
                    
                    cut=trgbcuts[i]
                    break
                
        else:
            
            cut=cut        
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
    def CM_cut(self,hkcut=0,jhcut=0):
        
        #variables set
        
        data=self.data
        mdata=data.copy()
        cdata=data.copy()
        
        #define cuts in each colour space for each galaxy
        
        #hkcuts=[0.44,0.44,0.60,0.57]
        #jhcuts=[0.82,0.82,0.77,0.93]
        
        hkcuts=[0.337,0.323,0.407,0.477]
        hksigs=[0.047,0.052,0.034,0]
        jhcuts=[0.883,0.857,0.930,0.913,0.91]
        jhsigs=[0.046,0.042,0.049,0]
        
        #match galaxy to cut, or make cut from argument if used
        
        if hkcut == 0 or jhcut == 0:
        
            for i in range(len(self.galaxies)):
                
                if self.galaxies[i]==self.galaxy:
                    
                    hkcut=hkcuts[i]+hksigs[i]
                    jhcut=jhcuts[i]+jhsigs[i]
                    
                    break
                
        else:
            
            hkcut=hkcut
            jhcut=jhcut
            
            
        
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
    def CM_cut_slant(self,slant_m='-'):
        
        #create copy to hold cut data for completeness
        
        mdata=self.data.copy()
        
        #set class attribute for frame holding data to variable for ease
        
        cdata=self.data.copy()
        data=self.data.copy()
        jh=data.jmag-data.hmag
        hk=data.hmag-data.kmag
        
        #set galaxy name list to variable for matching between cut and galaxy
        
        galaxies=self.galaxies
        
        #cuts for each galaxy placed in list. Defined from inspection of
        #j-k CMD
        
        oneforecuts=[(1.09,0.13),(0,0),(0,0),(0,0)]
        twoforecuts=[(0.78,0.48),(0,0),(0,0),(0,0)]
        
        #loop through galaxies to match galaxy with foreground cut
        
        for i in range(len(oneforecuts)):
            
            if self.galaxy==galaxies[i]:
                cut1=oneforecuts[i]
                cut2=twoforecuts[i]
                break
        
        result=self.linecut(data,jh,hk,cut1,cut2)
        
        if slant_m=='+':
            mdata=result[0]
            cdata=result[1]
            
        elif slant_m=='-':
            mdata=result[1]
            cdata=result[0]
        else:
            print('Gradient argument incorrect in forecut_slant method')
            
        
        
        #set foreground dataframe as attribute so it can be accessed
        
        self.mdata=mdata.dropna()
        self.cdata=cdata.dropna()
    #separate C and M stars based on defined triangle in 2D colour space
    #REDUNDANT
    def CM_polygon_cut(self):
        
        #create copies for holding separated data
        
        data=self.data
        mdata=self.data.copy()
        cdata=self.data.copy()
        
        #define triangle holding C-stars for each galaxy
        
        vertex1s=[(0.75,0.53),(0,0),(0,0),(0,0)]
        vertex2s=[(1.06,0.21),(0,0),(0,0),(0,0)]
        vertex3s=[(2.25,1.31),(0,0),(0,0),(0,0)]
        vertex4s=[(2.25,2.25),(0,0),(0,0),(0,0)]
        vertex5s=[(1.65,2.25),(0,0),(0,0),(0,0)]
        
        
        #match galaxy to appropriate triangle
        
        for i in range(len(self.galaxies)):
            
            if self.galaxies[i]==self.galaxy:
                
                vertex1=vertex1s[i]
                vertex2=vertex2s[i]
                vertex3=vertex3s[i]
                vertex4=vertex4s[i]
                vertex5=vertex5s[i]
                    
                break
        
        #colours defined
        
        hk=data.hmag-data.kmag
        
        jh=data.jmag-data.hmag
        
        #create triangle between 3 defined vertices surrounding C-stars
        
        carea=Polygon([vertex1,vertex2,vertex3,vertex4,vertex5])
        
        #loop through data, allocate data inside triangle to cdata dataframe
        #remainder held in mdata dataframe
        
        for i in data.index:
            
            if carea.contains(Point(jh[i],hk[i])):
                
                mdata.loc[i]=np.nan
                
            else:
                
                cdata.loc[i]=np.nan
                
        #wipe NaN values, set data to class attributes
        
        mdata=mdata.dropna()
        cdata=cdata.dropna()
        
        self.mdata=mdata
        self.cdata=cdata
    
    #graphing method for plotting j-k cmds
    
    def plot_kj_cmd(self,stars='all',marker='o',markersize=1,color='black',newfig=False):
        
        #conditional statements plot only c,m, or both sets depending on 
        #optional stars argument
        
        if stars=='c':
            
            data=self.cdata
            
        elif stars=='m':
            
            data=self.mdata
        
        elif stars=='c+m' :
        
            data=self.mdata.append(self.cdata)
            
        else:
            
            data=self.data
            
        
        
        #axes, figure set, CMD plotted
        
        if newfig==True: 
        
            plt.figure()        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.kmag,data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        if color!='red':
            plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
    
    #same as above method, but with h-k on x axis
    def plot_hk_cmd(self,stars='all',marker='o',markersize=1,color='blue'):
        
        #conditional statements plot only c,m, or both sets depending on 
        #optional stars argument
        
        if stars=='c':
            
            data=self.cdata
            
        elif stars=='m':
            
            data=self.mdata
        
        elif stars=='c+m' :
        
            data=self.mdata.append(self.cdata)
            
        else:
            
            data=self.data
            
        
        
        #axes, figure set, CMD plotted

        #plt.figure()        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.hmag - data.kmag,data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        if color!='red':
            plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$H_0$-$K_0$')
    
    #graphing method for plotting h-k/j-h 2CD
    
    def plot_cc(self,stars='all',marker='o',markersize=1,color='black',newfig=False):
        
        #stars argument works in same way as plot_kj_cmd
        
        if stars=='c':
            
            data=self.cdata
            
        elif stars=='m':
            
            data=self.mdata
        
        elif stars=='c+m' :
        
            data=self.mdata.append(self.cdata)
            
        else:
            
            data=self.data
            
        
        #axes, figure set, data plotted
        
        if newfig==True:
        
            plt.figure()        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.hmag,data.hmag-data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)

        plt.ylabel('$H_0$-$K_0$')
        plt.xlabel('$J_0$-$H_0$')
        
    
    #graphing method for plotting spatial distributions
    #can plot m, c, all agb stars, or m and c on separate subplots
    def plot_spatial(self,stars='all',marker='o',markersize=1,color='black'):
    #conditional statement based on stars defines data appropriately
        if stars=='c':
            
            data=self.cdata
            
        elif stars=='m':
            
            data=self.mdata
        
        elif stars=='c+m' :
        
            mdata=self.mdata
            cdata=self.cdata
            
        else:
            
            data=self.data
        
        #if c+m cnot chosen for stars argument, data plotted in single plot
        
        if stars !='c+m':
        
            plt.rc('axes',labelsize=20)
            plt.plot(data.xi,data.eta,linestyle='none',marker=marker,markersize=markersize,color=color)
            plt.gca().set_ylabel(r'$\eta$')
            plt.gca().set_xlabel(r'$\xi$')
            
            if color!='red':
                plt.gca().invert_xaxis()
                    
                plt.show()
        
        
        #if c+m chosen for stars argument, plot c and m stars in subplots
        else:
            
            #shared axes for clarity and ease of comparison
            
            fig,axs=plt.subplots(2,1,sharex=True,sharey=True)
            
            plt.rc('axes',labelsize=20)

            
            axs[0].plot(mdata.xi,mdata.eta,linestyle='none',marker=marker,markersize=markersize,color='blue')
            axs[1].plot(cdata.xi,cdata.eta,linestyle='none',marker=marker,markersize=markersize,color='red')
            
            axs[0].set_ylabel(r'$\eta$')
            axs[1].set_ylabel(r'$\eta$')
            
            axs[1].set_xlabel(r'$\xi$')
            
            axs[0].set_title('M-AGB')
            axs[1].set_title('C-AGB')
            
            axs[0].invert_xaxis()
            axs[1].invert_xaxis()
    
    #plot contour map of stars
    def plot_contour(self,overlay=False):
        
        #set graphing visuals
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'18',
        'axes.titlesize':'14','xtick.labelsize':'12',
        'ytick.labelsize':'12','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        #read in dataframe
        data=self.data
        #make plot
        sns.kdeplot(data.xi,data.eta,levels=np.logspace(-1,1,50))
        #don't reinvert x axis and relabel if being used as an overlay
        if overlay==False:
            
            
            
            plt.ylabel(r'$\eta$')
            plt.xlabel(r'$\xi$')
            plt.gca().invert_xaxis()
        
        
        
        
            
    #basically redundant method for finding trgb
    #far less elegant and robust than trgbfind
    
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
    def FEH_slices(self,a_width=0.03,outer_rad=0.3):
        
        #read in AGB data together, and individual C and M catalogues
        cdata=self.cdata
        mdata=self.mdata
        data=self.data
        #list of eccentricities of galaxies
        eccentricities=[0.44,0,0,0]
        #list of rotation angles of galaxies
        rotations=[56,0,0,0]
        #redundant, for estimating eccentricity
        atup=(0.075,0.089)
        btup=(0.032,-0.042)
        #initialise selection_utils instance for constructing ellipses
        check=selection_utils()
        #list for holding data slices
        slices=[]
        
        #fill slices with elliptical selections of decreasing size, from outer_rad to a_width
        for i in range(int((outer_rad*1000)/(a_width*1000))):
        
            
            slices.append(check.select_ellipse(self.data,afl=outer_rad-(a_width * i),eccentricity=eccentricities[0],clockrot=rotations[0]))
            
            
            
        #
        #remove overlapping data between ellipses, leaving only elliptical slices in the list
        #
        
        
        
        for i in range(len(slices)-1):
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
                

                        
            #slices[i]=slices[i].dropna()
            
            



        mslices=[]
        cslices=[]
        
        for i in slices:
            
            mslices.append(i.copy())
            cslices.append(i.copy())
        

        
        for i in range(len(slices)):
            
            clocs=np.where(mslices[i].orig_index==mdata.orig_index)
            mlocs=np.where(cslices[i].orig_index==cdata.orig_index)
            

            
            mindices=[]
            cindices=[]
            
            for j in mlocs:
                
                mindices.append(cslices[i].index[j])
            
            for j in clocs:
                
                cindices.append(mslices[i].index[j])
                

                
            for k in mindices:
                
                mslices[i].loc[k]=np.nan


            for k in cindices:
                
                cslices[i].loc[k]=np.nan


            mslices[i]=mslices[i].dropna()
            cslices[i]=cslices[i].dropna()
            

            
        mnum=[]
        cnum=[]            
        
        for i in range(len(slices)):
            mnum.append(len(mslices[i]))
            cnum.append(len(cslices[i]))
            
        mnum=np.array(mnum)
        cnum=np.array(cnum)
        
        
        cm=cnum/mnum
        
        FEH= self.CM_to_FEH(cm)
        
        xdata=np.linspace(outer_rad,0,num=(outer_rad*1000)/(a_width*1000))
        plt.plot(xdata,FEH,linestyle='none',marker='o',markersize='3',color='black')
        m,b=np.polyfit(xdata,FEH,1)
        plt.plot(xdata, m*xdata + b,color='red')
        plt.xlabel('a/deg')
        plt.ylabel('[Fe/H]')


        
             
        
        
        
        #make ellipses
        
        #cut redundant points
        
        #calculate C/M, FEH
        
            
            
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
        
        
        
        
        
        
                
        
        
        
        
        
        