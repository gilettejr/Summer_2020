#basic Python imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#imports for manipulating astronomical data
from astropy.coordinates import SkyCoord
from astropy.io import ascii
#extinction correction package
from dustmaps.sfd import SFDQuery

#class for processing starting from the raw WFCAM datasets
class data_load:
    #reads in and performs initial cuts on data from chosen galaxy
    #change optional arguments to false to skip initial cuts
    #path_to_file argument used to specify where WFCAM data is stored
    def __init__(self,galaxy, CLS=True, mag=True, ext=True, path_to_files='initial_data/'):
        #not used yet, could be useful later
        #chooses the most recent data based on the highest flag number by default
        #can also specify flag number to choose e.g data removed in foreground cut
        def plotting_select(data,plot='max'):
            
            #default,takes highest flag
            
            if plot=='max':
                #loop finds highest flag (most processed) data, wipes all other data for plotting
                
                for i in range(len(data.kmag)):
                
                    if data.flag[i]!=np.max(data.flag):
                        
                        data.loc[i]=np.nan
            #if specific flag number was chosen in optinal argument, only data with that flag is retained
            
            else:
                
                for i in range(len(data.kmag)):
                    
                
                    if data.flag[i]!=int(plot):
                        
                        data.loc[i]=np.nan
                    
            #returns flag specified data, wipes all NaN values
            
            return(data.dropna())
            
        #class method for use in plotting methods
        
        self.plotting_select=plotting_select
        
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
                
                
        #adds flag column with initial 0 value to DataFrame, for retaining and identifying data from cuts
        
        def make_flag_col(frame):
            
            #column of 0s created
            
            flag=np.zeros(len(frame))
            
            #added to dataframe
            
            frame['flag']=flag
        
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
            
            print('CLS cut reduced data by ' + str(decrease) + '%')
            
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
            
            print('Magerr cut reduced data by ' + str(decrease) + '%')
            
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
        
        infilenames=['lot_n147.unique','lot_n185.unique','NGC205_new_trimmed.unique','M32_new.asc']
        
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
        
        
    def plot_kj_cmd(self,marker='o',markersize=1,color='blue'):
        
        data=self.data
        
        print(data)

        plt.figure()        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.kmag,data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
        
    #method to save data as binary parquet file, to be used by data_read class    
    def save_to_parquet(self,fileloc):
        
        self.frame.to_parquet(fileloc)
        
        
        
        
        
        
        
        