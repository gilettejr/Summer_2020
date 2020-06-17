from data_load import data_load

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#class inherits all methods from data_loads, but takes data from defined binary file rather than loading from raw data

class data_read(data_load):
    
    def __init__(self,stage='cm',galaxy='ngc147'):
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
        
        
        if stage=='agb':
            path_to_folder='processed_data/agb_data/'
            
        elif stage=='cls_cut':
            path_to_folder='processed_data/cls_cut_data/'
        
        elif stage=='fore_cut':
            path_to_folder='processed_data/fore_cut_data/'
            
        elif stage=='cm':
            mpath_to_folder='processed_data/m_agb_data/'
            cpath_to_folder='processed_data/c_agb_data/'
        else:
            path_to_folder=stage
            
        if stage=='cm':
            mframe=pd.read_parquet(mpath_to_folder + galaxy)
            cframe=pd.read_parquet(cpath_to_folder + galaxy)
            
            self.mdata=mframe
            self.cdata=cframe
            self.data=mframe.append(cframe)
        else:
            frame=pd.read_parquet(path_to_folder + galaxy)
            
            self.data=frame


            
        self.galaxy=galaxy
        
        galaxies=['ngc147','ngc185','ngc205','m32']
        
        self.galaxies=galaxies
        
        
        
