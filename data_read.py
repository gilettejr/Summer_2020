from data_load import data_load

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#class inherits all methods from data_loads, but takes data from defined binary file rather than loading from raw data

class data_read(data_load):
    
    def __init__(self,stage='agb',galaxy='ngc147'):
        
        
        if stage=='agb':
            path_to_folder='processed_data/agb_data/'
        
        elif stage=='fore_cut':
            path_to_folder='processed_data/fore_cut_data/'
        
        else:
            path_to_folder=stage
        frame=pd.read_parquet(path_to_folder + galaxy)
        
        galaxies=['ngc147','ngc185','ngc205','m32']
        
        self.data=frame
        self.galaxy=galaxy
        self.galaxies=galaxies
        
