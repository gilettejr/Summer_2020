from data_load import data_load

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#class inherits all methods from data_loads, but takes data from defined binary file rather than loading from raw data

class data_read(data_load):
    
    def __init__(self,filename):
        
        frame=pd.read_parquet(filename)
        
        self.frame=frame
        
        
