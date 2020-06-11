import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from data_load import data_load

from data_read import data_read

class data_readall:
    
    def __init__(self,path_to_folder):
        
        self.n147=data_read(path_to_folder + 'ngc147',galaxy='ngc147')
        self.n185=data_read(path_to_folder + 'ngc185',galaxy='ngc185')
        self.n205=data_read(path_to_folder + 'ngc205',galaxy='ngc205')
        self.m32=data_read(path_to_folder + 'm32',galaxy='m32')
        
    def plot_kj_cmds(self):
        
        n147=self.n147
        n185=self.n185
        n205=self.n205
        m32=self.m32
        
        fig,ax=plt.subplots(2,2)
        
    

