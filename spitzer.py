import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
class spitzer:
    
    def __init__(self,file='final_M32_var.csv'):
        
        sdata=pd.read_csv(file)
        
        print(sdata)