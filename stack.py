from data_readall import data_readall
from data_read import data_read
import numpy as np
import pandas as pd

class stack(data_read):
    
    def __init__(self,stage):
        and1=data_read(stage,'and1')
        and2=data_read(stage,'and2')
        and3=data_read(stage,'and3')
        and6=data_read(stage,'and6')
        and7=data_read(stage,'and7')
        and10=data_read(stage,'and10',)
        and14=data_read(stage,'and14',)
        and15=data_read(stage,'and15')
        and16=data_read(stage,'and16')
        and17=data_read(stage,'and17')
        and18=data_read(stage,'and18')
        and19=data_read(stage,'and19')
        and20=data_read(stage,'and20')
        
        sphs=np.array([and1,and2,and3,and6,and7,and10,and14,and15,and16,and17,and18,and19,and20])
        
        

        nums=[1,2,3,6,7,10,14,15,16,17,18,19,20]
        for i in range(len(sphs)):
            
            col=np.full(len(sphs[i].data),nums[i])
            
            sphs[i].data['sph_index']=col
            
        
        frames=[]
        
        for i in sphs:
            frames.append(i.data)
            
        stacks=pd.concat(frames)
        
        print(stacks)
        
        self.data=stacks