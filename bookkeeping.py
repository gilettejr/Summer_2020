from data_load import data_load
from data_read import data_read

class bookkeeping:
    
    def __init__(self):
        
        self.n147=data_load('ngc147')
        self.n185=data_load('ngc185')
        self.n205=data_load('ngc205')
        self.m32=data_load('m32')
        
    
    def update(self):
        
        sets=[self.n147,self.n185,self.n205,self.m32]
        
        for i in sets:
            
            i.forecut()
            i.save_to_parquet('processed_data/fore_cut_data/' + i.galaxy)
            i.trgbcut()
            i.save_to_parquet('processed_data/agb_data/' + i.galaxy)
    

    
    