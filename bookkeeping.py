from data_load import data_load




#class to keep track of everything in intermediate data directories
class bookkeeping:
    #performs extinction, cls and magerr cuts for each galaxy,
    def __init__(self):
        
        self.n147=data_load('ngc147')
        self.n185=data_load('ngc185')
        self.n205=data_load('ngc205')
        self.m32=data_load('m32')
        
    #update class runs through all the processing steps and saves all
    #intermediate data from the most up to date parameters
    def update(self):
        
        #create list of objects to run through all galaxies
        sets=[self.n147,self.n185,self.n205,self.m32]
        
        
        #carry out processing for each galaxy, saving to binary file after
        #each step
        for i in sets:
            
            i.forecut()
            i.save_to_parquet('processed_data/fore_cut_data/' + i.galaxy)
            i.trgbcut()
            i.save_to_parquet('processed_data/agb_data/' + i.galaxy)
            i.CM_cut()
            i.cm_save_to_parquet('processed_data/m_agb_data/' + i.galaxy,'processed_data/c_agb_data/' + i.galaxy)
    

    
    