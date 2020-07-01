from data_load import data_load
from data_read import data_read



#class to keep track of everything in intermediate data directories
class bookkeeping:
    #performs extinction, cls and magerr cuts for each galaxy,
    def __init__(self,CLS_mags='all'):
        
        def load():
            n147=data_load('ngc147',CLS_mags=CLS_mags)
            n185=data_load('ngc185',CLS_mags=CLS_mags)
            n205=data_load('ngc205',CLS_mags=CLS_mags)
            m32=data_load('m32',CLS_mags=CLS_mags)
            
            return [n147,n185,n205,m32]
        
        def read(stage='agb'):
            n147=data_read(stage=stage,galaxy='ngc147')
            n185=data_read(stage=stage,galaxy='ngc185')
            n205=data_read(stage=stage,galaxy='ngc205')
            m32=data_read(stage=stage,galaxy='m32')
            
            return [n147,n185,n205,m32]
            
        self.load=load
        self.read=read

    #update class runs through all the processing steps and saves all
    #intermediate data from the most up to date parameters
    def update(self):
        
        #create list of objects to run through all galaxies

        sets=self.load()
        
        #carry out processing for each galaxy, saving to binary file after
        #each step
        for i in sets:
            i.save_to_parquet('processed_data/cls_cut_data/' + i.galaxy)
            i.forecut()
            i.save_to_parquet('processed_data/fore_cut_data/' + i.galaxy)
            i.trgbcut()
            i.save_to_parquet('processed_data/agb_data/' + i.galaxy)
            i.CM_cut()
            i.cm_save_to_parquet('processed_data/m_agb_data/' + i.galaxy,'processed_data/c_agb_data/' + i.galaxy)
            
    def prep_cross(self):
        
        sets=self.read()
        
        for i in sets:
            i.save_to_csv('crossmatching/ukirt_pre/' + i.galaxy)
    
    def crossmatch_update(self):
        
        for i in self.sets:
            i.gaia_remove('crossmatching/gaia/' + i.galaxy)
            
    

    
    