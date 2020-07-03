from data_load import data_load
from data_read import data_read
from crossmatch_stilts import crossmatch
from edge_detectors import edge_detectors

import numpy as np

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

        def edge(stage,galaxy,axis='none'):

            
            locs=[]
            sigs=[]
            

                
            
            edge=edge_detectors(stage=stage,galaxy=galaxy)
            
            for i in range(3):
                if stage=='cls_cut':
                
                    result=edge.forefind()
                
                elif stage=='fore_cut':
                    
                    result=edge.trgbfind()
                    
                elif stage=='agb_crossed' or stage=='agb':
                    
                    if axis=='hk':
                    
                        result=edge.hkcmfind()
                        
                    elif axis=='jh':
                        
                        result=edge.jhcmfind()
                
                locs.append(result[0])
                sigs.append(result[1])
                
            locs=np.array(locs)
            sigs=np.array(sigs)
            
            loc=np.mean(locs)
            sig=np.mean(sigs)
            
            loc=np.abs(loc)
            
            return [loc,sig]
            
        self.edge=edge

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
            
        crosssets=self.read() 
        for i in crosssets:
            i.save_to_csv('crossmatching/ukirt_pre/' + i.galaxy)
            crossmatch(i.galaxy)
        
        for i in crosssets:
            i.gaia_remove('crossmatching/gaia/' + i.galaxy)
            i.save_to_parquet('processed_data/agb_crossed_data/' + i.galaxy)
            i.CM_cut()
            i.cm_save_to_parquet('processed_data/m_agb_data/' + i.galaxy,'processed_data/c_agb_data/' + i.galaxy)
            
    def edge_update(self):
        
        sets=self.load()
        edge=self.edge
        for i in sets:
            i.save_to_parquet('processed_data/cls_cut_data/' + i.galaxy)
            forebord=edge(stage='cls_cut',galaxy=i.galaxy,axis='none')
            print(i.galaxy)
            print(forebord[0])
            print(forebord[1])
            i.forecut(cut=forebord[0])
            i.save_to_parquet('processed_data/fore_cut_data/' + i.galaxy)
            trgbord=edge(stage='fore_cut',galaxy=i.galaxy,axis='none')
            print(trgbord[0])
            print(trgbord[1])
            i.trgbcut(cut=trgbord[0])
            i.save_to_parquet('processed_data/agb_data/' + i.galaxy)
            hkbord=edge(stage='agb_crossed',galaxy=i.galaxy,axis='hk')
            print(hkbord[0])
            print(hkbord[1])
            jhbord=edge(stage='agb_crossed',galaxy=i.galaxy,axis='jh')
            print(jhbord[0])
            print(jhbord[1])
            i.CM_cut(hkcut=hkbord[0],jhcut=jhbord[0])
            i.cm_save_to_parquet('processed_data/m_agb_data/' + i.galaxy,'processed_data/c_agb_data/' + i.galaxy)

            
    

    
    