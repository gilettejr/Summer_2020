from data_load import data_loader
from data_read import data_reader
from crossmatch_stilts import crossmatch
from edge_detectors import edge_detectors

import numpy as np
import os

#class to keep track of everything in intermediate data directories
class bookkeeping:
    #performs extinction, cls and magerr cuts for each galaxy,
    def __init__(self,cls_bands='norm'):
        
        def load_all():
            n147=data_loader('ngc147',cls_bands=cls_bands)
            n185=data_loader('ngc185',cls_bands=cls_bands)
            n205=data_loader('ngc205',cls_bands=cls_bands)
            m32=data_loader('m32',cls_bands=cls_bands)
            
            return [n147,n185,n205,m32]
        
        def read_all(stage='agb'):
            n147=data_reader(stage=stage,galaxy='ngc147')
            n185=data_reader(stage=stage,galaxy='ngc185')
            n205=data_reader(stage=stage,galaxy='ngc205')
            m32=data_reader(stage=stage,galaxy='m32')
            
            return [n147,n185,n205,m32]
        
        def load_sph_all():
            and1=data_loader('and1',cls_bands=cls_bands)
            and2=data_loader('and2',cls_bands=cls_bands)
            and3=data_loader('and3',cls_bands=cls_bands)
            and6=data_loader('and6',cls_bands=cls_bands)
            and7=data_loader('and7',cls_bands=cls_bands)
            and10=data_loader('and10',cls_bands=cls_bands)
            and14=data_loader('and14',cls_bands=cls_bands)
            and15=data_loader('and15',cls_bands=cls_bands)
            and16=data_loader('and16',cls_bands=cls_bands)
            and17=data_loader('and17',cls_bands=cls_bands)
            and18=data_loader('and18',cls_bands=cls_bands)
            and19=data_loader('and19',cls_bands=cls_bands)
            and20=data_loader('and20',cls_bands=cls_bands)
            
            return [and1,and2,and3,and6,and7,and10,and14,and15,and16,and17,and18,and19,and20]
            
        def read_sph_all(stage='agb'):
            and1=data_reader(stage,'and1')
            and2=data_reader(stage,'and2')
            and3=data_reader(stage,'and3')
            and6=data_reader(stage,'and6')
            and7=data_reader(stage,'and7')
            and10=data_reader(stage,'and10',)
            and14=data_reader(stage,'and14',)
            and15=data_reader(stage,'and15')
            and16=data_reader(stage,'and16')
            and17=data_reader(stage,'and17')
            and18=data_reader(stage,'and18')
            and19=data_reader(stage,'and19')
            and20=data_reader(stage,'and20')
            
            return [and1,and2,and3,and6,and7,and10,and14,and15,and16,and17,and18,and19,and20]
            

            
        self.load_all=load_all
        self.read_all=read_all
        
        self.load_sph_all=load_sph_all
        self.read_sph_all=read_sph_all
        

    def update_from_scratch(self,galaxy_type='des'):
        
        #create list of objects to run through all galaxies
        if galaxy_type=='des':
            
            
            all_galaxy_datasets=self.load_all()
            
        elif galaxy_type=='dsphs':
            
            all_galaxy_datasets=self.load_sph_all()
       
        #carry out processing for each galaxy, saving to binary file after
        #each step
        for galaxy_object in all_galaxy_datasets:
            
            try:
            
                galaxy_object.save_to_parquet('processed_data/cls_cut_data/' + galaxy_object.galaxy)
                
            except:
                
                os.system('mkdir processed_data')
                os.system('mkdir processed_data/cls_cut_data')
                galaxy_object.save_to_parquet('processed_data/cls_cut_data/' + galaxy_object.galaxy)                
                
                
            galaxy_object.do_forecut()
            try:
                
                galaxy_object.save_to_parquet('processed_data/fore_cut_data/' + galaxy_object.galaxy)
            except:
                os.system('mkdir processed_data/fore_cut_data')
                galaxy_object.save_to_parquet('processed_data/fore_cut_data/' + galaxy_object.galaxy)
            galaxy_object.trgbcut()
            try:
                galaxy_object.save_to_parquet('processed_data/agb_data/' + galaxy_object.galaxy)
            except:
                os.system('mkdir processed_data/agb_data')
                galaxy_object.save_to_parquet('processed_data/agb_data/' + galaxy_object.galaxy)
            
        all_galaxy_datasets=self.read_all()
        for galaxy_object in all_galaxy_datasets:
            try:
                galaxy_object.save_to_csv('crossmatching/ukirt_pre/' + galaxy_object.galaxy)
            except:
                os.system('mkdir crossmatching')
                os.system('mkdir crossmatching/ukirt_pre')
                os.system('mkdir crossmatching/gaia')
                galaxy_object.save_to_csv('crossmatching/ukirt_pre/' + galaxy_object.galaxy)
                
            
            crossmatch(galaxy_object.galaxy)

        
        for galaxy_object in all_galaxy_datasets:
            galaxy_object.gaia_remove('crossmatching/gaia/' + galaxy_object.galaxy)
            try:
                galaxy_object.save_to_parquet('processed_data/agb_crossed_data/' + galaxy_object.galaxy)
            except:

                os.system('mkdir processed_data/agb_crossed_data')
                os.system('mkdir processed_data/m_agb_data')
                os.system('mkdir processed_data/c_agb_data')
                galaxy_object.save_to_parquet('processed_data/agb_crossed_data/' + galaxy_object.galaxy)
            galaxy_object.CM_cut()
            galaxy_object.cm_save_to_parquet('processed_data/m_agb_data/' + galaxy_object.galaxy,'processed_data/c_agb_data/' + galaxy_object.galaxy)
        
            
    def update_dsphs(self):
        

            
        sets=self.load_sph()
       
        #carry out processing for each galaxy, saving to binary file after
        #each step
        for i in sets:
            i.save_to_parquet('processed_data/cls_cut_data/' + i.galaxy)
            i.do_forecut()
            i.save_to_parquet('processed_data/fore_cut_data/' + i.galaxy)
            i.trgbcut()
            i.save_to_parquet('processed_data/agb_data/' + i.galaxy)
            
        crosssets=self.read_sph(stage='agb') 
        for i in crosssets:
            i.save_to_csv('crossmatching/ukirt_pre/' + i.galaxy)
            crossmatch(i.galaxy)
        
        for i in crosssets:
            i.gaia_remove('crossmatching/gaia/' + i.galaxy)
            i.save_to_parquet('processed_data/agb_crossed_data/' + i.galaxy)
            

        
    

            
    

    
    