from data_load import data_load
from data_read import data_read
from crossmatch_stilts import crossmatch
from edge_detectors import edge_detectors

import numpy as np
import os

#class to keep track of everything in intermediate data directories
class bookkeeping:
    #performs extinction, cls and magerr cuts for each galaxy,
    def __init__(self,CLS_mags='norm'):
        
        def load_all():
            n147=data_load('ngc147',CLS_mags=CLS_mags)
            n185=data_load('ngc185',CLS_mags=CLS_mags)
            n205=data_load('ngc205',CLS_mags=CLS_mags)
            m32=data_load('m32',CLS_mags=CLS_mags)
            
            return [n147,n185,n205,m32]
        
        def read_all(stage='agb'):
            n147=data_read(stage=stage,galaxy='ngc147')
            n185=data_read(stage=stage,galaxy='ngc185')
            n205=data_read(stage=stage,galaxy='ngc205')
            m32=data_read(stage=stage,galaxy='m32')
            
            return [n147,n185,n205,m32]
        
        def load_sph_all():
            and1=data_load('and1',CLS_mags=CLS_mags)
            and2=data_load('and2',CLS_mags=CLS_mags)
            and3=data_load('and3',CLS_mags=CLS_mags)
            and6=data_load('and6',CLS_mags=CLS_mags)
            and7=data_load('and7',CLS_mags=CLS_mags)
            and10=data_load('and10',CLS_mags=CLS_mags)
            and14=data_load('and14',CLS_mags=CLS_mags)
            and15=data_load('and15',CLS_mags=CLS_mags)
            and16=data_load('and16',CLS_mags=CLS_mags)
            and17=data_load('and17',CLS_mags=CLS_mags)
            and18=data_load('and18',CLS_mags=CLS_mags)
            and19=data_load('and19',CLS_mags=CLS_mags)
            and20=data_load('and20',CLS_mags=CLS_mags)
            
            return [and1,and2,and3,and6,and7,and10,and14,and15,and16,and17,and18,and19,and20]
            
        def read_sph_all(stage='agb'):
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
            
            return [and1,and2,and3,and6,and7,and10,and14,and15,and16,and17,and18,and19,and20]
            

            
        self.load_all=load_all
        self.read_all=read_all
        
        self.load_sph_all=load_sph_all
        self.read_sph_all=read_sph_all
        


        def edge_finder(stage,galaxy,axis='none'):

            
            edge=edge_detectors(stage=stage,galaxy=galaxy)
            

            if stage=='cls_cut':
            
                edge_location_array=edge.forefind()
            
            elif stage=='fore_cut':
                
                edge_location_array=edge.trgbfind()
                
            elif stage=='agb_crossed' or stage=='agb':
                
                if axis=='hk':
                
                    edge_location_array=edge.hkcmfind()
                    
                elif axis=='jh':
                    
                    edge_location_array=edge.jhcmfind()
                
            edge_location=np.abs(edge_location_array[0])
            edge_location_sig=np.abs(edge_location_array[1])
                

            
            return [edge_location,edge_location_sig]
            
        self.edge_finder=edge_finder

    #update class runs through all the processing steps and saves all
    #intermediate data from the most up to date parameters
    def update(self,galaxies='des'):
        
        #create list of objects to run through all galaxies
        if galaxies=='des':
            
        
            
            all_de_galaxy_datasets=self.load_all()
        
        elif galaxies=='dsphs':
            
            sets=self.load_sph()
       
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
                
                
            galaxy_object.forecut()
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
        
            
    def create_cls_crossed_data(self):
        
        sets=self.load()
        
        for i in sets:
            i.save_to_csv('crossmatching/ukirt_pre/' + i.galaxy)
            crossmatch(i.galaxy)
            i.gaia_remove('crossmatching/gaia/' + i.galaxy)
            i.save_to_parquet('processed_data/cls_crossed_data/' + i.galaxy)
            
    def update_dsphs(self):
        

            
        sets=self.load_sph()
       
        #carry out processing for each galaxy, saving to binary file after
        #each step
        for i in sets:
            i.save_to_parquet('processed_data/cls_cut_data/' + i.galaxy)
            i.forecut()
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
            
    def edge_update_dsphs(self):
        
        sets=self.load_sph()
        edge=self.edge
        for i in sets:
            i.save_to_parquet('processed_data/cls_cut_data/' + i.galaxy)
            forebord=edge(stage='cls_cut',galaxy=i.galaxy,axis='none')
            #print(i.galaxy)
            #print(forebord[0])
            #print(forebord[1])
            i.forecut(cut=forebord[0])
            i.save_to_parquet('processed_data/fore_cut_data/' + i.galaxy)

        
    

            
    

    
    