from data_processor import data_processor
from astropy.modeling import fitting
from astropy.modeling.models import Sersic1D

import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

import shapely.affinity

import os




class close_data_processor(data_processor):
    
    def get_close_FEH_slices(self):
        
        if self.galaxy=='ngc205':
            
            cut=5
            
        elif self.galaxy=='m32':
            
            cut=3
        
        c_slice_dataframe=pd.read_parquet('unfit_background_corrected_profiles/'+ self.galaxy + 'c')
        
        m_slice_dataframe=pd.read_parquet('unfit_background_corrected_profiles/'+self.galaxy +'m')
        
        full_slice_nums=c_slice_dataframe.slice_nums + m_slice_dataframe.slice_nums
        
        total_star_num=np.sum(full_slice_nums)

        slice_a=c_slice_dataframe.a
        
        c_slice_dataframe
        
        cm_slices=c_slice_dataframe.density/m_slice_dataframe.density
        print(m_slice_dataframe.density)
        print(c_slice_dataframe.density)



        cm_slices_uncs=cm_slices*np.sqrt((c_slice_dataframe.density_err/c_slice_dataframe.density)**2 + (m_slice_dataframe.density_err/m_slice_dataframe.density)**2)

        FEH_slices=self.CM_to_FEH(cm_slices)
        #print(cm_slices)
        #print(cm_slices_uncs)
        FEH_slice_uncs=self.CM_to_FEH_uncs(cm_slices,cm_slices_uncs)

        
        
        plt.errorbar(slice_a,FEH_slices,marker='o',linestyle='none',color='black',capsize=2,yerr=FEH_slice_uncs)
        plt.ylabel('[Fe/H]')
        plt.xlabel('Semi-major axis (arcmins)')
        
        weights=full_slice_nums[cut:]/total_star_num

        
        
        
        

        
        new_cm_slices=[]
        new_cm_slices_uncs=[]
        new_weights=[]
        for i in range(len(cm_slices)):
            
            if m_slice_dataframe.density[i]<0 or c_slice_dataframe.density[i]<0:
                
                continue
            
            else:
                
                new_cm_slices.append(cm_slices[i])
                new_cm_slices_uncs.append(cm_slices_uncs[i])
                new_weights.append(weights[i])

        avg_CM=np.average(new_cm_slices,weights=new_weights)        
        avg_CM_unc=np.sqrt(np.sum(np.square(new_cm_slices_uncs)))/len(new_cm_slices_uncs)
        
        FEH_avg=self.CM_to_FEH(avg_CM)
        
        FEH_avg_unc=self.CM_to_FEH_uncs(avg_CM,avg_CM_unc)
        print(avg_CM)
        print(avg_CM_unc)
        print(FEH_avg)
        print(FEH_avg_unc)
        
        self.FEH_x=slice_a
        self.FEH_y=FEH_slices
        self.FEH_yunc=FEH_slice_uncs
        
        

        #plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',capsize=2,marker='o',color='black')
        
    def fit_close_slice_count_profile(self,stars='agb',crowding_num=1):

        distribution=pd.read_parquet('unfit_background_corrected_profiles' + self.galaxy + stars)
        
        xdata=distribution.a
        ydata=distribution.density
        yerr=distribution.density_err
        
        
        plt.errorbar(xdata,ydata,yerr=yerr,marker='o',linestyle='none',color='black')
        
        plt.ylabel('Density (N/arcmins$^2$)')
        plt.xlabel('Semi-major axis (arcmins)')
        
        xdata_orig=xdata
            
        xdata=xdata[:len(xdata)-(crowding_num)]
        ydata=ydata[:len(ydata)-(crowding_num)]
        yerr=yerr[:len(yerr)-(crowding_num)]
        
        
        
        if self.galaxy=='ngc205':
        
            model=Sersic1D(r_eff=2.0,n=1.5)
            
        elif self.galaxy=='m32':
            
            model=Sersic1D(r_eff=1.7,n=1.00)
            
            
        
        fit=fitting.LevMarLSQFitter()
        
        s=fit(model,xdata,ydata,weights=1/yerr)

        xdata=np.linspace(min(xdata_orig),max(xdata_orig),num=1000)
                 
        
                    
        plt.plot(xdata,s(xdata),label='Sersic Fit with a$_{eff}$ = ' + str(round(s.r_eff.value,3)) + ', n = ' + str(round(s.n.value,3)))
        plt.legend()