from data_processor import data_processor
from astropy.modeling import fitting
from astropy.modeling.models import Sersic1D

import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

import shapely.affinity




class close_data_processor(data_processor):
    
    def close_FEH_slices(self):
        
        if self.galaxy=='ngc205':
            
            cut=5
            
        elif self.galaxy=='m32':
            
            cut=3
        
        c_slice_dataframe=pd.read_parquet(self.galaxy + '_radial_density_c')
        
        m_slice_dataframe=pd.read_parquet(self.galaxy +'_radial_density_m')
        
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
        
    def close_slice_count_profile(self):
        
        print('Make sure youve got a background dataframe for subtraction!')
        #remake alpha coordinates
        slice_shapes=self.slice_shapes
        background=self.background
        rotate_coords=self.rotate_coords
        stars=self.stars
        if stars=='agb':
            
            data=self.data
        elif stars=='m':
            
            
            data=self.mdata
            
        elif stars=='c':
            
            data=self.cdata
            
        areas=self.areas
        slices=self.slices
        bin_shapes=self.bin_shapes
        
        outer_rad=self.outer_rad
        a_width=self.a_width
        
        m31ra=10.68470833
        m31dec=41.26875
        
        if self.galaxy=='ngc205':
            tra=10.09189356
            tdec=41.68541564
        elif self.galaxy=='m32':
            tra=10.6742708
            tdec=40.8651694
        
        m31xi,m31eta=self.eq_to_tan(m31ra,m31dec,tra,tdec)
        

        
        #find angle of rotation
        theta=np.arctan((m31eta)/(m31xi))
        
        xi=data.xi.copy()
        eta=data.eta.copy()
        theta_deg=np.degrees(theta)
        #rotate
        
        
        
        for i in slices:
            i['alpha']=rotate_coords(xi,eta,theta)[0]
            i['beta']=rotate_coords(xi,eta,theta)[1]
        slice_shapes_rot=[]
        for i in slice_shapes:
        
            slice_shapes_rot.append(shapely.affinity.rotate(i,-theta_deg))
        
        
        slice_nums=[]
        
        for i in slices:
            
            slice_nums.append(len(i.dropna()))
            
        slice_nums=np.array(slice_nums)
        areas=np.array(areas)
        
        slice_densities_deg=slice_nums/areas
        slice_densities_unc=np.sqrt(slice_nums)
        slice_densities_min=slice_densities_deg/3600
        slice_densities_unc=(slice_densities_unc/areas)/3600
        #background columns: bin_densities, bin_uncs, bin_locs
        #areas in degrees for ease, background in mins
        new_density=[]
        new_uncs=[]
        for i in range(len(slice_shapes_rot)):
            slice_bin_densities=[]
            slice_bin_densities_uncs=[]
            slice_bin_areas=[]
            for j in background.index:
                

                
                if slice_shapes_rot[i].intersects(bin_shapes[j])==True:
                    
                    overlap_shape=slice_shapes_rot[i].intersection(bin_shapes[j])
                    overlap_area=overlap_shape.area
                    
                    
                    #subtract background
                    overlap_density=slice_densities_min[i]-background.density_fit[j]
                    overlap_density_unc=np.sqrt(slice_densities_unc[i]**2 + background.bin_uncs[j]**2)
                    
                    slice_bin_densities.append(overlap_density)
                    slice_bin_densities_uncs.append(overlap_density_unc)
                    slice_bin_areas.append(overlap_area)
                    
            slice_bin_areas=np.array(slice_bin_areas)
            slice_bin_densities=np.array(slice_bin_densities)
            
            weightings=slice_bin_areas/slice_shapes_rot[i].area
            
            avg_density=np.average(slice_bin_densities,weights=weightings)
            density_err=(np.sqrt(np.sum(np.square(slice_bin_densities_uncs))))/len(slice_bin_densities)
            new_density.append(avg_density)
            new_uncs.append(density_err)

            
        xdata=np.linspace(outer_rad-a_width/2,0+a_width/2,num=(outer_rad*1000)/(a_width*1000))
        xdata=xdata*60
            
        ydata=new_density
        yerr=new_uncs
            
        fig=plt.figure()
        

        
        density_distribution=pd.DataFrame({'a':xdata,'density':ydata,'density_err':yerr,'slice_area_deg':areas,'slice_nums':slice_nums})
        

            
        outfilename=self.galaxy+'_radial_density_' + stars
        
        density_distribution.to_parquet(outfilename)
            
        plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',capsize=2,marker='o',color='black')
        
    def fit_close_profiles(self,crowding_num=1):
        stars='agb'
        distribution=pd.read_parquet(self.galaxy + '_radial_density_' + stars)
        
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