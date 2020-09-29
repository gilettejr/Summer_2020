import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from data_processor import data_processor
from selection_utils import selection_utils
import os


from astropy.modeling import fitting
from astropy.modeling.models import Sersic1D

from shapely.geometry import Point, box, MultiPoint
from shapely.geometry.polygon import Polygon,LinearRing
import shapely.affinity

from matplotlib.patches import Ellipse,Rectangle


class background_constructor(data_processor):
    
   def find_background_density_border(self,stars,marker='o',markersize='1',color='black',borderwidth=0.07,show_figure=False):
        
        if stars=='agb':
            
            data=self.data
            
        elif stars=='m':
            
            data=self.mdata
            
        elif stars =='c':
            
            data=self.cdata
        
        border_data=data.copy()
        corner_locs=self.data_corners(self.data)
        
        for i in corner_locs:
            
            i=np.asarray(i)
            
        
        subtraction_identity=np.array([[-1,-1],[-1,1],[1,1],[1,-1]])
        
        spatial_subtraction_matrix=borderwidth * subtraction_identity
        
        inner_corner_locs=corner_locs + spatial_subtraction_matrix
            
            
        inner_corner_locs=map(tuple,inner_corner_locs)
        inner_corner_locs=tuple(inner_corner_locs)
        
        select_border=selection_utils()
        select_border.select_stars(border_data,inner_corner_locs[0],inner_corner_locs[2],unselect=True)
        
        inner_area=(inner_corner_locs[0][0]-inner_corner_locs[2][0]) * (inner_corner_locs[0][1]-inner_corner_locs[2][1])
        outer_area=(corner_locs[0][0]-corner_locs[2][0]) * (corner_locs[0][1]-corner_locs[2][1])
        
        
        border_area=outer_area-inner_area

        
        
        border_num=len(border_data.dropna())
        
        border_err=np.sqrt(border_num)
        
        
        if show_figure==True:
            
            borders=[Rectangle(xy=inner_corner_locs[2],width=inner_corner_locs[0][0]-inner_corner_locs[2][0],height=inner_corner_locs[0][1]-inner_corner_locs[2][1],fill=False)]
            
            
            plt.rc('axes',labelsize=20)
            fig,ax=plt.subplots()
            ax.plot(data.xi,data.eta,linestyle='none',marker=marker,markersize=markersize,color=color,zorder=1)
            for i in borders:
                
                ax.add_artist(i)
            ax.invert_xaxis()
            ax.set_ylabel(r'$\eta$')
            ax.set_xlabel(r'$\xi$')
            
         
            
        return np.array([border_num/border_area,border_err/border_area,border_num]) 
    
   def find_background_grad(self,stars,ellipse_a=0.15,ellipticity=0.43,clockrot=169.2,marker='o',markersize='1',color='black',show_figure=False,binwidth=0.02):
        
        self.stars=stars
       
        if self.galaxy=='m32':
            
            ellipticity=0.14
            clockrot=157.9
            ellipse_a=0.1
        
        def outer_boundary(data,x,y):
            
            corner1loc=np.where(data==np.min(data.alpha))[0][0]
            corner2loc=np.where(data==np.min(data.beta))[0][0]
            corner3loc=np.where(data==np.max(data.alpha))[0][0]
            corner4loc=np.where(data==np.max(data.beta))[0][0]
            
            corner1loc=data.index[data['alpha'] == np.min(data.alpha)].tolist()[0]
            corner2loc=data.index[data['beta'] == np.min(data.beta)].tolist()[0]
            corner3loc=data.index[data['alpha'] == np.max(data.alpha)].tolist()[0]
            corner4loc=data.index[data['beta'] == np.max(data.beta)].tolist()[0]
            
            
            
            corner1=(data.alpha[corner1loc],data.beta[corner1loc])
            corner2=(data.alpha[corner2loc],data.beta[corner2loc])
            corner3=(data.alpha[corner3loc],data.beta[corner3loc])
            corner4=(data.alpha[corner4loc],data.beta[corner4loc])
            
            boundary=LinearRing([corner1,corner2,corner3,corner4])
            
            return boundary
        
        
        if stars=='agb':
            
            data=self.data
            
        elif stars=='m':
            
            data=self.mdata
            
        elif stars=='c':
            
            data=self.cdata
        

        
        #wipe galaxy
        
        e=selection_utils()
        
        dE_stars=e.select_ellipse(data,afl=ellipse_a,eccentricity=ellipticity,clockrot=clockrot,unselect=True)
        dE_ellipse=dE_stars[2]
        
        data=dE_stars[0]

        
        #convert m31 coords to local tangent coords

        
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
        
        print(theta)
        
        xi=data.xi.copy()
        eta=data.eta.copy()

        #rotate
        data['alpha']=self.rotate_coords(xi,eta,theta)[0]
        data['beta']=self.rotate_coords(xi,eta,theta)[1]
        
        m31xi=np.degrees(m31xi)
        m31eta=np.degrees(m31eta)
        
        m31_alpha=self.rotate_coords(m31xi,m31eta,theta)[0]
        m31_beta=self.rotate_coords(m31xi,m31eta,theta)[1]
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        
        if show_figure==True:
        
            plt.plot(data.alpha,data.beta,marker=marker,markersize=markersize,linestyle='none',color=color,label=self.galaxy.upper())
            if self.galaxy=='m32':
                
                plt.xlim(0.3,-0.3)
                plt.ylim(-0.3,0.3)
                
            else:
                
                plt.gca().invert_xaxis()
            plt.gca().set_ylabel(r'$\beta$ (degrees)')
            plt.gca().set_xlabel(r'$\alpha$ (degrees)')
        
        corner1=np.array([np.min(data.xi),np.min(data.eta)])
        corner2=np.array([np.min(data.xi),np.max(data.eta)])
        corner3=np.array([np.max(data.xi),np.max(data.eta)])
        corner4=np.array([np.max(data.xi),np.min(data.eta)])
        
        corners=[corner1,corner2,corner3,corner4]
        corners_rot=[]
        for i in corners:
            
            corners_rot.append(self.rotate_coords(i[0],i[1],theta))
        
        
        


        border=LinearRing(corners_rot)
        s=Polygon(border)
        corners_x=np.array([corners_rot[0][0],corners_rot[1][0],corners_rot[2][0],corners_rot[3][0]])
        
        upper_bound=np.max(corners_x)-binwidth
        lower_bound=np.min(corners_x)
        
        print(upper_bound)
        print(lower_bound)
    
        
        #make bins
        


        #bin1
        
        bin_nums=[]
        bin_areas=[]
        bin_locs=[]
        #define masked ellipse
        mask=dE_ellipse
        theta_deg=np.degrees(theta)
        mask=shapely.affinity.rotate(dE_ellipse,-theta_deg)
        x,y=mask.exterior.xy

        bin_shapes=[]
        while upper_bound-lower_bound > binwidth:
            
        
            bin_slice_fill=box(upper_bound-binwidth,-1.0,upper_bound,1.0)

            bin_slice=bin_slice_fill.exterior
        
            #find intercepts
            bin_corners=border.intersection(bin_slice)
            
            swap0=bin_corners[0]
            swap1=bin_corners[1]
            swap2=bin_corners[3]
            swap3=bin_corners[2]
            
            corners=[swap0,swap1,swap2,swap3]
            bin_corners=MultiPoint([swap0,swap1,swap2,swap3])
            final_bin=Polygon(bin_corners)
            bin_shapes.append(final_bin)
            x,y=final_bin.exterior.xy
            if show_figure==True:
                
                
                
                plt.plot(x,y,color='orange')
                
                leg = plt.legend(handlelength=0, handletextpad=0, frameon=False,loc='upper right',markerscale=0.001)
                for item in leg.legendHandles:
                    item.set_visible(False)
                
                plt.savefig(self.galaxy + '_background_bins.png')
                plt.savefig(self.galaxy + '_background_bins.pdf')

            final_bin_area=final_bin.area
            
            box_data=data.copy()
            
            for i in box_data.index:
                
                if final_bin.contains(Point(box_data.alpha[i],box_data.beta[i]))==False:
                    
                    box_data.loc[i]=np.nan
                    
            bin_nums.append(len(box_data.dropna()))
            
            #subtract area overlapping with ellipse
            
            if bin_slice_fill.intersects(mask)==True:
                
                overlap_area=mask.intersection(bin_slice_fill).area

                final_bin_area=final_bin_area-overlap_area
            
            bin_locs.append(upper_bound-(binwidth/2))
            bin_areas.append(final_bin_area)
            upper_bound=upper_bound-binwidth
            
            if self.galaxy=='m32' and upper_bound-lower_bound < (binwidth * 2):
                
                break
            
        bin_nums=np.array(bin_nums)
        bin_locs=np.array(bin_locs)
        bin_areas=np.array(bin_areas)
        
        
        
        bin_uncs=np.sqrt(bin_nums)
        
        bin_densities=(bin_nums/bin_areas)/3600
        
        bin_uncs=(bin_uncs/bin_areas)/3600
        
        bin_locs=bin_locs * 60
        
        m31_alpha=m31_alpha*60

        m31_bin_locs=bin_locs-m31_alpha

        background=pd.DataFrame({'bin_densities':bin_densities,'bin_uncs':bin_uncs,'bin_locs':bin_locs,'m31_bin_locs':m31_bin_locs})
        
        self.background=background
        
        self.bin_shapes=bin_shapes
        self.stars=stars
        #fit with Sersic profile
        

#background=pd.read_parquet('205_c_tests')




   def fit_close_background(self,marker='o',markersize='3',color='black'):

        background = self.background.copy()
        
        
        xdata=-background.m31_bin_locs
        
        ydata=background.bin_densities
        
        yerr=background.bin_uncs
        
        
        
        model=Sersic1D(r_eff=2,n=10.00)
            
            
        fit=fitting.LevMarLSQFitter()
        
        if self.galaxy=='ngc205' and self.stars=='c':
            
        
            model=Sersic1D(r_eff=10,n=20)
            s=fit(model,xdata[2:],ydata[2:])
            
            
        elif self.galaxy=='m32' and self.stars=='c':
            model=Sersic1D(r_eff=20,n=10)
            s=fit(model,xdata[3:],ydata[3:])
        elif self.galaxy=='m32':
            
            model=Sersic1D(r_eff=20,n=2)
            s=fit(model,xdata[3:],ydata[3:],weights=1/yerr[3:])
        
        else:
            
            s=fit(model,xdata[2:],ydata[2:],weights=1/yerr[2:])
            
        plt.figure()
        
        sat_xdata=background.bin_locs
        
        plt.gca().set_xlabel(r'$\alpha$ (arcmins)')
        plt.gca().set_ylabel('Density (N/arcmin$^2$)')
        plt.gca().invert_xaxis()
        plt.errorbar(sat_xdata,ydata,yerr=yerr,capsize=2,linestyle='none',marker=marker,color=color)
        
        plt.plot(sat_xdata,s(xdata),label=self.galaxy.upper())
        
        leg = plt.legend(handlelength=0, handletextpad=0, frameon=False,loc='upper right',markerscale=0.001)
        for item in leg.legendHandles:
            item.set_visible(False)
        
        plt.savefig(self.galaxy + self.stars + '_background_fit.png')
        plt.savefig(self.galaxy + self.stars + '_background_fit.pdf')
        
        background['density_fit']=s(xdata)
        
        
        self.background=background
        
        print(self.stars)
        
        print(background)
        
   def find_close_slice_profile(self):
        
       
        stars=self.stars
        
        print(stars)
        print('again')
        
        print('Make sure youve got a background dataframe for subtraction!')
        print('You must have run the make_close_background function to create this' )
        #remake alpha coordinates
        slice_shapes=self.slice_shapes
        background=self.background
        rotate_coords=self.rotate_coords

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

            
        xdata=np.linspace(outer_rad-a_width/2,0+a_width/2,num=len(new_density))
        xdata=xdata*60
            
        ydata=new_density
        yerr=new_uncs
            

        

        
        density_distribution=pd.DataFrame({'a':xdata,'density':ydata,'density_err':yerr,'slice_area_deg':areas,'slice_nums':slice_nums})
        

        
        
        
        outfilename=self.galaxy + self.stars
        
        print(outfilename)
        
        try:
            density_distribution.to_parquet('unfit_background_corrected_profiles/' + outfilename)
            
        except:
            
            os.system('mkdir unfit_background_corrected_profiles')
            density_distribution.to_parquet('unfit_background_corrected_profiles/' + outfilename)
            