import matplotlib.pyplot as plt
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


class background_utilities(data_processor):
    
   def find_background_density_border(self,stars='agb',marker='o',markersize='1',color='black',borderwidth=0.07,show_figure=False):
        
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
    
   def find_background_grad(self,ellipse_a=0.15,ellipticity=0.43,clockrot=169.2,marker='o',markersize='1',color='black',show_figure=False,binwidth=0.02):
        
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
        
        stars=self.stars
        
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

        
        if show_figure==True:
        
            plt.plot(data.alpha,data.beta,marker=marker,markersize=markersize,linestyle='none',color=color)
            plt.gca().invert_xaxis()
            plt.gca().set_ylabel(r'$\beta$')
            plt.gca().set_xlabel(r'$\alpha$')
        
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
        #fit with Sersic profile
        

#background=pd.read_parquet('205_c_tests')


#stars='agb

   def fit_close_background(self,marker='o',markersize='3',color='black'):

        background =self.background.copy()
        
        
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
        
        plt.errorbar(xdata,ydata,yerr=yerr,capsize=2,marker=marker,color=color)
        
        plt.plot(xdata,s(xdata))
        
        background['density_fit']=s(xdata)
        
        try:
            pd.save_to_parquet('backgrounds/' + self.galaxy + self.stars)
            
        except:
            
            os.system('mkdir backgrounds')
            pd.save_to_parquet('backgrounds/' + self.galaxy + self.stars)
        
        self.background=background
        
        print(background)