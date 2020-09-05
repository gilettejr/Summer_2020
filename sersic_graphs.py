from data_read import data_read
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

from matplotlib.patches import Ellipse,Rectangle



class radial_graphs:
    
    def __init__(self,stage='cm',galaxy='ngc147'):
        
        galaxies=['ngc147','ngc185','ngc205','m32']
        
        PAs=[34.2,45.9,169.2,157.9]
        
        ellipticities=[0.46,0.22,0.43,0.14]
        
        tidals=[0.24,0.12,0.2,0.16]
        
        for i in range(len(galaxies)):
            
            if galaxies[i]==galaxy:
                
                self.PA=PAs[i]
                self.ellipticity=ellipticities[i]
                self.galaxy=galaxies[i]
                self.tidal=tidals[i]
                
                break
            
        
        
        galaxy_object=data_read(stage=stage,galaxy=galaxy)
        
        self.galaxy_object=galaxy_object
        
        def subplot_spatial_slices(galaxy_object,PA,ellipticity,ax,separate_stars=False):
            
            galaxy_data=galaxy_object.data
            galaxy_mdata=galaxy_object.mdata
            galaxy_cdata=galaxy_object.cdata
            
            a_width=galaxy_object.a_width
            outer_rad=galaxy_object.outer_rad
            
            majors=np.linspace(a_width,outer_rad,num=int((outer_rad*1000)/(a_width*1000)))
            minors=majors*(1-ellipticity)
        
            ells=[]
            
            if separate_stars==True:
                edgecolor='black'
                
            else:
                
                edgecolor='red'
        
            for i in range(len(majors)):
                ells.append(Ellipse(xy=[0,0],height=majors[i]*2,width=minors[i]*2,angle=360-PA,facecolor='none',edgecolor=edgecolor,linestyle='--',linewidth=1))
            
            if separate_stars==False:
            
                ax.plot(galaxy_data.xi,galaxy_data.eta,linestyle='none',marker='o',markersize=1,color='black',zorder=1)
            else:
                ax.plot(galaxy_mdata.xi,galaxy_mdata.eta,linestyle='none',marker='o',markersize=2,color='black',zorder=1,label='M-type')
                ax.plot(galaxy_cdata.xi,galaxy_cdata.eta,linestyle='none',marker='v',markersize=4,color='red',zorder=1,label='C-type')
                
                
            #axes[1].plot(0,0,linestyle='none',label='PA = ' + str(PA) + ' e = ' + str(ellipticity))
            for i in ells:
                ax.add_artist(i)
                
            ax.invert_xaxis()
            ax.set_ylabel(r'$\eta$')
            ax.set_xlabel(r'$\xi$')
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            
        self.subplot_spatial_slices=subplot_spatial_slices
    
    def plot_sersic_and_ellipses_far(self):
        
        PA=self.PA
        ellipticity=self.ellipticity
        tidal=self.tidal
        
        galaxy_name=self.galaxy
        
        galaxy_object=self.galaxy_object
        
        galaxy_data=galaxy_object.data.copy()
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        

        

        
        galaxy_object.make_slices(outer_rad=0.4)
        
        background=galaxy_object.find_background_density_border()
        

        
        xdata,ydata,yerr,xdata_fit,ydata_fit=galaxy_object.slice_count_profile(background_deg=background,crowding_num=1)
        



        fig,axes=plt.subplots(1,2)
        self.subplot_spatial_slices(galaxy_object,PA,ellipticity,axes[1])
        axes[0].errorbar(xdata,ydata,yerr=yerr,capsize=2,linestyle='none',color='black',marker='o',markersize='5')
        axes[0].plot(xdata_fit,ydata_fit,label = 'a$_{eff}$ = ' + str(round(galaxy_object.r_eff,3)) + ' n = ' + str(round(galaxy_object.n,3)))

        #axes[0].set_title(galaxy_name.upper())
        #axes[1].set_title(galaxy_name.upper())
        fig.suptitle(galaxy_name.upper(),x=0.51)
        axes[0].set_yscale('log')

        
        axes[0].set_xlabel('Semi-major axis (arcmins)')
        axes[0].set_ylabel('N$_{sources}$/arcmins$^2$')
        
        axes[0].xaxis.set_minor_locator(MultipleLocator(1))

        
        
        for i in axes:
            #i.label_outer()
            leg = i.legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
        
        
    def plot_radial_CM_far(self):
        
        PA=self.PA
        ellipticity=self.ellipticity
        tidal=self.tidal
        
        galaxy_name=self.galaxy
        
        galaxy_object=self.galaxy_object
        
        galaxy_data=galaxy_object.data.copy()
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        

        

        
        galaxy_object.make_slices(outer_rad=tidal)
        
        mbackground=galaxy_object.find_background_density_border(stars='m')
        cbackground=galaxy_object.find_background_density_border(stars='c')
        
        crowding_num=1
        xdata,ydata,yerr,avgy,avgyunc=galaxy_object.FEH_slices(mbackground,cbackground,crowding_num=crowding_num)
        
        fig,axes=plt.subplots(1,2)
        xdata=xdata*60
        self.subplot_spatial_slices(galaxy_object,PA,ellipticity,axes[1],separate_stars=True)
        
        axes[0].errorbar(xdata,ydata,yerr=yerr,capsize=2,color='black',linestyle='none',marker='o',markersize=3)
        axes[0].errorbar(xdata[-crowding_num],ydata[-crowding_num],yerr=yerr[-crowding_num],color='red',linestyle='none',capsize=2,marker='o',markersize=3)
        axes[0].plot([np.min(xdata),np.max(xdata)],[avgy,avgy],linestyle='--',color='black',linewidth=2,label='Mean [Fe/H] = ' + str(round(avgy,2)) + 'dex')
        axes[0].xaxis.set_minor_locator(MultipleLocator(1))
        axes[0].yaxis.set_minor_locator(MultipleLocator(0.05))
        
        axes[0].set_ylabel('[Fe/H]')
        axes[0].set_xlabel('Semi-major axis (arcmins)')
        fig.suptitle(galaxy_name.upper(),x=0.51)       
        

            #i.label_outer()
        leg = axes[0].legend(handletextpad=0, frameon=False)
        for item in leg.legendHandles:
            item.set_visible(False)

        
        plt.legend(markerscale=3,frameon=False)
        
        
        
        