
import numpy as np
from random import *
from scipy import *
from matplotlib import *
rcdefaults()

def create_hess_diagram(xdata,ydata,ax,label,xlims=[-1,4.5], ylims=[21,11],levels='None', colormap=cm.plasma,cbarr='None', cbarrtitle='', xlab='$J_0$-$K_0$', ylab='$K_0$'):
    rc('axes',labelsize=25)
    ### Plot all stars in CMD as small points
    ax.plot(xdata,ydata,linestyle='none',marker ='.',ms=1,color='black',zorder=1,label=label)



    ### Use color and magnitude data to split and average into a binned 2D histogram
    Z, xedges, yedges = np.histogram2d(xdata,ydata,bins=(300,300))
    
    ### Find the cell centers from the cell edges
    x = 0.5*(xedges[:-1] + xedges[1:])
    y = 0.5*(yedges[:-1] + yedges[1:])
    
    ### Promote to 2D arrays
    Y, X = np.meshgrid(y, x)
    
    ### Mask out values without data
    Z = np.ma.array(Z, mask=(Z==0))
    
    ### Set levels options as default or user-defined option
    if levels=='None':
        cntr=ax.contourf(X,Y,Z,cmap=colormap,zorder=2)
    else:
        cntr=ax.contourf(X,Y,Z,levels=levels,cmap=colormap,zorder=2)
    cntr.cmap.set_under('white')
    
    ### Set plot limits as default or user-defined option
    #xlim(xlims[0],xlims[1]);ylim(ylims[0],ylims[1])
    
    ### Determine whether to plot a color bar
    if cbarr!='None':
        cbar_ax = fig.add_axes([0.91, 0.14, 0.02, 0.7])
        d=fig.colorbar(cntr, cax=cbar_ax)#, title=cbarrtitle)
        d.set_label(label=cbarrtitle,size=10)

#matplotlib.rc('font',family='Bitstream Vera Serif')

class subhess:

    def plot_kj_cmd(galaxy_object,ax,label):
        """
        plotHess(mag,color) is a function intended to create a contour plot of stellar density in the canonical Color-Magnitude Diagram (CMD). plotHess() plots all the stars in the CMD and then overplots a contour of their density.
        mag = magnitudes of the stars
        color = colors of the stars
        mag and color are required to run plotHess() and need to be the same length arrays.
        Optional keyword arguments:
        
        =========   =======================================================
        Keyword     Description
        =========   =======================================================
        levels:    levels to be used for density contour; if 'None', the defaults are defined by the contour() function, which may not be optimal for the users dataset and may be changed if necessary.
        ylims:      define the y-range of the plotted values; defaults will need to be changed to fit user's data.
        xlims:      define the x-range of the plotted values; defaults will need to be changed to fit user's data.
        colormap:   allows the user to choose a Python color map to use for the contour; the default is grayscale.
        cbarr:      'None' means the colorbar will not be plot as a default. To change this, set cbarr='Yes'; in a true Hess diagram, the cbarr values represent the stellar point density in the CMD plot.
        cbarrtitle: sets the title for the colorbar; should be string
        xlabel:     sets the xlabel for the plot; should be string
        ylabel:     sets the ylabel for the plot; should be string
        saveas:     pathway for saving the output plot. The default is to save in the same folder as "HessCMD.png"
        
        
        """
        ### Set figure options
        
        galaxy_data=galaxy_object.data
        
        xdata=galaxy_data.jmag-galaxy_data.kmag
        ydata=galaxy_data.kmag
        
        create_hess_diagram(xdata,ydata,ax,label)
        
    
    def plot_cc(galaxy_object,ax,label):
        
        galaxy_data=galaxy_object.data
        
        xdata=galaxy_data.jmag-galaxy_data.hmag
        ydata=galaxy_data.hmag-galaxy_data.kmag
        
        create_hess_diagram(xdata,ydata,ax,label)
        
