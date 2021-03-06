from pylab import *
from numpy import *
from random import *
from scipy import *
from matplotlib import *
rcdefaults()
matplotlib.rc('font',family='Bitstream Vera Serif')



def plotsubhess(color,mag,ax,label,levels='None', xlims=[-1,4.5], ylims=[21,11], colormap=cm.plasma, cbarr='None', cbarrtitle='', xlab='$J_0$-$K_0$', ylab='$K_0$', saveas='HessCMD.png'):
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
    
    rc('axes',labelsize=25)
    ### Plot all stars in CMD as small points
    ax.plot(color,mag,linestyle='none',marker ='.',ms=1,color='black',zorder=1,label=label)



    ### Use color and magnitude data to split and average into a binned 2D histogram
    Z, xedges, yedges = np.histogram2d(color,mag,bins=(300,300))
    
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
    
    #plt.gca().invert_yaxis()
    ### Set final options and save figure to default (current folder) or to a user-defined location on disk
    #savefig(saveas,dpi=300)

    
    return