from data_read import data_reader

import seaborn as sns
import matplotlib.pyplot as plt

class agb_gaia_cleaning_plotter:
    
    def __init__(self,galaxy):
        
        self.cls_cut_galaxy_object=data_reader(galaxy=galaxy,stage='cls_cut')
        self.cls_crossed_galaxy_object=data_reader(galaxy=galaxy,stage='cls_crossed')
        
        self.agb_cut_galaxy_object=data_reader(galaxy=galaxy,stage='agb')
        self.agb_crossed_galaxy_object=data_reader(galaxy=galaxy,stage='agb_crossed')
        
        
        
    def plot_gaia_removal(self):
        
        n147_uncleaned=self.n147.data
        
        n147_cleaned=self.n147_crossed.data
        
        n147_fore=n147_uncleaned[~n147_uncleaned.isin(n147_cleaned)].dropna()
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})  
        
        fig,axes=plt.subplots(1,2)
        
        
        
        for i in axes:
            #i.label_outer()
            leg = i.legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)

        sns.scatterplot(n147_cleaned.jmag-n147_cleaned.kmag,n147_cleaned.kmag,s=10,marker='.',linewidth=0,color='black',ax=axes[0])
        sns.scatterplot(n147_fore.jmag-n147_fore.kmag,n147_fore.kmag,s=50,marker='o',color='red',linewidth=0,ax=axes[0])
        
        sns.scatterplot(n147_cleaned.jmag-n147_cleaned.hmag,n147_cleaned.hmag-n147_cleaned.kmag,s=10,marker='.',linewidth=0,color='black',ax=axes[1])
        sns.scatterplot(n147_fore.jmag-n147_fore.hmag,n147_fore.hmag-n147_fore.kmag,s=50,marker='o',color='red',linewidth=0,ax=axes[1])
        
        axes[0].invert_yaxis()
            
        #plt.subplots_adjust(wspace=0, hspace=0)
        
        axes[0].set_ylabel('K$_0$')
        axes[0].set_xlabel('(J-K)$_0$')
        
        axes[1].set_xlabel('(J-H)$_0$')
        axes[1].set_ylabel('(H-K)$_0$')