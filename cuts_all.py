from data_readall import data_readall

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import seaborn as sns

import numpy as np
#class for plotting populations in different colours to demonstrate selection cuts
class cuts_all:
    
    #initialise and set galaxy
    def __init__(self):
        #read in different selections of data by initialising data_read instances
        cis=data_readall(stage='cls_cut')
        fore=data_readall(stage='fore_cut')
        agb=data_readall(stage='agb')
        cm=data_readall(stage='cm')
        
        #set datasets as attributes
        self.cis=cis.dEs
        self.fore=fore.dEs
        self.agb=agb.dEs
        self.cm=cm.dEs

        
        #set figure size

    
    #plot and colour code selections on a cmd
    def plot_kj_cmds(self,marker='.'):
        

        
        #retrieve attributes
        cis=self.cis
        fore=self.fore
        agb=self.agb
        cm=self.cm
        
        #plotting stuff
        sns.set_context('paper')
        
        params={'legend.fontsize':'15','axes.labelsize':'16',
        'axes.titlesize':'18','xtick.labelsize':'14',
        'ytick.labelsize':'14','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        #plot each selection with a different colour and label on j-k CMD
        xsize=2
        ysize=2
        fig,axs=plt.subplots(ysize,xsize,sharex=True,sharey=True,figsize=[6,6])
        y=0
        x=0
        for i in range(len(cis)):
        
            cisplot=sns.scatterplot(cis[i].data.jmag-cis[i].data.kmag,cis[i].data.kmag,marker=marker,color='black',s=5,linewidth=0,ax=axs[y,x])
            forplot=sns.scatterplot(fore[i].data.jmag-fore[i].data.kmag,fore[i].data.kmag,marker=marker,color='orange',s=5,linewidth=0,ax=axs[y,x])
    
            cplot=sns.scatterplot(cm[i].cdata.jmag-cm[i].cdata.kmag,cm[i].cdata.kmag,marker=marker,color='red',s=5,linewidth=0,ax=axs[y,x])
            mplot=sns.scatterplot(cm[i].mdata.jmag-cm[i].mdata.kmag,cm[i].mdata.kmag,marker=marker,color='blue',s=5,linewidth=0,ax=axs[y,x],label=cm[i].galaxy.upper())
            
            if x!=xsize-1:
                x=x+1
            else:
                x=0
                y=y+1
                
        for i in range(2):
            for j in range(2):
                #i.label_outer()
                leg = axs[i,j].legend(handlelength=0, handletextpad=0, frameon=False)
                for item in leg.legendHandles:
                    item.set_visible(False)
                axs[i,j].set_ylim([19,12])
                #axs[i,j].xaxis.grid(True, which='minor')

                axs[i,j].xaxis.set_minor_locator(MultipleLocator(0.5))
                axs[i,j].yaxis.set_minor_locator(MultipleLocator(0.5))
                    
                #axs[i,j].invert_yaxis()
        #create legend
   
        
        
        #plt.legend(markerscale=3,frameon=False)
        #invert axis and set labels
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[1,0].set_xlabel('(J-K)$_0$')
        axs[1,1].set_xlabel('(J-K)$_0$')

        axs[1,0].set_ylabel('K$_0$')
        axs[0,0].set_ylabel('K$_0$')
        #plt.tight_layout()
        ##plt.savefig('report_images/147kj_cuts.pdf')
        
    #plot and colour code selections on a 2CD
    def plot_cc(self):

        a=plt.figure(figsize=[7,6])
        
        #retrieve attributes
        cisdata=self.cisdata
        foredata=self.foredata
        agbdata=self.agbdata
        cdata=self.cdata
        mdata=self.mdata
        
        #plotting stuff
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'15','axes.labelsize':'16',
        'axes.titlesize':'18','xtick.labelsize':'16',
        'ytick.labelsize':'16','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        #plot each selection with a different colour and label on 2CD
        cisplot=sns.scatterplot(cisdata.jmag-cisdata.hmag,cisdata.hmag-cisdata.kmag,marker='o',color='black',s=5,linewidth=0,label='Foreground Stars')
        forplot=sns.scatterplot(foredata.jmag-foredata.hmag,foredata.hmag-foredata.kmag,marker='o',color='orange',s=5,linewidth=0,label='RGB Stars')

        cplot=sns.scatterplot(cdata.jmag-cdata.hmag,cdata.hmag-cdata.kmag,marker='o',color='red',s=5,linewidth=0,label='Candidate C-Stars')
        mplot=sns.scatterplot(mdata.jmag-mdata.hmag,mdata.hmag-mdata.kmag,marker='o',color='blue',s=5,linewidth=0,label='Candidate M-Stars')
        
        
        #create legend
        plt.legend(markerscale=3,frameon=False,loc='upper left')
        
        #set labels
        plt.ylabel('(H-K)$_0$',labelpad=-8)
        plt.xlabel('(J-H)$_0$')
        #plt.tight_layout()
        plt.savefig('report_images/147cc_cuts.pdf')
    #def plot_cmd_cc(self):
