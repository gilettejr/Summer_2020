#imports
from data_read import data_read
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
#class for plotting populations in different colours to demonstrate selection cuts
class cuts:
    
    #initialise and set galaxy
    def __init__(self,galaxy='ngc147'):
        #read in different selections of data by initialising data_read instances
        cis=data_read(stage='cls_cut',galaxy=galaxy)
        fore=data_read(stage='fore_cut',galaxy=galaxy)
        agb=data_read(stage='agb',galaxy=galaxy)
        cm=data_read(stage='cm',galaxy=galaxy)
        
        #set datasets as attributes
        self.cisdata=cis.data
        self.foredata=fore.data
        self.agbdata=agb.data
        self.cdata=cm.cdata
        self.mdata=cm.mdata
        
        #set figure size

    
    #plot and colour code selections on a cmd
    def plot_kj_cmd(self):
        
        a=plt.figure(figsize=[7,6])
        
        #retrieve attributes
        cisdata=self.cisdata
        foredata=self.foredata
        agbdata=self.agbdata
        cdata=self.cdata
        mdata=self.mdata
        
        #plotting stuff
        sns.set_context('paper')
        
        params={'legend.fontsize':'15','axes.labelsize':'18',
        'axes.titlesize':'18','xtick.labelsize':'18',
        'ytick.labelsize':'18','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        #plot each selection with a different colour and label on j-k CMD
        cisplot=sns.scatterplot(cisdata.jmag-cisdata.kmag,cisdata.kmag,marker='o',color='black',s=5,linewidth=0,label='Foreground Stars')
        forplot=sns.scatterplot(foredata.jmag-foredata.kmag,foredata.kmag,marker='o',color='orange',s=5,linewidth=0,label='RGB Stars')

        cplot=sns.scatterplot(cdata.jmag-cdata.kmag,cdata.kmag,marker='o',color='red',s=5,linewidth=0,label='Candidate C-Stars')
        mplot=sns.scatterplot(mdata.jmag-mdata.kmag,mdata.kmag,marker='o',color='blue',s=5,linewidth=0,label='Candidate M-Stars')
        
        
        #create legend
        plt.legend(markerscale=3,frameon=False)
        #invert axis and set labels
        cisplot.invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
    #plot and colour code selections on a 2CD
    def plot_cc(self):

        a=plt.figure(figsize=[12,12])
        
        #retrieve attributes
        cisdata=self.cisdata
        foredata=self.foredata
        agbdata=self.agbdata
        cdata=self.cdata
        mdata=self.mdata
        
        #plotting stuff
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'15','axes.labelsize':'18',
        'axes.titlesize':'18','xtick.labelsize':'18',
        'ytick.labelsize':'18','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
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
        plt.ylabel('$(H-K)_0$')
        plt.xlabel('$(J-H)_0$')

    
    