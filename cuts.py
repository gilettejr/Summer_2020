from data_read import data_read
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np

class cuts:
    
    def __init__(self,galaxy='ngc147'):
    
        cis=data_read(stage='cls_cut',galaxy=galaxy)
        fore=data_read(stage='fore_cut',galaxy=galaxy)
        agb=data_read(stage='agb',galaxy=galaxy)
        cm=data_read(stage='cm',galaxy=galaxy)
        
        self.cisdata=cis.data
        self.foredata=fore.data
        self.agbdata=agb.data
        self.cdata=cm.cdata
        self.mdata=cm.mdata
        
    def plot_kj_cmd(self):
        
        cisdata=self.cisdata
        foredata=self.foredata
        agbdata=self.agbdata
        cdata=self.cdata
        mdata=self.mdata
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'18','axes.labelsize':'18',
        'axes.titlesize':'18','xtick.labelsize':'18',
        'ytick.labelsize':'18','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        
        cisplot=sns.scatterplot(cisdata.jmag-cisdata.kmag,cisdata.kmag,marker='o',color='black',s=5,linewidth=0,label='Foreground Stars')
        forplot=sns.scatterplot(foredata.jmag-foredata.kmag,foredata.kmag,marker='o',color='orange',s=5,linewidth=0,label='RGB Stars')

        cplot=sns.scatterplot(cdata.jmag-cdata.kmag,cdata.kmag,marker='o',color='red',s=5,linewidth=0,label='Candidate C-Stars')
        mplot=sns.scatterplot(mdata.jmag-mdata.kmag,mdata.kmag,marker='o',color='blue',s=5,linewidth=0,label='Candidate M-Stars')
        
        plt.legend(markerscale=5,frameon=False)
        cisplot.invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')
        
    def plot_cc(self):
        
        cisdata=self.cisdata
        foredata=self.foredata
        agbdata=self.agbdata
        cdata=self.cdata
        mdata=self.mdata
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'18','axes.labelsize':'18',
        'axes.titlesize':'18','xtick.labelsize':'18',
        'ytick.labelsize':'18','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        
        cisplot=sns.scatterplot(cisdata.jmag-cisdata.hmag,cisdata.hmag-cisdata.kmag,marker='o',color='black',s=5,linewidth=0,label='Foreground Stars')
        forplot=sns.scatterplot(foredata.jmag-foredata.hmag,foredata.hmag-foredata.kmag,marker='o',color='orange',s=5,linewidth=0,label='RGB Stars')

        cplot=sns.scatterplot(cdata.jmag-cdata.hmag,cdata.hmag-cdata.kmag,marker='o',color='red',s=5,linewidth=0,label='Candidate C-Stars')
        mplot=sns.scatterplot(mdata.jmag-mdata.hmag,mdata.hmag-mdata.kmag,marker='o',color='blue',s=5,linewidth=0,label='Candidate M-Stars')
        
        plt.legend(markerscale=5,frameon=False,loc='upper left')
        plt.ylabel('$(H-K)_0$')
        plt.xlabel('$(J-H)_0$')

    
    