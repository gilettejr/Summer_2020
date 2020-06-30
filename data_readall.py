import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from data_load import data_load

from data_read import data_read

class data_readall:
    
    def __init__(self,stage='agb'):
        
        self.n147=data_read(stage=stage,galaxy='ngc147')
        self.n185=data_read(stage=stage,galaxy='ngc185')
        self.n205=data_read(stage=stage,galaxy='ngc205')
        self.m32=data_read(stage=stage,galaxy='m32')
        
    def plot_kj_cmds(self,marker='o',markersize=3,color='black'):
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'18',
        'axes.titlesize':'14','xtick.labelsize':'12',
        'ytick.labelsize':'12','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True,figsize=[6,6])
        
        sns.scatterplot(n147.jmag-n147.kmag,n147.kmag,marker=marker,s=markersize,linewidth=0,color='black',ax=axs[0,0],label='NGC 147')
        sns.scatterplot(n185.jmag-n185.kmag,n185.kmag,marker=marker,s=markersize,linewidth=0,color='black',ax=axs[0,1],label='NGC 185')
        sns.scatterplot(n205.jmag-n205.kmag,n205.kmag,marker=marker,s=markersize,linewidth=0,color='black',ax=axs[1,0],label='NGC 205')
        sns.scatterplot(m32.jmag-m32.kmag,m32.kmag,marker=marker,s=markersize,linewidth=0,color='black',ax=axs[1,1],label='M32')
        

        
        plt.subplots_adjust(wspace=0, hspace=0)
        
        axes=[axs[0,0],axs[0,1],axs[1,0],axs[1,1]]
        for i in axes:
            #i.label_outer()
            leg = i.legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
        #axs[0,0].plot(n147.jmag-n147.kmag,n147.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        #axs[0,1].plot(n185.jmag-n185.kmag,n185.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        #axs[1,0].plot(n205.jmag-n205.kmag,n205.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        #axs[1,1].plot(m32.jmag-m32.kmag,m32.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        
        #axs[0,0].set_ylim(11,20)
        #axs[0,0].set_xlim(-0.6,4.5)
        
        
        for i in range(2):
            
            for j in range(2):
                
                axs[i,j].invert_yaxis()

        axs[1,0].set_xlabel('(J-K)$_0$')
        axs[1,1].set_xlabel('(J-K)$_0$')
        axs[1,0].set_ylabel('K$_0$')
        axs[0,0].set_ylabel('K$_0$')
        
        plt.savefig('report_images/cls_all.pdf')
        
        #axs[0,0].set_title('NGC147')
        #axs[0,1].set_title('NGC185')
        #axs[1,0].set_title('NGC205')
        #axs[1,1].set_title('M32')
        
    def plot_ccs(self,marker='o',markersize=1,color='black'):
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True)
        axs[0,0].plot(n147.jmag-n147.hmag,n147.hmag-n147.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[0,1].plot(n185.jmag-n185.hmag,n185.hmag-n185.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,0].plot(n205.jmag-n205.hmag,n205.hmag-n205.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,1].plot(m32.jmag-m32.hmag,m32.hmag-m32.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        
        axs[0,0].set_ylim(-0.6,3.3)
        axs[0,0].set_xlim(-0.3,3.3)

        axs[1,0].set_xlabel('$J_0$-$H_0$')
        axs[1,1].set_xlabel('$J_0$-$H_0$')
        axs[1,0].set_ylabel('$H_0$-$K_0$')
        axs[0,0].set_ylabel('$H_0$-$K_0$')
        axs[0,0].set_title('NGC147')
        axs[0,1].set_title('NGC185')
        axs[1,0].set_title('NGC205')
        axs[1,1].set_title('M32')
        
    def plot_spatials(self,marker='o',markersize=1,color='black',sharex=True,sharey=True):
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        fig,axs=plt.subplots(2,2)
        axs[0,0].plot(n147.xi,n147.eta,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[0,1].plot(n185.xi,n185.eta,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,0].plot(n205.xi,n205.eta,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,1].plot(m32.xi,m32.eta,linestyle='none',marker=marker,markersize=markersize,color='black')
        
        axs[0,0].invert_xaxis()
        axs[0,1].invert_xaxis()
        axs[1,0].invert_xaxis()
        axs[1,1].invert_xaxis()
        
        axs[1,0].set_xlabel(r'$\xi$')
        axs[1,1].set_xlabel(r'$\xi$')
        axs[1,0].set_ylabel(r'$\eta$')
        axs[0,0].set_ylabel(r'$\eta$')
        axs[0,0].set_title('NGC147')
        axs[0,1].set_title('NGC185')
        axs[1,0].set_title('NGC205')
        axs[1,1].set_title('M32')
        
        
        
        
        
    

