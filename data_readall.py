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
        
    def plot_kj_cmds(self,marker=',',markersize=1,color='black'):
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True)
        axs[0,0].plot(n147.jmag-n147.kmag,n147.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[0,1].plot(n185.jmag-n185.kmag,n185.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,0].plot(n205.jmag-n205.kmag,n205.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,1].plot(m32.jmag-m32.kmag,m32.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        
        
        for i in range(2):
            
            for j in range(2):
                
                axs[i,j].invert_yaxis()

        axs[1,0].set_xlabel('$J_0$-$K_0$')
        axs[1,1].set_xlabel('$J_0$-$K_0$')
        axs[1,0].set_ylabel('$K_0$')
        axs[0,0].set_ylabel('$K_0$')
        axs[0,0].set_title('NGC147')
        axs[0,1].set_title('NGC185')
        axs[1,0].set_title('NGC205')
        axs[1,1].set_title('M32')
        
    def plot_ccs(self,marker=',',markersize=1,color='black'):
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        fig,axs=plt.subplots(2,2)
        axs[0,0].plot(n147.jmag-n147.hmag,n147.hmag-n147.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[0,1].plot(n185.jmag-n185.hmag,n185.hmag-n185.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,0].plot(n205.jmag-n205.hmag,n205.hmag-n205.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        axs[1,1].plot(m32.jmag-m32.hmag,m32.hmag-m32.kmag,linestyle='none',marker=marker,markersize=markersize,color='black')
        

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
        
        
        
        
        
    

