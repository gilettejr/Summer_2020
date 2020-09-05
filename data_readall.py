import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from data_load import data_load

from data_read import data_read

from subhess import plotsubhess

from selection_utils import selection_utils

class data_readall:
    
    def __init__(self,stage='agb',arcmin_radius=False):
        
        if stage!='cls_crossed':
        
            self.n147=data_read(stage=stage,galaxy='ngc147')
            self.n185=data_read(stage=stage,galaxy='ngc185')
            self.n205=data_read(stage=stage,galaxy='ngc205')
            self.m32=data_read(stage=stage,galaxy='m32')
            self.dEs=np.array([self.n147,self.n185,self.n205,self.m32])
            
            self.n147_cm=data_read(stage='cm',galaxy='ngc147')
            self.n185_cm=data_read(stage='cm',galaxy='ngc185')
            self.n205_cm=data_read(stage='cm',galaxy='ngc205')
            self.m32_cm=data_read(stage='cm',galaxy='m32')
        if stage=='agb':
            
            self.n147_crossed=data_read(stage='agb_crossed',galaxy='ngc147')
            
        else:
            
            self.n147_crossed=data_read(stage='cls_crossed',galaxy='ngc147')
            
            
            
        if arcmin_radius!=False:
            e=selection_utils()
            for i in self.dEs:
                i.data = e.select_ellipse(i.data,afl=arcmin_radius/60)[0]
            
        if stage!='cm':
            and1=data_read(stage,'and1')
            and2=data_read(stage,'and2')
            and3=data_read(stage,'and3')
            and6=data_read(stage,'and6')
            and7=data_read(stage,'and7')
            and10=data_read(stage,'and10',)
            and14=data_read(stage,'and14',)
            and15=data_read(stage,'and15')
            and16=data_read(stage,'and16')
            and17=data_read(stage,'and17')
            and18=data_read(stage,'and18')
            and19=data_read(stage,'and19')
            and20=data_read(stage,'and20')
        
            self.sphs=np.array([and1,and2,and3,and6,and7,and10,and14,and15,and16,and17,and18,and19,and20])

        
    def plot_kj_cmds(self,marker='.',markersize=3,color='black'):
        
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
        
        
        #plt.gca().invert_yaxis()
        axs[1,0].set_xlabel('(J-K)$_0$')
        axs[1,1].set_xlabel('(J-K)$_0$')
        axs[1,0].set_ylabel('K$_0$')
        axs[0,0].set_ylabel('K$_0$')
        
        plt.savefig('report_images/cls_all.pdf')
        
        #axs[0,0].set_title('NGC147')
        #axs[0,1].set_title('NGC185')
        #axs[1,0].set_title('NGC205')
        #axs[1,1].set_title('M32')
    
    def plot_kj_cmd_gaia_removal(self):
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        n147_uncleaned=self.n147.data
        
        n147_cleaned=self.n147_crossed.data
        
        fig,axes=plt.subplots(1,2,sharex=True,sharey=True)
        
        sns.scatterplot(n147_uncleaned.jmag-n147_uncleaned.kmag,n147_uncleaned.kmag,ax=axes[0],linewidth=0,s=5,marker='.',color='black',label='Uncleaned')
        sns.scatterplot(n147_cleaned.jmag-n147_cleaned.kmag,n147_cleaned.kmag,ax=axes[1],marker='.',linewidth=0,s=5,color='black',label='Cleaned')
        

        
        for i in axes:
            #i.label_outer()
            leg = i.legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
        for i in axes:
            
            i.invert_yaxis()
            
        plt.subplots_adjust(wspace=0, hspace=0)
        
        axes[0].set_ylabel('K$_0$')
        axes[0].set_xlabel('(J-K)$_0$')
        
        axes[1].set_xlabel('(J-K)$_0$')
        
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
    
    def plot_kj_cmd_hess(self):
        
        n147_uncleaned=self.n147.data
        
        n147_cleaned=self.n147_crossed.data
        
        def jk_errors(data,binsize,start_loc):
            

            bin_locs=[]
            binsy=[]
            binsx=[]

            
            
            while start_loc > np.min(data.kmag.dropna()):
                biny=[]
                binx=[]
                for i in data.index:
                
                    if start_loc > data.kmag[i] > start_loc - binsize:
                        biny.append(data.kerr[i])
                        binx.append(np.sqrt(data.kerr[i]**2 + data.jerr[i]**2))
                biny=np.array(biny)
                binx=np.array(binx)
                avgbiny=np.average(biny)
                avgbinx=np.average(binx)
                binsy.append(avgbiny)
                binsx.append(avgbinx)
                bin_locs.append(start_loc-binsize/2)
                start_loc=start_loc-binsize
                
            return bin_locs,binsx,binsy
            
            
            
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        n147_c=self.n147_cm.cdata
        n185_c=self.n185_cm.cdata
        n205_c=self.n205_cm.cdata
        m32_c=self.m32_cm.cdata
        
        n147_m=self.n147_cm.mdata
        n185_m=self.n185_cm.mdata
        n205_m=self.n205_cm.mdata
        m32_m=self.m32_cm.mdata
        
        dim_data=[[n147,n185],[n205,m32]]
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        markersize=0.5
        marker='o'
        
        trgb_locs=[18.1,17.86,17.93,17.8]
        trgb_lims=(0,4.5)
        fore_locs=[0.992,0.964,1.00,0.92]
    
        fore_lims=(19,11)
        
        trgb_locs_dim=[[trgb_locs[0],trgb_locs[1]],[trgb_locs[2],trgb_locs[3]]]
        fore_locs_dim=[[fore_locs[0],fore_locs[1]],[fore_locs[2],fore_locs[3]]]
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True,figsize=[6,6])
        
        plotsubhess(n147.jmag-n147.kmag,n147.kmag,ax=axs[0,0],label='NGC 147')
        plotsubhess(n185.jmag-n185.kmag,n185.kmag,ax=axs[0,1],label='NGC 185')
        plotsubhess(n205.jmag-n205.kmag,n205.kmag,ax=axs[1,0],label='NGC 205')
        plotsubhess(m32.jmag-m32.kmag,m32.kmag,ax=axs[1,1],label='M32')
        start_loc=19
        binsize=1
        err_xpos=-0.5
        for i in range(2):
            for j in range(2):
                
                xtrgb=trgb_lims
                ytrgb=(trgb_locs_dim[i][j],trgb_locs_dim[i][j])
                
                xfore=(fore_locs_dim[i][j],fore_locs_dim[i][j])
                yfore=fore_lims
        
                err_ypos,xerr,yerr=jk_errors(dim_data[i][j],binsize=binsize,start_loc=start_loc)
                
                err_xpos=np.full(shape=len(err_ypos),fill_value=-0.5)
                
                axs[i,j].errorbar(err_xpos,err_ypos,xerr=xerr,yerr=yerr,linestyle='none',elinewidth=1,capsize=0,color='black')
                axs[i,j].plot(xtrgb,ytrgb,linestyle='--',color='black')
                axs[i,j].plot(xfore,yfore,linestyle='dashdot',color='black')
            
        
        #axs[0,0].plot(n147_m.jmag-n147_m.kmag,n147_m.kmag,color='blue',markersize=markersize,linestyle='none',marker='o',label='C-type')
        #axs[0,1].plot(n185_m.jmag-n185_m.kmag,n185_m.kmag,color='blue',markersize=markersize,linestyle='none',marker='o')
        #axs[1,0].plot(n205_m.jmag-n205_m.kmag,n205_m.kmag,color='blue',markersize=markersize,linestyle='none',marker='o')
        #axs[1,1].plot(m32_m.jmag-m32_m.kmag,m32_m.kmag,color='blue',markersize=markersize,linestyle='none',marker='o')
  
        #axs[0,0].plot(n147_c.jmag-n147_c.kmag,n147_c.kmag,color='red',markersize=markersize,linestyle='none',marker='o')
        #axs[0,1].plot(n185_c.jmag-n185_c.kmag,n185_c.kmag,color='red',markersize=markersize,linestyle='none',marker='o')
        #axs[1,0].plot(n205_c.jmag-n205_c.kmag,n205_c.kmag,color='red',markersize=markersize,linestyle='none',marker='o')
        #axs[1,1].plot(m32_c.jmag-m32_c.kmag,m32_c.kmag,color='red',markersize=markersize,linestyle='none',marker='o')
        

        
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

                #axs[i,j].invert_yaxis()
                axs[i,j].set_ylim(bottom=19,top=13)
                #axs[i,j].set_xlim(left=-1)
                axs[i,j].xaxis.set_major_locator(plt.MaxNLocator(5))
                axs[i,j].yaxis.set_major_locator(plt.MaxNLocator(5))
                axs[i,j].xaxis.set_minor_locator(MultipleLocator(0.5))
                axs[i,j].yaxis.set_minor_locator(MultipleLocator(0.5))
        axs[1,0].set_xlabel('(J-K)$_0$')
        axs[1,1].set_xlabel('(J-K)$_0$')
        axs[1,0].set_ylabel('K$_0$')
        axs[0,0].set_ylabel('K$_0$')
        
        plt.savefig('report_images/cls_all.pdf')
        
        #axs[0,0].set_title('NGC147')
        #axs[0,1].set_title('NGC185')
        #axs[1,0].set_title('NGC205')
        #axs[1,1].set_title('M32')
        
    def plot_cc_hess(self):
        
        n147=self.n147.data
        n185=self.n185.data
        n205=self.n205.data
        m32=self.m32.data
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True,figsize=[6,6])
        
        hkcut=[[0.337,0.323],[0.407,0.477]]
        hklim=2.3
        jhcut=[[0.883,0.857],[0.930,0.913]]
        jhlim=3.18
        
        plotsubhess(n147.jmag-n147.hmag,n147.hmag-n147.kmag,ax=axs[0,0],label='NGC 147')
        plotsubhess(n185.jmag-n185.hmag,n185.hmag-n185.kmag,ax=axs[0,1],label='NGC 185')
        plotsubhess(n205.jmag-n205.hmag,n205.hmag-n205.kmag,ax=axs[1,0],label='NGC 205')
        plotsubhess(m32.jmag-m32.hmag,m32.hmag-m32.kmag,ax=axs[1,1],label='M32')
        

        
        plt.subplots_adjust(wspace=0, hspace=0)
        

        for i in range(2):
            
            for j in range(2):

                #i.label_outer()
                leg = axs[i,j].legend(handlelength=0, handletextpad=0, frameon=False,loc='upper left',markerscale=0.001)
                for item in leg.legendHandles:
                    item.set_visible(False)
                #axs[i,j].xaxis.grid(True, which='minor')
                
                xhk=(jhcut[i][j],jhlim)
                yhk=(hkcut[i][j],hkcut[i][j])
                
                xjh=(jhcut[i][j],jhcut[i][j])
                yjh=(hkcut[i][j],hklim)
                
                axs[i,j].plot(xhk,yhk,linestyle='--',color='black')
                axs[i,j].plot(xjh,yjh,linestyle='--',color='black')
                
                axs[i,j].set_xlim(-0.93,jhlim)
                axs[i,j].set_ylim(0,hklim)
                
                axs[i,j].xaxis.set_minor_locator(MultipleLocator(0.2))
                axs[i,j].yaxis.set_minor_locator(MultipleLocator(0.2))
   