import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from data_load import data_load

from data_read import data_read

from sub_hess import plotsubhess

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
        
        
        plt.gca().invert_yaxis()
        axs[1,0].set_xlabel('(J-K)$_0$')
        axs[1,1].set_xlabel('(J-K)$_0$')
        axs[1,0].set_ylabel('K$_0$')
        axs[0,0].set_ylabel('K$_0$')
        
        plt.savefig('report_images/cls_all.pdf')
        
        #axs[0,0].set_title('NGC147')
        #axs[0,1].set_title('NGC185')
        #axs[1,0].set_title('NGC205')
        #axs[1,1].set_title('M32')
        
    def plot_kj_cmd_hess(self):
        
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
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True,figsize=[6,6])
        
        plotsubhess(n147.jmag-n147.kmag,n147.kmag,ax=axs[0,0],label='NGC 147')
        plotsubhess(n185.jmag-n185.kmag,n185.kmag,ax=axs[0,1],label='NGC 185')
        plotsubhess(n205.jmag-n205.kmag,n205.kmag,ax=axs[1,0],label='NGC 205')
        plotsubhess(m32.jmag-m32.kmag,m32.kmag,ax=axs[1,1],label='M32')
        start_loc=20
        binsize=1
        err_xpos=-1.0
        for i in range(2):
            for j in range(2):
        
                err_ypos,xerr,yerr=jk_errors(dim_data[i][j],binsize=binsize,start_loc=start_loc)
                
                err_xpos=np.full(shape=len(err_ypos),fill_value=-1.0)
                
                axs[i,j].errorbar(err_xpos,err_ypos,xerr=xerr,yerr=yerr,linestyle='none',elinewidth=1,capsize=0,color='black')
        
            
        
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

                axs[i,j].invert_yaxis()
                #axs[i,j].set_ylim(bottom=20,top=11)
                #axs[i,j].set_xlim(left=-1)
                axs[i,j].xaxis.set_major_locator(plt.MaxNLocator(6))
                axs[i,j].yaxis.set_major_locator(plt.MaxNLocator(6))
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
        
        params={'legend.fontsize':'12','axes.labelsize':'18',
        'axes.titlesize':'12','xtick.labelsize':'12',
        'ytick.labelsize':'12','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        fig,axs=plt.subplots(2,2,sharex=True,sharey=True,figsize=[6,6])
        
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
               
                axs[i,j].xaxis.set_minor_locator(MultipleLocator(0.5))
                axs[i,j].yaxis.set_minor_locator(MultipleLocator(0.5))
        axs[1,0].set_xlabel('(J-H)$_0$')
        axs[1,1].set_xlabel('(J-H)$_0$')
        axs[1,0].set_ylabel('(H-K)$_0$')
        axs[0,0].set_ylabel('(H-K)$_0$')
            #i.invert_yaxis()
        
    def plot_kj_cmds_dsphs(self,marker='.',markersize=5,color='black'):
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'18',
        'axes.titlesize':'14','xtick.labelsize':'12',
        'ytick.labelsize':'12','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        xdim=4
        ydim=3
        
        x=0
        y=0
        
        fig,axs=plt.subplots(ydim,xdim,sharex=True,sharey=True)

        for i in self.sphs:
            
            if i.galaxy=='and18':
                
                continue
            
            sns.scatterplot(i.data.jmag-i.data.kmag,i.data.kmag,marker=marker,linewidth=0,s=markersize,ax=axs[y,x],label=i.galaxy,color='black')
            
            leg = axs[y,x].legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
            #axs[y,x].invert_yaxis()
            
            if x==0:
                axs[y,x].set_ylabel('K$_0$')
                
            if y==ydim-1:
                
                axs[y,x].set_xlabel('(J-K)$_0$')
                
            
            if x==xdim-1:
                
                x=0
                y=y+1
            else:
                x=x+1
                
        for i in range(ydim):
            for j in range(xdim):
                axs[i,j].invert_yaxis()
                
    def plot_ccs_dsphs(self,marker='.',markersize=5,color='black'):
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'18',
        'axes.titlesize':'14','xtick.labelsize':'12',
        'ytick.labelsize':'12','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        xdim=4
        ydim=3
        
        x=0
        y=0
        
        fig,axs=plt.subplots(ydim,xdim,sharex=True,sharey=True)

        for i in self.sphs:
            
            if i.galaxy=='and18':
                
                continue
            
            sns.scatterplot(i.data.jmag-i.data.hmag,i.data.hmag-i.data.kmag,marker=marker,linewidth=0,s=markersize,ax=axs[y,x],label=i.galaxy,color='black')
            
            leg = axs[y,x].legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
            #axs[y,x].invert_yaxis()
            
            if x==0:
                axs[y,x].set_ylabel('(H-K)$_0$')
                
            if y==ydim-1:
                
                axs[y,x].set_xlabel('(J-H)$_0$')
                
            
            if x==xdim-1:
                
                x=0
                y=y+1
            else:
                x=x+1
                
    def plot_lums_dEs(self):
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'15',
        'axes.titlesize':'12','xtick.labelsize':'10',
        'ytick.labelsize':'10','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        xdim=2
        ydim=2
        
        x=0
        y=0
        
        fig,axs=plt.subplots(ydim,xdim,sharex=True,sharey=True)
        

        for i in self.dEs:
            
            bin_no=int((np.max(i.data.kmag.dropna())-np.min(i.data.kmag.dropna()))/0.2)
            
            sns.distplot(i.data.kmag,bins=bin_no,ax=axs[y,x],kde=False,norm_hist=False,hist_kws={'histtype':'step','linewidth':1.5,'alpha':1,'color':sns.xkcd_rgb["black"]},label=i.galaxy.upper()+ ' <2\'')
            
            leg = axs[y,x].legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
            #axs[y,x].invert_yaxis()


            if x==0:
                axs[y,x].set_ylabel('No. Sources')
                
            if y==ydim-1:
                
                axs[y,x].set_xlabel('K$_0$')
                
            
            if x==xdim-1:
                
                x=0
                y=y+1
            else:
                x=x+1
                
        for i in range(ydim):
            for j in range(xdim):
                axs[i,j].invert_xaxis()
                axs[i,j].xaxis.set_minor_locator(MultipleLocator(0.5))
        
        plt.yscale('log')
        plt.subplots_adjust(wspace=0, hspace=0)

    def plot_spatials_dsphs(self,marker='.',markersize=5,color='black'):
        
        sns.set_context('paper')
        
        params={'legend.fontsize':'12','axes.labelsize':'18',
        'axes.titlesize':'14','xtick.labelsize':'12',
        'ytick.labelsize':'12','lines.linewidth':2,'axes.linewidth':2,'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        
        xdim=4
        ydim=3
        
        x=0
        y=0
        
        fig,axs=plt.subplots(ydim,xdim,sharex=True,sharey=True)

        for i in self.sphs:
            
            if i.galaxy=='and18':
                
                continue
            
            sns.scatterplot(i.data.xi,i.data.eta,marker=marker,linewidth=0,s=markersize,ax=axs[y,x],color='black')
            
            leg = axs[y,x].legend(handlelength=0, handletextpad=0, frameon=False)
            for item in leg.legendHandles:
                item.set_visible(False)
            #axs[y,x].invert_yaxis()
            axs[y,x].set_title(i.galaxy)
            if x==0:
                axs[y,x].set_ylabel(r'$\eta$')
                
            if y==ydim-1:
                
                axs[y,x].set_xlabel(r'$\xi$')
                
            
            if x==xdim-1:
                
                x=0
                y=y+1
            else:
                x=x+1
                
        for i in range(ydim):
            for j in range(xdim):
                axs[i,j].invert_xaxis()
                
        
        #plt.subplots_adjust(wspace=0, hspace=0)
        
        
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
        
        
        
        
        
    

