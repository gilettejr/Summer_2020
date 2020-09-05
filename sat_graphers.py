from graphing_class import graphs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
class sat(graphs):

    def __init__(self,galaxy,pandas=False,spitzer=False):
        
        def quadrat(a,b):
            
            c=np.sqrt(a**2+b**2)
            
            return c
        
        def data_select(data):
            if data=='all':
                outdata=self.data
            elif data=='c':
                outdata=self.cdata
            elif data=='m':
                outdata=self.mdata
                
            else:
                print('Invalid dataset chosen in method')
        
            return(outdata)
        
        self.mag=quadrat
        self.data_select=data_select
        infilenames=['ngc147','ngc185','ngc205','m32']
        
        for i in infilenames:
            if galaxy==i:
                infile=('pm33_4rem/' + i)
                break
        
        data=pd.read_parquet(infile)

        if pandas == True:
            
            for i in infilenames:
                if galaxy==i:
                    pinfile=('pandas/' + i + '_pand' )
                    break
                pdata=pd.read_parquet(pinfile)
                self.pdata=pdata
                
                
                
        self.galaxy=galaxy
        self.data=data.dropna()
        self.infilenames=infilenames
    
    def plot_kj_cmd(self,marker='o',markersize=1,color='blue',stars='all'):
        plt.figure()        
        data=self.data_select(stars)
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.kmag,data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        plt.gca().invert_yaxis()
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')    
        
    def plot_cc(self,marker='o',markersize=1,color='blue',stars='all'):
        data=self.data_select(stars)
        plt.rc('axes',labelsize = 15)
        plt.plot(data.hmag - data.kmag,data.jmag-data.hmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        plt.ylabel('$J_0$-$H_0$')
        plt.xlabel('$H_0$-$K_0$')
        
    def plot_jk_cmd(self,marker='o',markersize=1,color='blue',stars='all'):
        plt.figure()
        data=self.data_select(stars)
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.kmag,data.jmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        plt.gca().invert_yaxis()
        plt.ylabel('$J_0$')
        plt.xlabel('$J_0$-$K_0$')          
        
    def plot_tan(self,marker='o',markersize=1,color='black',stars='all'):
        
        data=self.data_select(stars)
        

        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.xi,data.eta,linestyle='none',marker = marker,markersize=markersize,color=color)

        plt.gca().invert_xaxis()
        plt.gca().set_ylabel(r'$\eta$/degrees')
        plt.gca().set_xlabel(r'$\xi$/degrees')
    
    def plot_panda_cmd(self,marker='o',markersize=1,color='black'):
        
        pdata=self.pdata
        
        plt.rc('axes',labelsize = 15)
        plt.plot(pdata.g-pdata.i,pdata.i,linestyle='none',marker=marker,markersize=markersize,color=color)
        
        plt.gca().invert_yaxis()
        plt.ylabel('$i_0$')
        plt.xlabel('$g_0$-$i_0$')

                    
                
            
            
        
