import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

class graphing_utils:
    
        
    def plot_kj_cmd(self,data,markersize='1',marker='o',color='black',overlay=False):
        
        #conditional statements plot only c,m, or both sets depending on 
        #optional stars argument
        #axes, figure set, CMD plotted

        #plt.figure()        
        plt.rc('axes',labelsize = 15)
        plt.plot(data.jmag - data.kmag,data.kmag,linestyle='none',markersize=markersize,marker=marker,color=color)
        if overlay==False:
            plt.gca().invert_yaxis()
            
        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')        
    
    