from sat_graphers import sat

class hists(sat):

    def plot_hist_slices(self,maxm=20,minm=12,slicemag=0.3,bins=30):
        
        slices=np.linspace(minm,maxm,int((maxm-minm)/slicemag),endpoint=False)

        colours=[]
        kmags=[]
        tdata=self.data.copy()
        
        for i in slices:
            
            data=tdata.copy()
            colour=data.jmag-data.kmag

            length=[]
            for j in range(len(data.kmag)):

                
                if np.isnan(data.kmag[j])==True:
                    continue
                    
                
                elif data.kmag[j] > (i + slicemag) or data.kmag[j] < i:

                    data.loc[j]=np.nan
                    #data.jmag[j]=np.nan
                    

            

            colours.append(data.jmag-data.kmag)
            kmags.append(data.kmag)
        for i in range(len(kmags)):
            
            fig,axs=plt.subplots(1,2)
            xdata=colours[i].copy()

        
            ax1=sns.distplot(xdata.dropna(),bins=bins,kde=True,ax=axs[0])
            ax2=sns.scatterplot(tdata.jmag-tdata.kmag,tdata.kmag,s=10,ax=axs[1])
            ax3=sns.scatterplot(colours[i],kmags[i],s=10,ax=axs[1])
            ax2.invert_yaxis()