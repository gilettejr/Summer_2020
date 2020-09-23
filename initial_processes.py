import numpy as np

class initial_processes:
    
    def __init__(self,galaxy_object):
        
        def get_data_and_name_attributes():
            
            galaxy_data=self.galaxy_data
            galaxy_name=self.galaxy_name
            
            return galaxy_data,galaxy_name
        
        self.galaxy_data=galaxy_object.data
        self.galaxy_name=galaxy_object.galaxy
        
        self.get_data_and_name_attributes=get_data_and_name_attributes

        
        
        
    def CLS_cut(self,cls_bands='norm'):
        
        galaxy_data,galaxy_name=self.get_data_and_name_attributes()
        
        #set of cuts used for normal data processing
        if cls_bands=='norm':
            #m32 treated differently to prevent gradient issues in spatial distribution
            if galaxy_name=='m32':
                for i in range(len(galaxy_data.jcis)):
                    #standard cls cuts made for the rest
                    if (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0 and galaxy_data.jcis[i]!=-3.0) or (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0 and galaxy_data.kcis[i]!=-3.0) or (galaxy_data.hcis[i]==-8.0) or galaxy_data.jmag[i]-galaxy_data.hmag[i] > 17:
                        galaxy_data.loc[i]=np.nan
            
            else:
            

            #standard cls cuts made for the rest
                nines=[]
                for i in range(len(galaxy_data.jcis)):

                    if galaxy_data.kcis[i]==-9.0 or galaxy_data.jcis[i]==-9.0 or galaxy_data.hcis[i]==-9.0:
                        nines.append(0)
                    
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0) or (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0) or (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
        #different combinations of cls cuts included below
        #can be used to test effect of different cls cuts

            
        elif cls_bands=='jh' or cls_bands=='hj':
            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0 and galaxy_data.jcis[i]!=-3.0) or (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0 and galaxy_data.hcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0) or (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan      
                        
        elif cls_bands=='hk':
            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0 and galaxy_data.kcis[i]!=-3.0) or (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0 and galaxy_data.hcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0) or (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
                        
        elif cls_bands =='jk':
            

            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0 and galaxy_data.kcis[i]!=-3.0) or (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0 and galaxy_data.jcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0) or (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
                        
        elif cls_bands == 'j':
            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0 and galaxy_data.jcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
                        
        elif cls_bands == 'h':
            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0 and galaxy_data.hcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
                        
        elif cls_bands == 'k':
            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0 and galaxy_data.kcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
                        
        elif cls_bands=='all':
            
            if self.galaxy=='m32':
                for i in range(len(galaxy_data.jcis)):
                
                    if (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0 and galaxy_data.jcis[i]!=-3.0) or (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0 and galaxy_data.kcis[i]!=-3.0) or (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0 and galaxy_data.hcis[i]!=-3.0):
                        galaxy_data.loc[i]=np.nan
            
            else:
            
            #removes any non-stellar sources from data
                for i in range(len(galaxy_data.jcis)):
                    if (galaxy_data.kcis[i] != -1.0 and galaxy_data.kcis[i]!=-2.0) or (galaxy_data.hcis[i] != -1.0 and galaxy_data.hcis[i]!=-2.0) or (galaxy_data.jcis[i] != -1.0 and galaxy_data.jcis[i]!=-2.0): 
                        galaxy_data.loc[i]=np.nan
                        
        else:
            #failure if invalid cls argument chosen upon class initialisation
            print('Invalid CLS argument chosen, no CLS cut made')
        #percentage decrease of catalogue length calculated and printed
        
        k=len(galaxy_data.jcis)
        l=len(galaxy_data.jcis.dropna())            
        decrease=100-(l/k)*100
    
        print('CLS cut reduced data by ' + str(int(decrease)) + '%')
        
            #function to cut data with large magnitude error
    #takes datagalaxy_data as input, only used internally in this class
    
    def mag_err_cut(self,magerr):
        
        galaxy_data,galaxy_name=self.get_data_and_name_attributes()
        
        #have to start being careful about referencing index column of datagalaxy_data now since dropna() will have been used
        
        for i in galaxy_data.index:
            if galaxy_data.kerr[i] > magerr or galaxy_data.herr[i] > magerr or galaxy_data.kerr[i] > magerr:
            
                galaxy_data.loc[i]=np.nan
        
        #percentage decrease of catalogue length calculated and printed
        
        k=len(galaxy_data.jcis)
        l=len(galaxy_data.jcis.dropna())
        
        decrease=100-(l/k)*100
        
        print('Magerr cut reduced data by ' + str(int(decrease)) + '%')
        
    def make_tan_coord(self):
        
        galaxy_data,galaxy_name=self.get_data_and_name_attributes()
        
        #intermediate function to do the heavy lifting
        def create_tangent_coords(galaxy_data,tangentra,tangentdec):

            #ra and dec attributes converted to variables in radians

            ra = np.radians(galaxy_data.RA)
            dec = np.radians(galaxy_data.DEC)
            
            #tangent co-ordinates also converted to radians
            
            tanra = np.radians(tangentra)
            tandec = np.radians(tangentdec)
            
            #conversion for xi carried out
            
            xi = (np.cos(dec)*np.sin(ra-tanra))/(np.sin(dec)*np.sin(tandec) + np.cos(dec)*np.cos(tandec)*np.cos(ra-tanra))
            
            #conversion for eta carried out
            
            eta = (np.sin(dec)*np.cos(tandec)-np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))/(np.sin(dec)*np.cos(tandec)+np.cos(dec)*np.sin(tandec)*np.cos(ra-tanra))
            
            #co-ordinates converted to degrees and set as attributes
            
            xi = xi * (180/np.pi)
            eta = eta * (180/np.pi)
            
            return (np.array([xi,eta]))
        #ra and dec positions in decimal degrees for each galaxy
        if galaxy_name=='ngc147':
            tra=8.300500
            tdec=48.50850
        elif galaxy_name=='ngc185':
            tra=9.7415417
            tdec=48.3373778
        elif galaxy_name=='ngc205':
            tra=10.09189356
            tdec=41.68541564
        elif galaxy_name=='m32':
            tra=10.6742708
            tdec=40.8651694
        
        elif galaxy_name=='m31':
            tra=10.68470833
            tdec=41.26875
            
        elif galaxy_name=='and1':
            tra=11.41583333
            tdec=38.04111111
            
        elif galaxy_name=='and2':
            tra=19.12416667
            tdec=33.41916667
            
        elif galaxy_name=='and3':
            tra=8.84458333
            tdec=36.50472222
            
        elif galaxy_name=='and6':
            tra=357.94333333
            tdec=24.58638889
            
        elif galaxy_name=='and7':
            tra=351.63208333
            tdec=50.67583333
            
        elif galaxy_name=='and10':
            tra=16.64041667
            tdec=44.80438889
            
            
        elif galaxy_name=='and14':
            tra=12.89583333
            tdec=29.69694444
            
        elif galaxy_name=='and15':
            tra=18.57791667
            tdec=38.1175
            
        elif galaxy_name=='and16':
            tra=14.87416667
            tdec=32.37666667
            
        elif galaxy_name=='and17':
            tra=9.27916667
            tdec=44.32222222
            
        elif galaxy_name=='and18':
            tra=0.56041667
            tdec=45.08888889
            
            
        elif galaxy_name=='and19':
            tra=4.88375
            tdec=35.04363889
            
        elif galaxy_name=='and20':
            tra=1.8779166700000003
            tdec=35.13233333
        #function run dependent on galaxy
        coords=create_tangent_coords(galaxy_data,tra,tdec)
        #output produced and xi and eta columns added to datagalaxy_data
        galaxy_data['xi']=coords[0]
        galaxy_data['eta']=coords[1]
            
#function to cut Datagalaxy_data based on cls index


    #for some reason if I put galaxy_data=galaxy_data.dropna() here it doesn't do anything, so I have to wipe the NaNs outside the function
    
        
    def make_deg_coord(self):
        galaxy_data=self.galaxy_data
        #takes relevant columns as arguements, converts into single ra and single dec column
        
        def coordtransform(rah,ram,ras,dd,dm,ds):
        
            #transformations carried out
            
            #ra
            
            long = (rah + np.divide(ram,60) + np.divide(ras,3600))*15
            
            #dec
            
            lat = dd + np.divide(dm,60) + np.divide(ds,3600)
            
            #array containing transformed co-ordinates returned
    
            return(np.array([long,lat]))
            
        #function called    
        
        coords=coordtransform(galaxy_data.rah,galaxy_data.ram,galaxy_data.ras,galaxy_data.dd,galaxy_data.dm,galaxy_data.ds)
            
        #columns added to Datagalaxy_data
            
        galaxy_data['RA']=coords[0]
        
        galaxy_data['DEC']=coords[1]
        
        #applies extinction correction to data in Datagalaxy_data format    
        #again used in the class only
    def ext_corr(galaxy_data):
        

        #attributes set as variables for ease
        
        ra = galaxy_data.RA
        dec = galaxy_data.DEC
        
        #astropy skycoord called, units and galaxy_data set
        
        coords=SkyCoord(ra,dec,unit='deg',galaxy_data='icrs')
        
        #sfd map chosen
        
        sfd=SFDQuery()
        
        #E(B-V) for each star loaded into array from map
        
        sfdred=sfd(coords)
           
        #corrections from table 6 of Schlafly and Finkbeiner(2011) made
        
        jext = sfdred * 0.709
        hext = sfdred * 0.449
        kext = sfdred * 0.302
        
        #extinction corrections carried out on class attributes
        
        galaxy_data.jmag=galaxy_data.jmag - jext
        galaxy_data.hmag=galaxy_data.hmag - hext
        galaxy_data.kmag=galaxy_data.kmag - kext
        
        print('UKIRT extinction corrections done')
            
            
            
    
    
    
    
    