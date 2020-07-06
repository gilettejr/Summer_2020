import numpy as np

#libraries for seperating C and M with non horizontal/perpendicular lines
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#libraries for defining ellipses
import shapely.affinity
#from descartes import PolygonPatch

class selection_utils:
    
    def select_ellipse(self,data,afl=False,atup=False,btup=False,eccentricity=0,clockrot=0):
        newdata=data.copy()
        if afl==False:
            a=np.sqrt(atup[0]**2 + atup[1]**2)
            
        else:
            a=afl*1000
        
        if eccentricity==False:
            b=np.sqrt(btup[0]**2 + btup[1]**2)
            b=b*1000
            eccentricity=np.sqrt(1-(b/a)**2)
            print('Eccentricity= ' + str(eccentricity))
            
        elif btup==False:
            
            b=a * np.sqrt(1-eccentricity**2)
        
        # 1st elem = center point (x,y) coordinates
        # 2nd elem = the two semi-axis values (along x, along y)
        # 3rd elem = angle in degrees between x-axis of the Cartesian base
        #            and the corresponding semi-axis
        ellipse = ((0, 0),(a, b),clockrot)
        
        # Let create a circle of radius 1 around center point:
        circ = shapely.geometry.Point(ellipse[0]).buffer(1)
        
        # Let create the ellipse along x and y:
        ell  = shapely.affinity.scale(circ, int(ellipse[1][0]), int(ellipse[1][1]))
        
        # Let rotate the ellipse (clockwise, x axis pointing right):
        ellr = shapely.affinity.rotate(ell,ellipse[2])
        
        for i in newdata.index:
            
            if ellr.contains(Point(newdata.xi[i]*1000,newdata.eta[i]*1000)) == False:
                
                newdata.loc[i]=np.nan
            
        return newdata
    
    #select subset of stars, defined by corners of a square given in tuple format
    def select_stars(self,data,corner1,corner3,graph='spatial'):
        
        
        corner2=(corner3[0],corner1[1])
        corner4=(corner1[0],corner3[1])
        
        select=Polygon([corner1,corner2,corner3,corner4])
        
        if graph=='spatial':
            xdata=data.xi
            ydata=data.eta
            
        elif graph=='kj_cmd':
            
            xdata=data.jmag-data.kmag
            ydata=data.kmag
            
        elif graph =='cc':
            
            xdata=data.jmag-data.hmag
            ydata=data.hmag-data.kmag
            
        else:
            
            print('Invalid graph format')
            
        
        for i in data.index:
            
            if select.contains(Point(xdata[i],ydata[i]))==False:
                
                data.loc[i]=np.nan
        
        return data.dropna()