#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""
from sat_graphers import sat
from runners import run_both,run_rgb,run_cross,kde_separator
from graphing_class import basic_graphs,graphs
from crossmatching_utils import topcatstuff
from iso_utils import import_isos
from designations import ratio_utils
from data_load import data_load
from data_read import data_read
from data_readall import data_readall
from cuts import cuts
from cuts_all import cuts_all
from bookkeeping import bookkeeping
from edge_detectors import edge_detectors
from spitzer import spitzer
from isochrones import isochrones
from crossmatch_stilts import crossmatch
from selection_utils import selection_utils
from stack import stack
from HESSCMD import plotHess
import numpy as np
import matplotlib.pyplot as plt


#class for reading in and plotting ascii file data

##185 centre: 9.7415417, 48.3373778
##147 centre: 8.3005, 48.5087389
##205 centre : 10.09189356, 41.68541564
####32 centre 10.6742708, 40.8651694

#andromeda, roughly 8 deg field

#for kde method, start at 18, go upwards in bins of mag=0.4

def main():

    #e=data_readall(stage='cls_cut')
    #e.plot_kj_cmds()
    

    
    f=data_read(stage='agb')
    f.plot_spatial()
    e=selection_utils()
    f.data = e.select_ellipse(f.data,afl=0.2,eccentricity=0.9,clockrot=0)
    f.plot_spatial(color='red')
  





    #f.select_stars((1.5,12.8),(1.15,18.8),graph='kj_cmd')
    
    #f.plot_spatial()



    #g=ratio_utils('ngc205')
    #g.define_marginals()
    #g.CM_marginals()

    #mdata=g.hkMdata
    #cdata=g.hkCdata
    
    #plt.scatter(cdata.hmag-cdata.kmag,cdata.jmag-cdata.hmag)
    

    
        
    

    

    

    
    
    

    

    

    
#0.25
    
    

    
main()


#radial gradiant: C/M = 0.22 <70'', C/M = 0.26 >70''

#ho Nhung paper