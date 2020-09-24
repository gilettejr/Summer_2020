#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""
from sat_graphers import sat
#from runners import run_both,run_rgb,run_cross,kde_separator
#from iso_utils import import_isos
#from designations import ratio_utils
from data_load import data_loader
from data_read import data_reader
from data_readall import data_readall
from cuts import cuts
from cuts_all import cuts_all
from bookkeeping import bookkeeping
from edge_detectors import edge_detectors
from spitzer import spitzer
from isochrones import isochrones
from crossmatch_stilts import crossmatch
from selection_utils import selection_utils
#from stack import stack
from HESSCMD import plotHess
from sersic_graphs import radial_graphs
from reload_all_data import reload_all_data
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
    
    reload_all_data()
    
    #e.make_slices(outer_rad=0.4)
    #e.overplot_ellipses(0.46,0.02,0.4,34.2)
    #e.slice_count_profile(background_deg=background,crowding_n'Semi-majorum=1)
    #e.slice_count_profile(background,crowding_num=1)
    #e.plot_gaia_removal()
    
    #e=data_readall(stage='agb')
    #e.plot_cc_hess()
    #outer_rad n147, 0.24, n185 0.12 n205 0.2 m32=0.16

    #f.make_slices(stars='agb',a_width=0.02,outer_rad=0.20)
    #m_background=f.find_background_density_border(stars='m')
   # c_background=f.find_background_density_border(stars='c')
    #f.FEH_slices(m_background,c_background)

    #f.make_slices(a_width=0.02,outer_rad=0.2)
    #f.make_slices(a_width=0.02,outer_rad=0.2)

    #background_corner=f.find_background_density_boxes()
    #background_border=f.find_background_density_border()
    

    #n147 coords
    #f.overplot_ellipse(0.46,0.3,34.2)
    #f.overplot_ellipses(0.46,0.02,0.4,34.2)
    #f.make_slices(a_width=0.05,outer_rad=0.4)
    #f.overplot_ellipse
    ###n185 coords
    #f.overplot_ellipse(0.22,0.3,45.9)
    #f.overplot_ellipses(0.22,0.02,0.2,45.9)
    
    #n205 coords

    #f.close_FEH_slices()
    
    #f.make_slices(stars='c',a_width=0.02,outer_rad=0.3)
    #f.find_background_grad()
    #f.fit_close_background()

    
    #f.fit_close_profiles()


    #m32
    #f.make_slices(stars='agb',a_width=0.02,outer_rad=0.16)
    #f.close_FEH_slices()
    #f.fit_close_profiles()
    #f.overplot_ellipse(ellipticity=0.14,a=0.1,PA=157.9)
    

  
    


    #f.slice_mag_profile()


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