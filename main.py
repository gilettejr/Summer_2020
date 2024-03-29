
"""
Created on Sat Sep 28 15:14:03 2019

@author: cameronrobertson
"""
from FEH_finders import FEH_finders
from basic_agb_plotter import basic_agb_plotter
from four_agb_plotter import four_agb_plotter
from close_data_processor import close_data_processor
# from runners import run_both,run_rgb,run_cross,kde_separator
# from iso_utils import import_isos
# from designations import ratio_utils
from data_load import data_loader
from data_read import data_reader
from data_readall import data_readall
from cuts import cuts
from cuts_all import cuts_all
from bookkeeping import bookkeeping
from edge_detectors import edge_detectors

from background_constructor import background_constructor


from isochrones import isochrones
from crossmatch_stilts import crossmatch
from selection_utils import selection_utils

from background_runner import background_runner
# from stack import stack
from HESSCMD import plotHess

from reload_all_data import reload_all_data
import numpy as np
import matplotlib.pyplot as plt


# class for reading in and plotting ascii file data

# 185 centre: 9.7415417, 48.3373778
# 147 centre: 8.3005, 48.5087389
# 205 centre : 10.09189356, 41.68541564
# 32 centre 10.6742708, 40.8651694

# andromeda, roughly 8 deg field

# for kde method, start at 18, go upwards in bins of mag=0.4

def run_interactive():

    print('Enter the number for what you would like to do')
    # run this initially - saves everything you need for the other options
    print('1 - Process raw UKIRT data from scratch, create intermediate files and background profiles')
    # everything else is in here
    print('2 - Visualise data from intermediate files')

    stage1 = input()

    if stage1 == '1':

        print("Are you sure you wish to process the UKIRT data from scratch (Y/n)")

        checkrerun = input()

        if checkrerun == "Y":
            reload_all_data()
            runner = background_runner()
            runner.make_all_close_backgrounds()
            print('Done, exiting interactive program')
        else:
            run_interactive()

    elif stage1 == '2':

        print('Which stage of processing would you like to visualise')

        print('1 - No colour cuts - only CLS and Mag err cuts have been applied')
        print('2 - Foreground colour cut, leaving only AGB and RGB stars')
        print('3 - Foreground + TRGB cut, leaving only AGB stars')
        print('4 - AGB stars only, plus Gaia crossmatching cut (final data)')

        stages = ['cls_cut', 'fore_cut', 'agb', 'cm']

        stage2 = input()

        stage = stages[int(stage2)-1]

        print('Which of the following are you wanting to do?')

        print('1 - View basic diagrams')
        if stage2 == '4':
            print(
                '2 - View background fitting profiles for NGC205 or M32 or view rotated coordinates')
            print(
                '3 - View [Fe/H] distributions and elliptical bins from saved background subtracted data')

        stage3 = input()

        if stage3 != '2':

            print('And for which galaxy?')

            print('1 - NGC147')
            print('2 - NGC185')
            print('3 - NGC205')
            print('4 - M32')
            print('5 - ALL')

            stage4 = input()

            galaxies = ['ngc147', 'ngc185', 'ngc205', 'm32']

            if stage4 != '5':

                galaxy = galaxies[int(stage4)-1]

            else:

                galaxy = 'all'

        else:
            print('Which galaxy would you like to generate a graph for?')
            print('This will likely take a while')
            print('1 - NGC205')
            print('2 - M32')
            stage4 = input()

            galaxies = ['ngc205', 'm32']
            if stage4 != '3':
                galaxy = galaxies[int(stage4)-1]

        if stage3 == '1':

            if stage4 != '5':

                galaxy_object = basic_agb_plotter(galaxy=galaxy, stage=stage)

                print('Which of the following would you like to view?')
                print('1 - J-K CMD')
                print('2 - H-K CMD')
                print('3 - 2CD')
                print('4 - Spatial distribution')
                print(
                    '5 - Spatial distribution for ngc147 or ngc185 with background border regions indicated')

                stage5 = input()

                if stage5 == '1':

                    galaxy_object.plot_kj_cmd()

                elif stage5 == '2':

                    galaxy_object.plot_hk_cmd()

                elif stage5 == '3':

                    galaxy_object.plot_cc()

                elif stage5 == '4':

                    galaxy_object.plot_spatial()
                elif stage5 == '5':
                    galaxy_object = background_constructor(galaxy=galaxy)
                    galaxy_object.find_background_density_border(
                        stars='agb', show_figure=True)
                    print(
                        'Would you like to also generate a plot of the flat background densities in this region? y/n')
                    choice = input()
                    if choice == 'y' or choice == 'yes':
                        galaxy_object = background_constructor(galaxy=galaxy)
                        plt.figure()
                        agb = galaxy_object.find_background_density_border(
                            stars='agb', show_figure=False)
                        m = galaxy_object.find_background_density_border(
                            stars='m', show_figure=False)
                        c = galaxy_object.find_background_density_border(
                            stars='c', show_figure=False)

                        plt.bar(['agb', 'm', 'c'], [
                                agb[0], m[0], c[0]], color='black')
                        plt.ylabel('sources/degrees$^2$')

            else:

                galaxy_object = four_agb_plotter(stage=stage)

                print(
                    'Which of the following would you like to view? The result will be saved')
                print('1 - J-K Hess CMD with foreground and TRGB cuts indicated')
                print('2 - Hess 2CD with C-M cuts indicated')
                print('3 - Spatial Distributions')

                stage6 = input()

                if stage6 == '1':

                    galaxy_object.plot_kj_cmd_hess(save=True)

                elif stage6 == '2':

                    galaxy_object.plot_cc_hess(save=True)

                elif stage6 == '3':

                    galaxy_object.plot_spatial(save=True)
        elif stage3 == '2':

            if stage4 != '3':

                print(
                    'For which star background? Alternatively, press 4 to view rotated coordinates')

                print('1 - Total AGB background')
                print('2 - M-star background')
                print('3 - C-star background')
                print('4 - View AGB stars in rotated coordinate system')

                stage7 = input()

                if stage7 == '4':
                    runner = background_runner()
                    runner.show_rotated_coords(galaxy=galaxy, stars='agb')

                else:

                    starss = ['agb', 'm', 'c']

                    stars = starss[int(stage7)-1]

                    runner = background_runner()
                    runner.make_close_backgrounds(galaxy=galaxy, stars=stars)

        elif stage3 == '3':

            if galaxy == 'ngc147' or galaxy == 'ngc185':

                FEH_finders.find_far_FEH(galaxy=galaxy)

            elif galaxy == 'ngc205' or galaxy == 'm32':

                FEH_finders.find_close_FEH(galaxy=galaxy)


def main():

    run_interactive()

    # e.make_slices(outer_rad=0.4)
    # e.overplot_ellipses(0.46,0.02,0.4,34.2)
    # e.slice_count_profile(background_deg=background,crowding_n'Semi-majorum=1)
    # e.slice_count_profile(background,crowding_num=1)
    # e.plot_gaia_removal()

    # e=data_readall(stage='agb')
    # e.plot_cc_hess()
    # outer_rad n147, 0.24, n185 0.12 n205 0.2 m32=0.16

    # f.make_slices(stars='agb',a_width=0.02,outer_rad=0.20)
    # m_background=f.find_background_density_border(stars='m')
   # c_background=f.find_background_density_border(stars='c')
    # f.FEH_slices(m_background,c_background)

    # f.make_slices(a_width=0.02,outer_rad=0.2)
    # f.make_slices(a_width=0.02,outer_rad=0.2)

    # background_corner=f.find_background_density_boxes()
    # background_border=f.find_background_density_border()

    # n147 coords
    # f.overplot_ellipse(0.46,0.3,34.2)
    # f.overplot_ellipses(0.46,0.02,0.4,34.2)
    # f.make_slices(a_width=0.05,outer_rad=0.4)
    # f.overplot_ellipse
    # n185 coords
    # f.overplot_ellipse(0.22,0.3,45.9)
    # f.overplot_ellipses(0.22,0.02,0.2,45.9)

    # n205 coords

    # f.close_FEH_slices()

    # f.make_slices(stars='c',a_width=0.02,outer_rad=0.3)
    # f.find_background_grad()
    # f.fit_close_background()

    # f.fit_close_profiles()

    # m32
    # f.make_slices(stars='agb',a_width=0.02,outer_rad=0.16)
    # f.close_FEH_slices()
    # f.fit_close_profiles()
    # f.overplot_ellipse(ellipticity=0.14,a=0.1,PA=157.9)

    # f.slice_mag_profile()

    # f.select_stars((1.5,12.8),(1.15,18.8),graph='kj_cmd')

    # f.plot_spatial()

    # g=ratio_utils('ngc205')
    # g.define_marginals()
    # g.CM_marginals()

    # mdata=g.hkMdata
    # cdata=g.hkCdata

    # plt.scatter(cdata.hmag-cdata.kmag,cdata.jmag-cdata.hmag)


# 0.25


main()


# radial gradiant: C/M = 0.22 <70'', C/M = 0.26 >70''

# ho Nhung paper
