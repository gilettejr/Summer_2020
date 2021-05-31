from data_read import data_reader
from selection_utils import selection_utils

from astropy.modeling.models import Sersic1D
from astropy.modeling import fitting

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import pandas as pd

import numpy as np


class data_processor(data_reader):

    def construct_slices(self, stars='agb', a_width=0.02):

        # read in AGB data together, and individual C and M catalogues
        if stars == 'agb':
            data = self.data

        elif stars == 'c':

            data = self.cdata

        elif stars == 'm':

            data = self.mdata

        outer_rads = [0.24, 0.12, 0.2, 0.16]

        self.a_width = a_width

        self.stars = stars

        # list of eccentricities of galaxies
        eccentricities = [0.46, 0.22, 0.43, 0.14]
        # list of rotation angles of galaxies
        rotations = [34.2, 45.9, 169.2, 157.9]

        for i in range(len(self.galaxies)):
            if self.galaxy == self.galaxies[i]:

                eccentricity = eccentricities[i]
                rotation = rotations[i]
                outer_rad = outer_rads[i]
        # redundant, for estimating eccentricity
        self.outer_rad = outer_rad
        # initialise selection_utils instance for constructing ellipses
        check = selection_utils()
        # list for holding data slices
        slices = []
        areas = []
        ellipse_shapes = []
        # fill slices with elliptical selections of decreasing size, from outer_rad to a_width
        for i in range(int((outer_rad*1000)/(a_width*1000))):

            ellipse = check.select_ellipse(
                data, afl=outer_rad-(a_width * i), eccentricity=eccentricity, clockrot=rotation)
            slices.append(ellipse[0])
            areas.append(ellipse[1])
            ellipse_shapes.append(ellipse[2])

        self.ellipticity = eccentricity

        for i in range(len(slices)-1):

            areas[i] = areas[i]-areas[i+1]
            intersection = ellipse_shapes[i].intersection(ellipse_shapes[i+1])
            ellipse_shapes[i] = ellipse_shapes[i].difference(intersection)
            # looks ugly but far faster than using more nested loops

            # produce list of common points between outer and inner ellipse at each element in slices list
            locs = np.where(slices[i].orig_index == slices[i+1].orig_index)
            # np.where produces extra array for some reason
            # remove that
            locs = locs[0]
            # match index in array to original pandas index
            indices = []
            for j in locs:
                indices.append(slices[i].index[j])

            # wipe all common values in outer slice
            for k in indices:
                slices[i].loc[k] = np.nan

        for i in slices:

            i = i.dropna()

        self.areas = areas
        self.slices = slices
        self.slice_shapes = ellipse_shapes

    def find_FEH_slices(self, m_background, c_background):

        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        # read in AGB data together, and individual C and M catalogues
        cdata = self.cdata
        mdata = self.mdata
        slices = self.slices
        outer_rad = self.outer_rad
        a_width = self.a_width

        if self.galaxy == 'ngc147':

            crowding_num = 1

        elif self.galaxy == 'ngc185':

            crowding_num = 1

        # list of eccentricities of galaxies

        mslices = []
        cslices = []

        # need to create copies to prevent dataframes being
        # continuously altered

        for i in slices:

            mslices.append(i.copy())
            cslices.append(i.copy())

        # find locations of m/c data in slices

        for i in range(len(slices)):

            clocs = np.where(mslices[i].orig_index == mdata.orig_index)
            mlocs = np.where(cslices[i].orig_index == cdata.orig_index)

            # match array index to pandas index

            mindices = []
            cindices = []

            for j in mlocs:

                mindices.append(cslices[i].index[j])

            for j in clocs:

                cindices.append(mslices[i].index[j])

            # wipe C stars from m data slices

            for k in mindices:

                mslices[i].loc[k] = np.nan

            # wipe M stars from c data slices

            for k in cindices:

                cslices[i].loc[k] = np.nan

            mslices[i] = mslices[i].dropna()
            cslices[i] = cslices[i].dropna()

        # construct new lists to find c/m ratio

        mnum = []
        cnum = []

        # find number of c and m stars in each slices

        for i in range(len(slices)):
            mnum.append(len(mslices[i]))
            cnum.append(len(cslices[i]))

        mnum = np.array(mnum)
        cnum = np.array(cnum)

        # create arrray of c/m ratios in each slices

        cm = cnum/mnum

        c_unc = np.sqrt(cnum)
        m_unc = np.sqrt(mnum)

        back_num = c_background[2]/m_background[2]

        back_cm = c_background[0]/m_background[0]

        slice_all = cnum+mnum

        all_stars = np.sum(slice_all)

        weights = (slice_all/all_stars)

        back_m_uncs = np.sqrt(m_background[2])
        back_c_uncs = np.sqrt(c_background[2])

        back_cm_unc = back_cm * \
            np.sqrt((back_m_uncs/m_unc)**2 + (back_c_uncs/c_unc)**2)

        cm_slice_unc = cm * np.sqrt((c_unc/cnum)**2 + (m_unc/mnum)**2)

        cm_slice_unc = np.sqrt(cm_slice_unc**2 + back_cm_unc**2)

        cm = cm-(back_cm * (back_num/all_stars))

        avg_cm = np.average(cm[:(len(cm)-crowding_num)],
                            weights=weights[:(len(cm)-crowding_num)])

        avg_cm_unc = np.sqrt(np.sum(cm_slice_unc**2))/len(cm_slice_unc)

        print(avg_cm)

        print(avg_cm_unc)

        # convert to [Fe/H]

        FEH = self.CM_to_FEH(cm)

        FEH_unc = self.CM_to_FEH_uncs(cm, cm_slice_unc)

        avgFEH = self.CM_to_FEH(avg_cm)

        avgFEH_unc = self.CM_to_FEH_uncs(avg_cm, avg_cm_unc)

        print('Average [Fe/H]:')
        print(avgFEH)
        print('Uncertainty:')
        print(avgFEH_unc)

        # plot radial distribution

        xdata = np.linspace(outer_rad-a_width/2, 0+a_width/2, num=len(FEH))
        xdata = xdata*60
        plt.errorbar(xdata, FEH, yerr=FEH_unc, capsize=2,
                     linestyle='none', marker='o', markersize='5', color='black')
        plt.plot(xdata, FEH, linestyle='none', markersize='1',
                 marker='o', label=self.galaxy.upper(), color='black')
        FEH_tab = pd.DataFrame({'FEH': FEH, 'FEH_unc': FEH_unc, 'a': xdata})
        FEH_tab.to_csv(self.galaxy+'_FEH_table')
        # i.label_outer()
        leg = plt.legend(handlelength=0, handletextpad=0,
                         frameon=False, loc='upper left', markerscale=0.001)
        for item in leg.legendHandles:
            item.set_visible(False)
        plt.ylabel('[Fe/H]')
        plt.xlabel('Semi-major axis (arcmins)')

        plt.savefig(self.galaxy + '_radial_FEH')

        m, b = np.polyfit(xdata, FEH, 1)
        #plt.plot(xdata, m*xdata + b,color='red')

        return xdata, FEH, FEH_unc, avgFEH, avgFEH_unc

    def fit_slice_count_profile(self, background_deg, crowding_num=0):

        slices = self.slices
        areas = self.areas
        outer_rad = self.outer_rad
        a_width = self.a_width
        ellipticity = self.ellipticity
        slice_nums = []

        for i in slices:
            slice_nums.append(len(i.dropna()))

        slice_nums = np.array(slice_nums)

        areas = np.array(areas)

        slice_count_uncs = np.sqrt(slice_nums)

        slice_star_densities = slice_nums/areas
        slice_star_densities_uncs = (
            np.sqrt((slice_count_uncs/areas)**2 + background_deg[1]**2))

        slice_star_densities = slice_star_densities-background_deg[0]

        slice_star_densities = slice_star_densities/3600
        slice_star_densities_uncs = slice_star_densities_uncs/3600

        xdata = np.linspace(outer_rad-a_width/2, 0+a_width/2,
                            num=(outer_rad*1000)/(a_width*1000))
        xdata = xdata*60
        # convert from a to r_eff
        ydata = slice_star_densities
        yerr = slice_star_densities_uncs

        #plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',marker='o',markersize=5,capsize=2,color='black',label='AGB data')
        #plt.xlabel('Semi-major axis/arcmins')
        # plt.ylabel('Nstars/arcmins$^2$')

        # n147 params

        xdata_orig = xdata
        ydata_orig = ydata
        yerr_orig = yerr

        xdata = xdata[:len(xdata)-(crowding_num)]
        ydata = ydata[:len(ydata)-(crowding_num)]
        yerr = yerr[:len(yerr)-(crowding_num)]

        if self.galaxy == 'ngc147':

            model = Sersic1D(r_eff=5.2, n=1.8)

        elif self.galaxy == 'ngc185':

            model = Sersic1D(r_eff=1.7, n=1.00)

        fit = fitting.LevMarLSQFitter()

        s = fit(model, xdata, ydata, weights=1/yerr)

        xdata = np.linspace(min(xdata_orig), max(xdata_orig), 1000)

        #plt.plot(xdata,s(xdata),label='Sersic Fit with r$_{eff}$ = ' + str(round(s.r_eff.value,3)) + ', n = ' + str(round(s.n.value,3)))

        # plt.legend()
        print(s)
        self.r_eff = s.r_eff.value
        self.n = s.n.value

        return xdata_orig, ydata_orig, yerr_orig, xdata, s(xdata)

        # plt.xlabel('a/deg')
        # plt.ylabel('[Fe/H]')

    def overplot_ellipses(self, ellipticity, a_inner, a_outer, PA, marker='o', markersize=1, color='black'):

        data = self.data
        majors = np.linspace(a_inner, a_outer, num=int(
            (a_outer*1000)/(a_inner*1000)))
        minors = majors*(1-ellipticity)

        ells = []

        for i in range(len(majors)):
            ells.append(Ellipse(xy=[0, 0], height=majors[i]*2, width=minors[i]*2, angle=360 -
                        PA, facecolor='none', edgecolor='red', linestyle='--', linewidth=1))

        plt.rc('axes', labelsize=20)
        fig, ax = plt.subplots()
        ax.plot(data.xi, data.eta, linestyle='none', marker=marker,
                markersize=markersize, color=color, zorder=1)
        for i in ells:
            ax.add_artist(i)
        ax.invert_xaxis()
        ax.set_ylabel(r'$\eta$')
        ax.set_xlabel(r'$\xi$')
from data_read import data_reader
from selection_utils import selection_utils

from astropy.modeling.models import Sersic1D
from astropy.modeling import fitting

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import pandas as pd

import numpy as np


class data_processor(data_reader):

    def construct_slices(self, stars='agb', a_width=0.02):

        # read in AGB data together, and individual C and M catalogues
        if stars == 'agb':
            data = self.data

        elif stars == 'c':

            data = self.cdata

        elif stars == 'm':

            data = self.mdata

        outer_rads = [0.24, 0.12, 0.2, 0.16]

        self.a_width = a_width

        self.stars = stars

        # list of eccentricities of galaxies
        eccentricities = [0.46, 0.22, 0.43, 0.14]
        # list of rotation angles of galaxies
        rotations = [34.2, 45.9, 169.2, 157.9]

        for i in range(len(self.galaxies)):
            if self.galaxy == self.galaxies[i]:

                eccentricity = eccentricities[i]
                rotation = rotations[i]
                outer_rad = outer_rads[i]
        # redundant, for estimating eccentricity
        self.outer_rad = outer_rad
        # initialise selection_utils instance for constructing ellipses
        check = selection_utils()
        # list for holding data slices
        slices = []
        areas = []
        ellipse_shapes = []
        # fill slices with elliptical selections of decreasing size, from outer_rad to a_width
        for i in range(int((outer_rad*1000)/(a_width*1000))):

            ellipse = check.select_ellipse(
                data, afl=outer_rad-(a_width * i), eccentricity=eccentricity, clockrot=rotation)
            slices.append(ellipse[0])
            areas.append(ellipse[1])
            ellipse_shapes.append(ellipse[2])

        self.ellipticity = eccentricity

        for i in range(len(slices)-1):

            areas[i] = areas[i]-areas[i+1]
            intersection = ellipse_shapes[i].intersection(ellipse_shapes[i+1])
            ellipse_shapes[i] = ellipse_shapes[i].difference(intersection)
            # looks ugly but far faster than using more nested loops

            # produce list of common points between outer and inner ellipse at each element in slices list
            locs = np.where(slices[i].orig_index == slices[i+1].orig_index)
            # np.where produces extra array for some reason
            # remove that
            locs = locs[0]
            # match index in array to original pandas index
            indices = []
            for j in locs:
                indices.append(slices[i].index[j])

            # wipe all common values in outer slice
            for k in indices:
                slices[i].loc[k] = np.nan

        for i in slices:

            i = i.dropna()

        self.areas = areas
        self.slices = slices
        self.slice_shapes = ellipse_shapes

    def find_FEH_slices(self, m_background, c_background):

        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        # read in AGB data together, and individual C and M catalogues
        cdata = self.cdata
        mdata = self.mdata
        slices = self.slices
        outer_rad = self.outer_rad
        a_width = self.a_width

        if self.galaxy == 'ngc147':

            crowding_num = 1

        elif self.galaxy == 'ngc185':

            crowding_num = 1

        # list of eccentricities of galaxies

        mslices = []
        cslices = []

        # need to create copies to prevent dataframes being
        # continuously altered

        for i in slices:

            mslices.append(i.copy())
            cslices.append(i.copy())

        # find locations of m/c data in slices

        for i in range(len(slices)):

            clocs = np.where(mslices[i].orig_index == mdata.orig_index)
            mlocs = np.where(cslices[i].orig_index == cdata.orig_index)

            # match array index to pandas index

            mindices = []
            cindices = []

            for j in mlocs:

                mindices.append(cslices[i].index[j])

            for j in clocs:

                cindices.append(mslices[i].index[j])

            # wipe C stars from m data slices

            for k in mindices:

                mslices[i].loc[k] = np.nan

            # wipe M stars from c data slices

            for k in cindices:

                cslices[i].loc[k] = np.nan

            mslices[i] = mslices[i].dropna()
            cslices[i] = cslices[i].dropna()

        # construct new lists to find c/m ratio

        mnum = []
        cnum = []

        # find number of c and m stars in each slices

        for i in range(len(slices)):
            mnum.append(len(mslices[i]))
            cnum.append(len(cslices[i]))

        mnum = np.array(mnum)
        cnum = np.array(cnum)

        # create arrray of c/m ratios in each slices

        cm = cnum/mnum

        c_unc = np.sqrt(cnum)
        m_unc = np.sqrt(mnum)

        back_num = c_background[2]/m_background[2]

        back_cm = c_background[0]/m_background[0]

        slice_all = cnum+mnum

        all_stars = np.sum(slice_all)

        weights = (slice_all/all_stars)

        back_m_uncs = np.sqrt(m_background[2])
        back_c_uncs = np.sqrt(c_background[2])

        back_cm_unc = back_cm * \
            np.sqrt((back_m_uncs/m_unc)**2 + (back_c_uncs/c_unc)**2)

        cm_slice_unc = cm * np.sqrt((c_unc/cnum)**2 + (m_unc/mnum)**2)

        cm_slice_unc = np.sqrt(cm_slice_unc**2 + back_cm_unc**2)

        cm = cm-(back_cm * (back_num/all_stars))

        avg_cm = np.average(cm[:(len(cm)-crowding_num)],
                            weights=weights[:(len(cm)-crowding_num)])

        avg_cm_unc = np.sqrt(np.sum(cm_slice_unc**2))/len(cm_slice_unc)

        print(avg_cm)

        print(avg_cm_unc)

        # convert to [Fe/H]

        FEH = self.CM_to_FEH(cm)

        FEH_unc = self.CM_to_FEH_uncs(cm, cm_slice_unc)

        avgFEH = self.CM_to_FEH(avg_cm)

        avgFEH_unc = self.CM_to_FEH_uncs(avg_cm, avg_cm_unc)

        print('Average [Fe/H]:')
        print(avgFEH)
        print('Uncertainty:')
        print(avgFEH_unc)

        # plot radial distribution

        xdata = np.linspace(outer_rad-a_width/2, 0+a_width/2, num=len(FEH))
        xdata = xdata*60
        plt.errorbar(xdata, FEH, yerr=FEH_unc, capsize=2,
                     linestyle='none', marker='o', markersize='5', color='black')
        plt.plot(xdata, FEH, linestyle='none', markersize='1',
                 marker='o', label=self.galaxy.upper(), color='black')
        FEH_tab = pd.DataFrame({'FEH': FEH, 'FEH_unc': FEH_unc, 'a': xdata})
        FEH_tab.to_csv(self.galaxy+'_FEH_table')
        # i.label_outer()
        leg = plt.legend(handlelength=0, handletextpad=0,
                         frameon=False, loc='upper left', markerscale=0.001)
        for item in leg.legendHandles:
            item.set_visible(False)
        plt.ylabel('[Fe/H]')
        plt.xlabel('Semi-major axis (arcmins)')

        plt.savefig(self.galaxy + '_radial_FEH')

        m, b = np.polyfit(xdata, FEH, 1)
        #plt.plot(xdata, m*xdata + b,color='red')

        return xdata, FEH, FEH_unc, avgFEH, avgFEH_unc

    def fit_slice_count_profile(self, background_deg, crowding_num=0):

        slices = self.slices
        areas = self.areas
        outer_rad = self.outer_rad
        a_width = self.a_width
        ellipticity = self.ellipticity
        slice_nums = []

        for i in slices:
            slice_nums.append(len(i.dropna()))

        slice_nums = np.array(slice_nums)

        areas = np.array(areas)

        slice_count_uncs = np.sqrt(slice_nums)

        slice_star_densities = slice_nums/areas
        slice_star_densities_uncs = (
            np.sqrt((slice_count_uncs/areas)**2 + background_deg[1]**2))

        slice_star_densities = slice_star_densities-background_deg[0]

        slice_star_densities = slice_star_densities/3600
        slice_star_densities_uncs = slice_star_densities_uncs/3600

        xdata = np.linspace(outer_rad-a_width/2, 0+a_width/2,
                            num=(outer_rad*1000)/(a_width*1000))
        xdata = xdata*60
        # convert from a to r_eff
        ydata = slice_star_densities
        yerr = slice_star_densities_uncs

        #plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',marker='o',markersize=5,capsize=2,color='black',label='AGB data')
        #plt.xlabel('Semi-major axis/arcmins')
        # plt.ylabel('Nstars/arcmins$^2$')

        # n147 params

        xdata_orig = xdata
        ydata_orig = ydata
        yerr_orig = yerr

        xdata = xdata[:len(xdata)-(crowding_num)]
        ydata = ydata[:len(ydata)-(crowding_num)]
        yerr = yerr[:len(yerr)-(crowding_num)]

        if self.galaxy == 'ngc147':

            model = Sersic1D(r_eff=5.2, n=1.8)

        elif self.galaxy == 'ngc185':

            model = Sersic1D(r_eff=1.7, n=1.00)

        fit = fitting.LevMarLSQFitter()

        s = fit(model, xdata, ydata, weights=1/yerr)

        xdata = np.linspace(min(xdata_orig), max(xdata_orig), 1000)

        #plt.plot(xdata,s(xdata),label='Sersic Fit with r$_{eff}$ = ' + str(round(s.r_eff.value,3)) + ', n = ' + str(round(s.n.value,3)))

        # plt.legend()
        print(s)
        self.r_eff = s.r_eff.value
        self.n = s.n.value

        return xdata_orig, ydata_orig, yerr_orig, xdata, s(xdata)

        # plt.xlabel('a/deg')
        # plt.ylabel('[Fe/H]')

    def overplot_ellipses(self, ellipticity, a_inner, a_outer, PA, marker='o', markersize=1, color='black'):

        data = self.data
        majors = np.linspace(a_inner, a_outer, num=int(
            (a_outer*1000)/(a_inner*1000)))
        minors = majors*(1-ellipticity)

        ells = []

        for i in range(len(majors)):
            ells.append(Ellipse(xy=[0, 0], height=majors[i]*2, width=minors[i]*2, angle=360 -
                        PA, facecolor='none', edgecolor='red', linestyle='--', linewidth=1))

        plt.rc('axes', labelsize=20)
        fig, ax = plt.subplots()
        ax.plot(data.xi, data.eta, linestyle='none', marker=marker,
                markersize=markersize, color=color, zorder=1)
        for i in ells:
            ax.add_artist(i)
        ax.invert_xaxis()
        ax.set_ylabel(r'$\eta$')
        ax.set_xlabel(r'$\xi$')
