from data_processor import data_processor
from astropy.modeling import fitting
from astropy.modeling.models import Sersic1D

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import pandas as pd

import shapely.affinity

import os


class close_data_processor(data_processor):

    def get_close_FEH_slices(self):

        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        c_slice_dataframe = pd.read_parquet(
            'unfit_background_corrected_profiles/' + self.galaxy + 'c')

        m_slice_dataframe = pd.read_parquet(
            'unfit_background_corrected_profiles/'+self.galaxy + 'm')

        full_slice_nums = c_slice_dataframe.slice_nums + m_slice_dataframe.slice_nums

        total_star_num = np.sum(full_slice_nums)

        slice_a = c_slice_dataframe.a

        cm_slices = c_slice_dataframe.density/m_slice_dataframe.density
        print(m_slice_dataframe.density)
        print(c_slice_dataframe.density)

        cm_slices_uncs = cm_slices*np.sqrt((c_slice_dataframe.density_err/c_slice_dataframe.density)**2 + (
            m_slice_dataframe.density_err/m_slice_dataframe.density)**2)

        FEH_slices = self.CM_to_FEH(cm_slices)
        # print(cm_slices)
        # print(cm_slices_uncs)
        FEH_slice_uncs = self.CM_to_FEH_uncs(cm_slices, cm_slices_uncs)

        plt.errorbar(slice_a, FEH_slices, marker='o', linestyle='none',
                     color='black', capsize=2, yerr=FEH_slice_uncs)
        plt.plot(slice_a, FEH_slices, marker='o', linestyle='none',
                 color='black', label=self.galaxy.upper(), markersize=1)

        plt.ylabel('[Fe/H]')
        plt.xlabel('Semi-major axis (arcmins)')
        FEH_tab = pd.DataFrame(
            {'FEH': FEH_slices, 'FEH_unc': FEH_slice_uncs, 'a': slice_a, 'C/M': cm_slices, 'C/M_unc': cm_slices_uncs})
        FEH_tab.to_csv(self.galaxy+'_FEH_table')
        leg = plt.legend(handlelength=0, handletextpad=0,
                         frameon=False, loc='upper left', markerscale=0.001)
        for item in leg.legendHandles:
            item.set_visible(False)

        weights = full_slice_nums/total_star_num

        # m_slice_dataframe.density=m_slice_dataframe.density[5:]

        # c_slice_dataframe.density=c_slice_dataframe.density[5:]

        new_cm_slices = []
        new_cm_slices_uncs = []
        new_weights = []
        for i in range(len(cm_slices)):

            if m_slice_dataframe.density[i] < 0 or c_slice_dataframe.density[i] < 0:

                continue

            else:

                new_cm_slices.append(cm_slices[i])
                new_cm_slices_uncs.append(cm_slices_uncs[i])
                new_weights.append(weights[i])

        avg_CM = np.average(new_cm_slices, weights=new_weights)
        avg_CM_unc = np.sqrt(
            np.sum(np.square(new_cm_slices_uncs)))/len(new_cm_slices_uncs)

        FEH_avg = self.CM_to_FEH(avg_CM)

        FEH_avg_unc = self.CM_to_FEH_uncs(avg_CM, avg_CM_unc)

        print('Average [Fe/H]:')
        print(FEH_avg)
        print('Uncertainty:')
        print(FEH_avg_unc)

        plt.savefig(self.galaxy + '_radial_FEH')

        self.FEH_x = slice_a
        self.FEH_y = FEH_slices
        self.FEH_yunc = FEH_slice_uncs

        # plt.errorbar(xdata,ydata,yerr=yerr,linestyle='none',capsize=2,marker='o',color='black')

    def fit_close_slice_count_profile(self, stars='agb', crowding_num=1):

        distribution = pd.read_parquet(
            'unfit_background_corrected_profiles' + self.galaxy + stars)

        xdata = distribution.a
        ydata = distribution.density
        yerr = distribution.density_err

        plt.errorbar(xdata, ydata, yerr=yerr, marker='o',
                     linestyle='none', color='black')

        plt.ylabel('Density (N/arcmins$^2$)')
        plt.xlabel('Semi-major axis (arcmins)')

        xdata_orig = xdata

        xdata = xdata[:len(xdata)-(crowding_num)]
        ydata = ydata[:len(ydata)-(crowding_num)]
        yerr = yerr[:len(yerr)-(crowding_num)]

        if self.galaxy == 'ngc205':

            model = Sersic1D(r_eff=2.0, n=1.5)

        elif self.galaxy == 'm32':

            model = Sersic1D(r_eff=1.7, n=1.00)

        fit = fitting.LevMarLSQFitter()

        s = fit(model, xdata, ydata, weights=1/yerr)

        xdata = np.linspace(min(xdata_orig), max(xdata_orig), num=1000)

        plt.plot(xdata, s(xdata), label='Sersic Fit with a$_{eff}$ = ' + str(
            round(s.r_eff.value, 3)) + ', n = ' + str(round(s.n.value, 3)))
        plt.legend()
